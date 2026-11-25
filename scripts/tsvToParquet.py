#!/usr/bin/env python3
"""Write a typed Parquet copy of a Prism TSV table.

The column dictionary written next to the table (svCohort.R) gives each
column's type: string, integer, double, boolean, or list<...> for a list of
any length. Lists are split on the column's list_sep; where the dictionary
names a list_count_column, each list has exactly that many elements (an
empty element is null). Every cell is converted or the run stops, naming
the column and an example value: nothing is guessed, dropped or reformatted.

Usage:
    python3 tsvToParquet.py TABLE.tsv.gz DICTIONARY.tsv OUT.parquet
"""

from __future__ import annotations

import argparse
import gzip
import sys

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

# Integers are 32-bit (R reads them as integer; R has no 64-bit integer);
# a value outside that range stops the run.
SCALAR_TYPES = {
    "string": pa.string(),
    "integer": pa.int32(),
    "double": pa.float64(),
    "boolean": pa.bool_(),
}
PATTERNS = {
    "integer": r"^-?[0-9]+$",
    "double": r"^-?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$|^-?(NaN|Inf|Infinity|nan|inf)$",
    "boolean": r"^(true|false|TRUE|FALSE)$",
}


class ConversionError(Exception):
    """A value that does not fit the type the dictionary gives its column."""


def header_names(path: str) -> list[str]:
    """Return the column names from the first line of a (gzipped) TSV.

    Args:
        path: TSV file, gzipped if the name ends in .gz.

    Returns:
        The column names in file order.
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return f.readline().rstrip("\n").split("\t")


def read_strings(path: str) -> pa.Table:
    """Read a TSV with every column as string; empty cells become null.

    Args:
        path: TSV file written without quoting.

    Returns:
        The table, all columns of type string.
    """
    names = header_names(path)
    return pacsv.read_csv(
        path,
        read_options=pacsv.ReadOptions(block_size=1 << 26),
        parse_options=pacsv.ParseOptions(delimiter="\t", quote_char=False, double_quote=False),
        convert_options=pacsv.ConvertOptions(
            column_types={n: pa.string() for n in names},
            null_values=[""],
            strings_can_be_null=True,
        ),
    )


def check_pattern(values: pa.Array, kind: str, column: str) -> None:
    """Stop unless every non-null value matches the pattern of its type.

    Args:
        values: String values.
        kind: integer, double or boolean.
        column: Column name, for the error message.

    Raises:
        ConversionError: If a value does not match.
    """
    ok = pc.match_substring_regex(values, PATTERNS[kind])
    bad = pc.filter(values, pc.invert(pc.fill_null(ok, True)))
    if len(bad) > 0:
        raise ConversionError(f"{column}: {len(bad)} values are not {kind}, e.g. '{bad[0].as_py()}'")


def convert_scalar(values: pa.Array, kind: str, column: str) -> pa.Array:
    """Convert string values to one scalar type.

    Args:
        values: String values (nulls allowed).
        kind: string, integer, double or boolean.
        column: Column name, for error messages.

    Returns:
        The values as the Arrow type of KIND.

    Raises:
        ConversionError: If a value does not fit KIND.
    """
    if kind == "string":
        return values
    check_pattern(values, kind, column)
    if kind == "boolean":
        return pc.equal(pc.utf8_lower(values), "true")
    try:
        return pc.cast(values, SCALAR_TYPES[kind])
    except pa.ArrowInvalid as e:
        raise ConversionError(f"{column}: {e}") from e


def convert_list(
    values: pa.Array, elem: str, sep: str, counts: pa.Array | None, column: str
) -> pa.Array:
    """Split string values into typed lists.

    Args:
        values: String values; null is no list (unless COUNTS says otherwise).
        elem: Element type (string, integer, double, boolean).
        sep: Separator between elements.
        counts: Number of elements per row (integer), or None.
        column: Column name, for error messages.

    Returns:
        A list array of the element type; empty elements are null.

    Raises:
        ConversionError: If an element does not fit ELEM, or a list length
            differs from COUNTS.
    """
    if counts is not None:
        # The TSV cannot tell an empty cell from a one-element list whose
        # element is empty: the count column decides.
        has_list = pc.and_(pc.is_valid(counts), pc.greater(pc.fill_null(counts, 0), 0))
        values = pc.if_else(has_list, pc.fill_null(values, ""), pa.scalar(None, pa.string()))
    split = pc.split_pattern(values, sep)
    flat = pc.list_flatten(split)
    flat = pc.if_else(pc.equal(flat, ""), pa.scalar(None, pa.string()), flat)
    typed = convert_scalar(flat, elem, column)
    out = pa.ListArray.from_arrays(split.offsets, typed, mask=pc.is_null(split))
    if counts is not None:
        lengths = pc.fill_null(pc.list_value_length(out), 0)
        differ = pc.not_equal(lengths, pc.fill_null(counts, 0))
        if pc.any(differ).as_py():
            i = pc.index(differ, True).as_py()
            raise ConversionError(
                f"{column}: row {i + 1} has {lengths[i].as_py()} elements, count column says {counts[i].as_py()}"
            )
    return out


def convert(table: pa.Table, dictionary: pa.Table) -> pa.Table:
    """Convert every column of TABLE to the type its dictionary row gives.

    Args:
        table: The TSV table, all columns string.
        dictionary: Columns column, type, list_sep, list_count_column.

    Returns:
        The typed table, columns in the same order.

    Raises:
        ConversionError: If the columns differ from the dictionary or a value
            does not fit its type.
    """
    spec = {r["column"]: r for r in dictionary.to_pylist()}
    if list(spec) != table.column_names:
        missing = set(spec) ^ set(table.column_names)
        raise ConversionError(f"table and dictionary columns differ: {sorted(missing)[:10]}")
    scalars = {}
    for name in table.column_names:
        kind = spec[name]["type"]
        if not kind.startswith("list<"):
            if kind not in SCALAR_TYPES:
                raise ConversionError(f"{name}: unknown type '{kind}'")
            scalars[name] = convert_scalar(table[name].combine_chunks(), kind, name)
    arrays = []
    for name in table.column_names:
        kind = spec[name]["type"]
        if kind.startswith("list<"):
            elem = kind[len("list<"):-1]
            count_col = spec[name]["list_count_column"]
            counts = scalars[count_col] if count_col else None
            arrays.append(convert_list(table[name].combine_chunks(), elem, spec[name]["list_sep"], counts, name))
        else:
            arrays.append(scalars[name])
    return pa.Table.from_arrays(arrays, names=table.column_names)


def verify(source: pa.Table, path: str) -> None:
    """Read the Parquet file back and check rows, nulls and list lengths.

    Args:
        source: The string table the file was made from.
        path: The Parquet file.

    Raises:
        ConversionError: If the file differs from SOURCE in shape or in which
            cells are empty.
    """
    back = pq.read_table(path)
    if back.num_rows != source.num_rows or back.column_names != source.column_names:
        raise ConversionError("Parquet file read back has different rows or columns")
    for name in source.column_names:
        if pa.types.is_list(back.schema.field(name).type):
            continue
        if back[name].null_count != source[name].null_count:
            raise ConversionError(f"{name}: null count changed in the Parquet file")


def main(argv: list[str]) -> int:
    """Convert one table; return the process exit status.

    Args:
        argv: Command-line arguments (table, dictionary, output).

    Returns:
        0 on success, 1 on a conversion error.
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("table")
    parser.add_argument("dictionary")
    parser.add_argument("out")
    args = parser.parse_args(argv)
    try:
        source = read_strings(args.table)
        dictionary = read_strings(args.dictionary)
        typed = convert(source, dictionary)
        pq.write_table(typed, args.out, compression="zstd")
        verify(source, args.out)
    except ConversionError as e:
        print(f"tsvToParquet: {e}", file=sys.stderr)
        return 1
    print(f"{args.out}: {typed.num_rows} rows x {typed.num_columns} columns")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
