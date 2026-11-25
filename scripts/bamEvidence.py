#!/usr/bin/env python3
"""Count the reads of one BAM at every SV breakend of the cohort.

Prism recounts the evidence for each breakend in the BAM files, without
the SV callers: the depth and mapping quality at the breakend, the reads
that are clipped there, the clipped reads whose other part aligns at the
mate breakend, and the read pairs that span the two breakends. Every BAM
of the cohort (each tumor and each normal) is counted at every breakend
of every group, so svCohort.R can report the counts in the group's own
tumor, in its own normal, and in the other samples.

Breakends: every record of every group's PURPLE SV VCF
(OUT/<group>/purple/<tumor>.purple.sv.vcf.gz; all FILTER values). The
mate of a breakend is the record its MATEID names; single and inferred
breakends have none.

Orientation (as ESVEE and LINX): 1 if ALT starts with the reference base,
else -1. A breakend with orientation 1 keeps the reference up to POS (the
supporting reads align up to POS and continue elsewhere); orientation -1
keeps the reference from POS on.

Reads: unmapped, secondary, duplicate and QC-fail records are skipped.
Primary records give the depth and the read pairs; primary and
supplementary records give the clipped reads. Counts are of distinct
read names (fragments). No mapping quality threshold is applied to the
counts; the mean mapping quality of the counted reads is reported
instead.

The counts are a recount from the alignments, not a copy of ESVEE's:
ESVEE assembles the reads, so it also counts reads that extend into an
inserted sequence or a repeat without a clip at POS. clip_n and split_n
are therefore often lower than PURPLE's SF.

Columns of the output (one row per breakend; NA where not measured):
    group_id, id        the breakend (group and VCF ID)
    depth               primary reads aligned over POS
    mapq_mean           their mean mapping quality
    mapq0_frac          the fraction of them with mapping quality 0
    softclip_frac_w50   among primary reads that overlap POS +/- 50,
                        the fraction with a soft clip of at least
                        MIN_CLIP bases at either end; the 3' clip of a
                        read in a proper pair whose fragment (|TLEN|) is
                        no longer than the read (adapter read-through)
                        is not counted
    clip_n              reads clipped (soft or hard, at least MIN_CLIP
                        bases) on the breakend's side, with the clip at
                        POS + CIPOS, widened by CLIP_TOL bases
    split_n             of those, the reads with an SA alignment at the
                        mate breakend (same rule at the mate's POS,
                        CIPOS and orientation). NA without a mate
    disc_n              read pairs with one read on the breakend's side
                        within DISC_WIN bases and strand, and the mate
                        on the mate breakend's side within DISC_WIN and
                        strand. NA without a mate, and NA when the two
                        breakends are on the same chromosome less than
                        DISC_MIN_SPAN apart (normal pairs fit that
                        geometry)
    improper_frac_w1000 among primary paired reads on the breakend's side
                        within DISC_WIN bases and strand, the fraction
                        whose pair the aligner did not flag proper
                        (mate unmapped, on another chromosome, wrong
                        orientation, or insert size outside the
                        library's distribution). A measure of the local
                        background: in these libraries many pairs are
                        chimeric (mates on random chromosomes), so a
                        count of improper pairs is not support
    support_mapq_mean   mean mapping quality of the reads counted in
                        clip_n and disc_n

Usage:
    python3 bamEvidence.py --out-dir OUT --samples SAMPLES.tsv --bam BAM \\
        --sample-id SAMPLE --output OUT.tsv.gz [--threads N]

SAMPLES.tsv is res/prism/db/samples.tsv (postOA.R): its group_id and
sample_id columns name the groups and their PURPLE VCFs.

The run stops if the BAM's read groups do not all name SAMPLE, if a
PURPLE VCF is missing, if an ALT has no recognizable orientation, or if a
MATEID names no record.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import pysam

MIN_CLIP = 5
CLIP_TOL = 5
SOFTCLIP_WIN = 50
DISC_WIN = 1000
DISC_MIN_SPAN = 1000

COLUMNS = [
    "group_id", "id", "depth", "mapq_mean", "mapq0_frac", "softclip_frac_w50",
    "clip_n", "split_n", "disc_n", "improper_frac_w1000", "support_mapq_mean",
]

# CIGAR operations: 0 M, 1 I, 2 D, 3 N, 4 S, 5 H, 6 P, 7 =, 8 X.
REF_OPS = {"M", "D", "N", "=", "X"}
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass
class Breakend:
    """One breakend of a PURPLE SV VCF record.

    Attributes:
        group_id: oncoanalyser group.
        id: VCF ID.
        chrom: VCF CHROM.
        pos: VCF POS (1-based).
        orient: 1 or -1 (see the module docstring).
        ci_lower: first CIPOS value (0 when absent).
        ci_upper: second CIPOS value (0 when absent).
        mate_id: VCF ID of the mate, or None.
        mate: the mate Breakend, or None.
    """

    group_id: str
    id: str
    chrom: str
    pos: int
    orient: int
    ci_lower: int
    ci_upper: int
    mate_id: str | None
    mate: Breakend | None = None


def orientation(ref: str, alt: str) -> int:
    """Return the breakend orientation from REF and ALT.

    Args:
        ref: VCF REF.
        alt: VCF ALT.

    Returns:
        1 if ALT starts with the reference base, else -1.

    Raises:
        ValueError: if ALT neither starts nor ends with the reference base.
    """
    if alt.startswith(ref):
        return 1
    if alt.endswith(ref):
        return -1
    raise ValueError(f"no orientation in ALT '{alt}' (REF '{ref}')")


def read_groups(samples_file: Path) -> list[tuple[str, str]]:
    """Return (group_id, tumor sample_id) of every group in the run.

    Args:
        samples_file: TSV with the columns group_id and sample_id.

    Returns:
        The groups in file order.
    """
    with open(samples_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        gi, si = header.index("group_id"), header.index("sample_id")
        return [(r[gi], r[si]) for r in (line.rstrip("\n").split("\t") for line in f)]


def read_breakends(out_dir: Path, groups: list[tuple[str, str]]) -> list[Breakend]:
    """Read every record of every group's PURPLE SV VCF.

    Args:
        out_dir: oncoanalyser output folder.
        groups: (group_id, tumor sample_id) pairs.

    Returns:
        The breakends, with mates linked.

    Raises:
        FileNotFoundError: if a VCF is missing.
        ValueError: on an ALT without orientation or a MATEID without record.
    """
    breakends = []
    for group_id, sample_id in groups:
        vcf = out_dir / group_id / "purple" / f"{sample_id}.purple.sv.vcf.gz"
        if not vcf.exists():
            raise FileNotFoundError(f"PURPLE SV VCF not found: {vcf}")
        by_id = {}
        with pysam.VariantFile(str(vcf)) as f:
            for rec in f:
                if rec.id is None or rec.ref is None or rec.alts is None or len(rec.alts) != 1:
                    raise ValueError(f"{vcf}: record at {rec.chrom}:{rec.pos} lacks one ID, REF or ALT")
                cipos = rec.info.get("CIPOS", (0, 0))
                mate_id = rec.info.get("MATEID")
                b = Breakend(
                    group_id=group_id, id=rec.id, chrom=rec.chrom, pos=rec.pos,
                    orient=orientation(rec.ref, rec.alts[0]),
                    ci_lower=int(cipos[0]), ci_upper=int(cipos[1]),
                    mate_id=mate_id if isinstance(mate_id, str) else None,
                )
                if b.id in by_id:
                    raise ValueError(f"{vcf}: ID {b.id} appears twice")
                by_id[b.id] = b
                breakends.append(b)
        for b in by_id.values():
            if b.mate_id is not None:
                if b.mate_id not in by_id:
                    raise ValueError(f"{vcf}: MATEID {b.mate_id} of {b.id} names no record")
                b.mate = by_id[b.mate_id]
    return breakends


def ref_length(cigar: str) -> int:
    """Return the number of reference bases a CIGAR string covers.

    Args:
        cigar: CIGAR string, e.g. "60S90M".

    Returns:
        The sum of the M, D, N, = and X lengths.
    """
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in REF_OPS)


def clip_lengths(read: pysam.AlignedSegment) -> tuple[int, int]:
    """Return the clipped bases (soft plus hard) at the left and right end.

    Args:
        read: an aligned read.

    Returns:
        (left clip, right clip) in bases.
    """
    cig = read.cigartuples
    left = right = 0
    for op, n in cig:
        if op in (4, 5):
            left += n
        else:
            break
    for op, n in reversed(cig):
        if op in (4, 5):
            right += n
        else:
            break
    return left, right


def adapter_end(read: pysam.AlignedSegment) -> str | None:
    """Return the end of the read that runs into the adapter, if any.

    In a proper pair whose fragment (|TLEN|) is no longer than the read,
    the read runs past the fragment end into the adapter, and the aligner
    clips its 3' end. That end is "right" for a forward read and "left"
    for a reverse read.

    Args:
        read: an aligned primary read.

    Returns:
        "left", "right" or None.
    """
    if read.is_proper_pair and 0 < abs(read.template_length) <= read.query_length:
        return "left" if read.is_reverse else "right"
    return None


def soft_clipped(read: pysam.AlignedSegment) -> bool:
    """Return True if an end has a soft clip of at least MIN_CLIP bases.

    The 3' clip of a read that runs into the adapter (adapter_end) is not
    counted.

    Args:
        read: an aligned read.

    Returns:
        Whether the read is soft-clipped.
    """
    cig = read.cigartuples
    adapter = adapter_end(read)
    left = cig[0][0] == 4 and cig[0][1] >= MIN_CLIP and adapter != "left"
    right = cig[-1][0] == 4 and cig[-1][1] >= MIN_CLIP and adapter != "right"
    return left or right


def junction_window(b: Breakend) -> tuple[int, int]:
    """Return the 1-based positions where a clip supports the breakend.

    Args:
        b: the breakend.

    Returns:
        (first, last) position, POS + CIPOS widened by CLIP_TOL.
    """
    return b.pos + b.ci_lower - CLIP_TOL, b.pos + b.ci_upper + CLIP_TOL


def clip_at(read: pysam.AlignedSegment, b: Breakend) -> bool:
    """Return True if the read is clipped at the breakend, on its side.

    Orientation 1: the alignment ends at the breakend and the right end is
    clipped. Orientation -1: the alignment starts there and the left end is
    clipped.

    Args:
        read: an aligned read.
        b: the breakend.

    Returns:
        Whether the read supports the breakend by its clip.
    """
    lo, hi = junction_window(b)
    left, right = clip_lengths(read)
    if b.orient == 1:
        return right >= MIN_CLIP and lo <= read.reference_end <= hi
    return left >= MIN_CLIP and lo <= read.reference_start + 1 <= hi


def sa_at_mate(read: pysam.AlignedSegment, mate: Breakend) -> bool:
    """Return True if an SA alignment of the read ends at the mate breakend.

    Orientation 1 at the mate: the SA alignment's last reference base is in
    the mate's junction window; orientation -1: its first base is.

    Args:
        read: an aligned read.
        mate: the mate breakend.

    Returns:
        Whether a supplementary part of the read aligns at the mate.
    """
    if not read.has_tag("SA"):
        return False
    lo, hi = junction_window(mate)
    for entry in read.get_tag("SA").rstrip(";").split(";"):
        rname, pos, _strand, cigar, _mapq, _nm = entry.split(",")
        if rname != mate.chrom:
            continue
        start = int(pos)
        end = start + ref_length(cigar) - 1
        edge = end if mate.orient == 1 else start
        if lo <= edge <= hi:
            return True
    return False


def on_side(read: pysam.AlignedSegment, b: Breakend) -> bool:
    """Return True if a primary read lies on the breakend's side and strand.

    Orientation 1: forward strand, alignment start in [POS - DISC_WIN, POS].
    Orientation -1: reverse strand, alignment end in [POS, POS + DISC_WIN].

    Args:
        read: an aligned read.
        b: the breakend.

    Returns:
        Whether the read can belong to a pair spanning the breakend.
    """
    if b.orient == 1:
        return not read.is_reverse and b.pos - DISC_WIN <= read.reference_start + 1 <= b.pos
    return read.is_reverse and b.pos <= read.reference_end <= b.pos + DISC_WIN


def mate_on_side(read: pysam.AlignedSegment, mate: Breakend) -> bool:
    """Return True if the read's mate lies on the mate breakend's side.

    The mate's end comes from its MC tag (mate CIGAR); without MC its
    start plus the read length is used.

    Args:
        read: an aligned primary read with a mapped mate.
        mate: the mate breakend.

    Returns:
        Whether the mate fits the mate breakend's side and strand.
    """
    if read.next_reference_name != mate.chrom:
        return False
    start = read.next_reference_start + 1
    if mate.orient == 1:
        return not read.mate_is_reverse and mate.pos - DISC_WIN <= start <= mate.pos
    if not read.mate_is_reverse:
        return False
    if read.has_tag("MC"):
        end = start + ref_length(read.get_tag("MC")) - 1
    else:
        end = start + read.query_length - 1
    return mate.pos <= end <= mate.pos + DISC_WIN


def normal_pair(read: pysam.AlignedSegment) -> bool:
    """Return True if a primary read and its mate form a normal pair.

    Normal: the aligner flagged the pair proper (bwa-mem2: FR orientation
    and an insert size within the library's distribution).

    Args:
        read: an aligned primary paired read.

    Returns:
        Whether the pair is normal.
    """
    return read.is_proper_pair


def mean(values: list[int]) -> float | None:
    """Return the mean of a list, or None if it is empty.

    Args:
        values: numbers.

    Returns:
        The mean or None.
    """
    return sum(values) / len(values) if values else None


def count(bam: pysam.AlignmentFile, b: Breakend) -> dict[str, object]:
    """Count the reads of the BAM at one breakend.

    Args:
        bam: the open BAM.
        b: the breakend.

    Returns:
        The output columns for the breakend (None for NA).
    """
    mate = b.mate
    disc_measured = mate is not None and not (
        mate.chrom == b.chrom and abs(mate.pos - b.pos) < DISC_MIN_SPAN
    )
    pad = max(DISC_WIN, SOFTCLIP_WIN, CLIP_TOL - b.ci_lower, b.ci_upper + CLIP_TOL) + 1
    depth_mapq = []
    win_reads = win_clipped = 0
    side_reads = side_improper = 0
    clip, split, disc = set(), set(), set()
    support_mapq = {}
    if b.chrom not in bam.references:
        raise ValueError(f"{b.group_id} {b.id}: chromosome {b.chrom} not in the BAM")
    for read in bam.fetch(b.chrom, max(0, b.pos - 1 - pad), b.pos + pad):
        if read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail:
            continue
        name = read.query_name
        primary = not read.is_supplementary
        start1, end1 = read.reference_start + 1, read.reference_end
        if primary and start1 <= b.pos <= end1:
            depth_mapq.append(read.mapping_quality)
        if primary and start1 <= b.pos + SOFTCLIP_WIN and end1 >= b.pos - SOFTCLIP_WIN:
            win_reads += 1
            win_clipped += soft_clipped(read)
        if clip_at(read, b):
            clip.add(name)
            support_mapq.setdefault(name, read.mapping_quality)
            if mate is not None and sa_at_mate(read, mate):
                split.add(name)
        if primary and read.is_paired and on_side(read, b):
            side_reads += 1
            side_improper += not normal_pair(read)
            if disc_measured and not read.mate_is_unmapped and mate_on_side(read, mate):
                disc.add(name)
                support_mapq.setdefault(name, read.mapping_quality)
    return {
        "group_id": b.group_id,
        "id": b.id,
        "depth": len(depth_mapq),
        "mapq_mean": mean(depth_mapq),
        "mapq0_frac": (sum(q == 0 for q in depth_mapq) / len(depth_mapq)) if depth_mapq else None,
        "softclip_frac_w50": (win_clipped / win_reads) if win_reads else None,
        "clip_n": len(clip),
        "split_n": len(split) if mate is not None else None,
        "disc_n": len(disc) if disc_measured else None,
        "improper_frac_w1000": (side_improper / side_reads) if side_reads else None,
        "support_mapq_mean": mean(list(support_mapq.values())),
    }


def check_sample(bam: pysam.AlignmentFile, sample_id: str) -> None:
    """Stop unless every read group of the BAM names the sample.

    Args:
        bam: the open BAM.
        sample_id: the expected SM.

    Raises:
        ValueError: if a read group names another sample, or there is none.
    """
    names = {rg.get("SM") for rg in bam.header.to_dict().get("RG", [])}
    if names != {sample_id}:
        raise ValueError(f"BAM read groups name {sorted(map(str, names))}, expected {sample_id}")


def fmt(value: object) -> str:
    """Format a value for the TSV: NA for None, 6 significant digits for floats.

    Args:
        value: a column value.

    Returns:
        The text.
    """
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def main() -> int:
    """Parse the arguments, count every breakend and write the table.

    Returns:
        The exit status.
    """
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--out-dir", required=True, type=Path, help="oncoanalyser output folder")
    p.add_argument("--samples", required=True, type=Path, help="samples.tsv (group_id, sample_id)")
    p.add_argument("--bam", required=True, type=Path, help="BAM to count")
    p.add_argument("--sample-id", required=True, help="sample the BAM holds (its SM)")
    p.add_argument("--output", required=True, type=Path, help="output TSV (.gz)")
    p.add_argument("--limit", type=int, default=None, help="count only the first N breakends (tests)")
    p.add_argument("--threads", type=int, default=1, help="BAM decompression threads")
    a = p.parse_args()

    breakends = read_breakends(a.out_dir, read_groups(a.samples))
    print(f"{len(breakends)} breakends in {len({b.group_id for b in breakends})} groups", flush=True)
    if a.limit is not None:
        breakends = breakends[: a.limit]
    part = a.output.with_name(a.output.name + ".part")
    with pysam.AlignmentFile(str(a.bam), "rb", threads=a.threads) as bam:
        check_sample(bam, a.sample_id)
        with gzip.open(part, "wt") as out:
            out.write("\t".join(COLUMNS) + "\n")
            for i, b in enumerate(breakends, 1):
                row = count(bam, b)
                out.write("\t".join(fmt(row[c]) for c in COLUMNS) + "\n")
                if i % 5000 == 0:
                    print(f"  {i} breakends", flush=True)
    part.rename(a.output)
    print(f"{len(breakends)} breakends counted in {a.bam} -> {a.output}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
