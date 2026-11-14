export NXF_VER=25.10.4
export NXF_HOME=/scratch/core001/bic/socci/Oncoanalyser/nextflow_home
# New work folder each time this file is sourced (format 20261112_193000).
# To -resume a run, set NXF_WORK back to that run's folder after sourcing.
export NXF_WORK=/scratch/core001/bic/socci/Oncoanalyser/$(date +%Y%m%d_%H%M%S)
export NXF_SINGULARITY_CACHEDIR=/scratch/core001/bic/socci/singularity-cache
export SINGULARITY_CACHEDIR=/scratch/core001/bic/socci/singularity-cache
export APPTAINER_CACHEDIR=/scratch/core001/bic/socci/singularity-cache
export NXF_OPTS='-Xms1g -Xmx4g'
