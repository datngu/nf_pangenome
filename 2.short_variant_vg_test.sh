#!/bin/bash
#SBATCH --job-name=deep_variant_test
#SBATCH --account=nn9114k
#SBATCH --time=8:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=deep_variant_test_%j.log

# Load Nextflow
module load Nextflow/24.04.2



## set Singularity cache and tmp directories
mkdir -p /cluster/projects/nn9114k/datngu/singularity/tmp

export SINGULARITY_CACHEDIR="/cluster/projects/nn9114k/datngu/singularity"
export SINGULARITY_TMPDIR="/cluster/projects/nn9114k/datngu/singularity/tmp"

## set Nextflow Singularity cache directory
export NXF_SINGULARITY_CACHEDIR="/cluster/projects/nn9114k/datngu/singularity"

## singularity build quay.io-vgteam-vg-v1.65.0.img docker://quay.io/vgteam/vg:v1.65.0


nextflow run main_short_var.nf \
    --input_csv samples.csv \
    --reference test_data/GRCh38.fa \
    --graph test_data/hprc-v1.1-mc-grch38.gbz \
    --outdir test_results_deep_variants \
    -profile saga,singularity \
    -resume

echo "✓ Pipeline finished!"
