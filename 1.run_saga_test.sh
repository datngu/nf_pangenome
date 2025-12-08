#!/bin/bash
#SBATCH --job-name=pangenie_test
#SBATCH --account=nn9114k
#SBATCH --time=8:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=pangenie_test_%j.log

# Load Nextflow
module load Nextflow/24.04.2



## set Singularity cache and tmp directories

export SINGULARITY_CACHEDIR="/cluster/projects/nn9114k/datngu/singularity"
export SINGULARITY_TMPDIR="/cluster/projects/nn9114k/datngu/singularity/tmp"

## singularity build vg_v1.65.0.sif docker://quay.io/vgteam/vg:v1.65.0




nextflow run main.nf \
    --input_csv samples.csv \
    --reference reference.fa \
    --phased_vcf data/1000GP_ONT_shapeit5-phased-callset_final-vcf.phased.vcf.gz \
    --outdir results_pangenie \
    -profile saga \
    -resume

echo "Pipeline completed!"
