#!/bin/bash
# Quick test script for SAGA

set -euo pipefail

# Source environment
source run_env.sh

# Load Nextflow
module load Nextflow/23.04.0

# Run pipeline
nextflow run main.nf \
    --input_csv samples.csv \
    --reference test_data/GRCh38_chr22.fa \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --outdir results_test \
    -profile saga \
    -resume

echo ""
echo "Pipeline complete! Check results in: results_test/"
