#!/bin/bash
#SBATCH --job-name=nf_pangenome_test
#SBATCH --account=nn9114k
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=nf_pangenome_test_%j.log

# SLURM sbatch script for running pangenome pipeline test on SAGA
# Usage: sbatch run_test_saga.sh

set -o errexit  # Exit on errors
set -o nounset  # Treat unset variables as errors

# Load required modules
module --quiet purge
module load Nextflow/24.04.2

# Set Singularity cache directory (avoid filling $HOME)
export SINGULARITY_CACHEDIR=${SLURM_SUBMIT_DIR}/singularity_cache
export NXF_SINGULARITY_CACHEDIR=${SLURM_SUBMIT_DIR}/singularity_cache

echo "=========================================="
echo "Starting Pangenome Pipeline Test on SAGA"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $SLURM_SUBMIT_DIR"
echo "Singularity cache: $SINGULARITY_CACHEDIR"
echo "=========================================="
echo ""

# Change to submission directory
cd ${SLURM_SUBMIT_DIR}

# Verify test data exists
if [ ! -d "test_data" ]; then
    echo "ERROR: test_data directory not found!"
    echo "Please run the download steps first (see README.md)"
    exit 1
fi

if [ ! -f "samples.csv" ]; then
    echo "ERROR: samples.csv not found!"
    echo "Please run the preparation steps first (see README.md)"
    exit 1
fi

# Run Nextflow pipeline
echo "Running Nextflow pipeline with full genome..."
echo ""

nextflow run main.nf \
    --input_csv samples.csv \
    --reference test_data/GRCh38.fa \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --outdir results_test \
    -profile saga \
    -resume

echo ""
echo "=========================================="
echo "Pipeline Test Complete!"
echo "=========================================="
echo "Results directory: results_test/"
echo "Check log file: nf_pangenome_test_${SLURM_JOB_ID}.log"
echo "=========================================="
