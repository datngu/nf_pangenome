#!/bin/bash
set -euo pipefail

echo "=============================================="
echo "  Preparing Pipeline Test Run"
echo "=============================================="

# Check if test data exists
if [ ! -d "test_data" ]; then
    echo "Error: test_data directory not found"
    echo "Please run: bash bin/download_test_data.sh first"
    exit 1
fi

# Create samples CSV with NA12878 test data
echo ""
echo "Creating samples.csv..."
cat > samples.csv << 'EOF'
sample,read1,read2
NA12878_200M,test_data/A006850052_NA12878_200M_R1.fq.gz,test_data/A006850052_NA12878_200M_R2.fq.gz
NA12878_75M,test_data/A006850052_NA12878_75M_R1.fq.gz,test_data/A006850052_NA12878_75M_R2.fq.gz
EOF

echo "✓ samples.csv created with NA12878 test samples"
cat samples.csv

# Create Singularity cache directory
echo ""
echo "Creating Singularity cache directory..."
mkdir -p singularity_cache
echo "✓ singularity_cache/ created"

# Set environment variables
echo ""
echo "Setting environment variables..."
export SINGULARITY_CACHEDIR="$(pwd)/singularity_cache"
echo "export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache" > run_env.sh
echo "✓ run_env.sh created"

echo ""
echo "=============================================="
echo "  Setup Complete!"
echo "=============================================="
echo ""
echo "To run the pipeline on SAGA:"
echo ""
echo "  # Load environment"
echo "  source run_env.sh"
echo "  module load Nextflow/23.04.0"
echo ""
echo "  # Run pipeline"
echo "  nextflow run main.nf \\"
echo "    --input_csv samples.csv \\"
echo "    --reference test_data/GRCh38_chr22.fa \\"
echo "    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \\"
echo "    --outdir results \\"
echo "    -profile saga"
echo ""
echo "To run the pipeline on TSD:"
echo ""
echo "  # Load environment"
echo "  source run_env.sh"
echo "  module load Nextflow/23.04.0"
echo ""
echo "  # Run pipeline"
echo "  nextflow run main.nf \\"
echo "    --input_csv samples.csv \\"
echo "    --reference test_data/GRCh38_chr22.fa \\"
echo "    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \\"
echo "    --outdir results \\"
echo "    -profile tsd"
echo ""
