#!/bin/bash
set -euo pipefail

echo "=============================================="
echo "  Downloading Test Data"
echo "=============================================="

# Create test data directory
mkdir -p test_data
cd test_data

echo ""
echo "1. Downloading reference genome (chr22 only for testing)..."
if [ ! -f "GRCh38_chr22.fa" ]; then
    # Download full reference
    wget -O GRCh38_full.fa.gz \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    
    echo "Extracting chr22..."
    gunzip GRCh38_full.fa.gz
    
    # Extract chr22 only (requires samtools)
    if command -v samtools &> /dev/null; then
        samtools faidx GRCh38_full.fa chr22 > GRCh38_chr22.fa
        samtools faidx GRCh38_chr22.fa
        rm GRCh38_full.fa GRCh38_full.fa.fai
    else
        echo "Warning: samtools not found. Keeping full reference."
        mv GRCh38_full.fa GRCh38_chr22.fa
    fi
    
    echo "✓ Reference genome ready"
else
    echo "✓ Reference genome already exists"
fi

echo ""
echo "2. Downloading HPRC v1.1 pangenome graph..."
if [ ! -f "hprc-v1.1-mc-grch38.gbz" ]; then
    echo "Downloading HPRC v1.1 graph (~50GB)..."
    echo "This may take a while..."
    
    wget --progress=bar:force \
        -O hprc-v1.1-mc-grch38.gbz \
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz"
    
    echo "✓ HPRC v1.1 graph downloaded"
else
    echo "✓ HPRC v1.1 graph already exists"
fi

echo ""
echo "3. Downloading NA12878 test reads from Zenodo..."

if [ ! -f "A006850052_NA12878_200M_R1.fq.gz" ]; then
    echo "Downloading NA12878 200M coverage reads..."
    wget --progress=bar:force \
        -O A006850052_NA12878_200M_R1.fq.gz \
        "https://zenodo.org/records/6513789/files/A006850052_NA12878_200M_R1.fq.gz?download=1"
    
    wget --progress=bar:force \
        -O A006850052_NA12878_200M_R2.fq.gz \
        "https://zenodo.org/records/6513789/files/A006850052_NA12878_200M_R2.fq.gz?download=1"
    
    echo "✓ NA12878 200M reads downloaded"
else
    echo "✓ NA12878 200M reads already exist"
fi

if [ ! -f "A006850052_NA12878_75M_R1.fq.gz" ]; then
    echo "Downloading NA12878 75M coverage reads..."
    wget --progress=bar:force \
        -O A006850052_NA12878_75M_R1.fq.gz \
        "https://zenodo.org/records/6513789/files/A006850052_NA12878_75M_R1.fq.gz?download=1"
    
    wget --progress=bar:force \
        -O A006850052_NA12878_75M_R2.fq.gz \
        "https://zenodo.org/records/6513789/files/A006850052_NA12878_75M_R2.fq.gz?download=1"
    
    echo "✓ NA12878 75M reads downloaded"
else
    echo "✓ NA12878 75M reads already exist"
fi

echo ""
echo "=============================================="
echo "  Download Complete!"
echo "=============================================="
echo ""
echo "Test data location: $(pwd)"
echo ""
echo "Files:"
ls -lh GRCh38_chr22.fa* 2>/dev/null || echo "  - Reference: Not found"
ls -lh hprc-v1.1-mc-grch38.gbz 2>/dev/null || echo "  - HPRC v1.1 graph: Not found"
ls -lh A006850052_NA12878_*_R*.fq.gz 2>/dev/null || echo "  - NA12878 test reads: Not found"
echo ""
