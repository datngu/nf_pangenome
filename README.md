## Quick Start - Test Pipeline

Follow these steps to download test data and run the pipeline on SAGA:

### Step 1: Clone Repository

```bash
git clone https://github.com/datngu/nf_pangenome.git
cd nf_pangenome
```

### Step 2: Download Test Data (~55 GB total)

Copy and paste these commands to download reference genome, HPRC graph, and NA12878 test reads:

```bash
# Create test data directory
mkdir -p test_data
cd test_data

# 1. Download GRCh38 full reference genome
echo "Downloading GRCh38 reference genome..."
wget -O GRCh38.fa.gz \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

gunzip GRCh38.fa.gz

# Index reference (optional but recommended)
module load SAMtools/1.17-GCC-12.2.0
samtools faidx GRCh38.fa

echo "✓ GRCh38 reference genome ready"

# 2. Download HPRC v1.1 pangenome graph (~50GB - this takes time!)
echo "Downloading HPRC v1.1 graph (~50GB)..."
wget --progress=bar:force \
    -O hprc-v1.1-mc-grch38.gbz \
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz"

echo "✓ HPRC v1.1 graph downloaded"

# 3. Download NA12878 test reads from Zenodo
echo "Downloading NA12878 200M coverage reads..."
wget --progress=bar:force \
    -O A006850052_NA12878_200M_R1.fq.gz \
    "https://zenodo.org/records/6513789/files/A006850052_NA12878_200M_R1.fq.gz?download=1"

wget --progress=bar:force \
    -O A006850052_NA12878_200M_R2.fq.gz \
    "https://zenodo.org/records/6513789/files/A006850052_NA12878_200M_R2.fq.gz?download=1"

echo "Downloading NA12878 75M coverage reads..."
wget --progress=bar:force \
    -O A006850052_NA12878_75M_R1.fq.gz \
    "https://zenodo.org/records/6513789/files/A006850052_NA12878_75M_R1.fq.gz?download=1"

wget --progress=bar:force \
    -O A006850052_NA12878_75M_R2.fq.gz \
    "https://zenodo.org/records/6513789/files/A006850052_NA12878_75M_R2.fq.gz?download=1"

echo "✓ NA12878 reads downloaded"

# Return to main directory
cd ..

echo "=========================================="
echo "Download Complete!"
echo "=========================================="
ls -lh test_data/

```

### Step 3: Prepare Sample Sheet

```bash
# Create samples.csv with test data (use absolute paths)
cat > samples.csv << EOF
sample,read1,read2
NA12878_200M,${PWD}/test_data/A006850052_NA12878_200M_R1.fq.gz,${PWD}/test_data/A006850052_NA12878_200M_R2.fq.gz
NA12878_75M,${PWD}/test_data/A006850052_NA12878_75M_R1.fq.gz,${PWD}/test_data/A006850052_NA12878_75M_R2.fq.gz
EOF

echo "✓ samples.csv created with absolute paths"
cat samples.csv


echo "✓ Setup complete! Ready to submit job."
```