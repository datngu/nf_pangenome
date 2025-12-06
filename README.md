# Pangenome Variant Calling Pipeline

A state-of-the-art Nextflow pipeline for pangenome-based variant calling using HPRC v1.1 graphs, pangenome-aware DeepVariant, and PanGenie. Designed for SAGA and TSD HPC environments with Singularity containers.

## Overview

### Why Pangenome-Based Variant Calling?

Traditional variant calling pipelines align sequencing reads to a single linear reference genome (e.g., GRCh38), which introduces **reference bias** and limits our ability to detect variants in structurally complex regions. This pipeline represents a **paradigm shift** by using pangenome graphs that incorporate genetic diversity from multiple individuals, enabling more accurate and comprehensive variant discovery.

### State-of-the-Art Technology Stack

This pipeline integrates three cutting-edge tools, each representing the latest advances in genomic analysis:

1. **vg Giraffe** - Ultra-fast pangenome graph alignment
   - Latest generation of the vg toolkit for variation graph analysis
   - Dramatically faster than traditional linear aligners while maintaining high accuracy
   - Leverages haplotype information from 94 diverse human genomes in HPRC v1.1

2. **Pangenome-Aware DeepVariant** - Graph-aware deep learning variant caller
   - Revolutionary extension of Google's DeepVariant using deep neural networks
   - First variant caller specifically designed to work with pangenome graphs
   - Significantly improves accuracy in complex genomic regions compared to linear-reference approaches

3. **PanGenie** - Kmer-based pangenome genotyping
   - Specialized for structural variant (SV) genotyping using known haplotypes
   - Exploits known haplotype paths through the pangenome graph
   - Provides accurate genotypes for large and complex structural variants

### Key Advantages Over Linear-Reference Pipelines

✨ **Improved Accuracy**: Reduces reference bias by incorporating diverse human genetic variation  
✨ **Better SV Calling**: Dramatically improved detection and genotyping of structural variants  
✨ **Complex Region Resolution**: Enhanced performance in segmental duplications, HLA, and other challenging loci  
✨ **Population-Aware**: Leverages population-scale haplotype information for more informed variant calls  
✨ **State-of-the-Art Methods**: Combines the latest advances in graph genomics and machine learning

**Key Features:**
- Multi-sample processing via CSV input
- Pangenome graph alignment with vg giraffe
- Pangenome-aware variant calling with DeepVariant
- SV genotyping with PanGenie
- Automatic HPRC v1.1 graph download
- Singularity containers for reproducibility
- Optimized for SAGA and TSD HPC systems
- **Ready-to-submit SLURM script** (`run_test_saga.sh`)

## Repository Structure

```
nf_pangenome/
├── run_test_saga.sh         # ⭐ Ready-to-submit SLURM script
├── main.nf                   # Main workflow
├── nextflow.config           # Configuration (SAGA/TSD profiles)
├── modules/                  # Pipeline processes
│   ├── download_graph.nf
│   ├── augment_graph.nf
│   ├── index_graph.nf
│   ├── giraffe_align.nf
│   ├── project_bam.nf
│   ├── call_snps_indels.nf
│   ├── call_svs.nf
│   └── postprocess.nf
├── bin/                      # Helper scripts (optional)
│   ├── download_test_data.sh
│   ├── prepare_test.sh
│   ├── run_test_saga.sh
│   └── run_test_tsd.sh
├── README.md                 # This file
└── PIPELINE_SUMMARY.md       # Technical documentation
```

## Quick Test on SAGA (TL;DR)

If you just want to run a quick test:

```bash
# 1. Clone repository
git clone https://github.com/datngu/nf_pangenome.git
cd nf_pangenome

# 2. Download test data (see Step 2 below - copy/paste commands) - ~55 GB
# 3. Prepare samples.csv (see Step 3 below - copy/paste commands)

# 4. Submit job
sbatch run_test_saga.sh

# Done! Check results in results_test/
```

See detailed step-by-step instructions below.

## Pipeline Workflow

```
┌─────────────────────────────────────────────────────────────────────┐
│                        INPUT: samples.csv                            │
│           (sample, read1.fq.gz, read2.fq.gz) × N samples            │
└────────────────────────────┬────────────────────────────────────────┘
                             │
                             ▼
         ┌───────────────────────────────────────┐
         │  1. DOWNLOAD_GRAPH (if not provided)  │
         │     Download HPRC v1.1 graph          │
         │     (~50GB, cached)                   │
         └───────────────┬───────────────────────┘
                         │
                         ▼
         ┌───────────────────────────────────────┐
         │  2. AUGMENT_GRAPH (optional)          │
         │     Add population SVs (1KGP)         │
         │     to pangenome graph                │
         └───────────────┬───────────────────────┘
                         │
                         ▼
         ┌───────────────────────────────────────┐
         │  3. INDEX_GRAPH                       │
         │     Build distance & minimizer        │
         │     indexes for fast alignment        │
         └───────────────┬───────────────────────┘
                         │
         ┌───────────────┴───────────────┐
         │    Process each sample        │
         └───────────────┬───────────────┘
                         │
        ┌────────────────┴────────────────┐
        │ PER-SAMPLE PARALLEL PROCESSING  │
        └────────────────┬────────────────┘
                         │
         ┌───────────────┴───────────────────────┐
         │  4. GIRAFFE_ALIGN                     │
         │     Align reads to pangenome graph    │
         │     Output: sample.gam                │
         └───────────────┬───────────────────────┘
                         │
         ┌───────────────┴───────────────────────┐
         │  5. PROJECT_BAM                       │
         │     Project GAM → BAM (hg38 linear)   │
         │     Sort & index BAM                  │
         └───────────────┬───────────────────────┘
                         │
         ┌───────────────┴───────────────────────┐
         │                                       │
         │  ┌─────────────────────────────────┐ │
         │  │  6. CALL_SNPS_INDELS            │ │
         │  │     Pangenome-aware DeepVariant │ │
         │  │     Output: sample.dv.vcf.gz    │ │
         │  └─────────────┬───────────────────┘ │
         │                │                     │
         │  ┌─────────────┴───────────────────┐ │
         │  │  7. CALL_SVS                    │ │
         │  │     PanGenie SV genotyping      │ │
         │  │     Output: sample.pg.vcf.gz    │ │
         │  └─────────────┬───────────────────┘ │
         │                │                     │
         └────────────────┴─────────────────────┘
                          │
         ┌────────────────┴───────────────────────┐
         │  8. POSTPROCESS                        │
         │     Merge SNPs/indels + SVs            │
         │     Output: sample.merged.vcf.gz       │
         └────────────────┬───────────────────────┘
                          │
                          ▼
         ┌────────────────────────────────────────┐
         │         OUTPUT PER SAMPLE              │
         │  - sample.merged.vcf.gz (all variants) │
         │  - sample.snps_indels.vcf.gz           │
         │  - sample.svs.vcf.gz                   │
         │  - sample.bam (optional)               │
         └────────────────────────────────────────┘
```

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

# Create Singularity cache directory
mkdir -p singularity_cache

echo "✓ Setup complete! Ready to submit job."
```

**Note:** The sbatch script (`run_test_saga.sh`) handles all environment setup automatically - no manual exports needed!

### Step 4: Submit Test Job to SAGA

Simply submit the included sbatch script:

```bash
sbatch run_test_saga.sh
```

The sbatch script will automatically:
- Load Nextflow/24.04.2 module
- Set Singularity cache to `singularity_cache/` (avoids filling $HOME)
- Verify test data and samples.csv exist
- Run pipeline with full genome reference
- Output results to `results_test/`
- Create log file: `nf_pangenome_test_<jobid>.log`

**Monitor the job:**
```bash
squeue -u $USER
```

**Check log file:**
```bash
tail -f nf_pangenome_test_*.log
```

### Alternative: Run Interactively on SAGA

If you prefer to run interactively (not recommended for full runs):

```bash
# Request interactive session
salloc --account=nn9114k --time=24:00:00 --mem=8G --cpus-per-task=2

# Load modules
module load Nextflow/24.04.2

# Set cache
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache

# Run pipeline
nextflow run main.nf \
    --input_csv samples.csv \
    --reference test_data/GRCh38.fa \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --outdir results_test \
    -profile saga \
    -resume
```

### Step 5: Check Results

After the pipeline completes:

```bash
# Check output structure
ls -lh results_test/

# Check final variants
ls -lh results_test/final/

# View variant counts
zcat results_test/final/NA12878_200M.merged.vcf.gz | grep -v "^#" | wc -l
```

---

## Production Data Setup

> **Note:** If you've completed the test run above, you already have the HPRC graph and reference downloaded in `test_data/`. You can reuse these files for production runs!

### Step-by-Step: Preparing Your Own Data

#### 1. Prepare Your Sequencing Reads

Organize your FASTQ files:
```bash
data/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
├── sample2_R2.fastq.gz
└── sample3_R1.fastq.gz
    sample3_R2.fastq.gz
```

#### 2. Download Reference Genome

```bash
# Create reference directory
mkdir -p reference
cd reference

# Download GRCh38 (no ALT contigs)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Extract
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fa

# Index (optional, but recommended)
samtools faidx GRCh38.fa

cd ..
```

#### 3. Download HPRC v1.1 Graph (Optional)

The pipeline will auto-download if not provided, but you can download manually to save time:

```bash
# Create graph directory
mkdir -p graphs
cd graphs

# Download HPRC v1.1 graph (~50GB)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz

# Download pre-built indexes (optional, pipeline can build them)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.dist
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.min

cd ..
```

#### 4. Create Sample Sheet

Create `samples.csv` with your samples:

```bash
cat > samples.csv << EOF
sample,read1,read2
patient001,/full/path/to/data/sample1_R1.fastq.gz,/full/path/to/data/sample1_R2.fastq.gz
patient002,/full/path/to/data/sample2_R1.fastq.gz,/full/path/to/data/sample2_R2.fastq.gz
patient003,/full/path/to/data/sample3_R1.fastq.gz,/full/path/to/data/sample3_R2.fastq.gz
EOF
```

**Important:** Use absolute paths in the CSV file!

#### 5. Setup Singularity Cache

```bash
# Create cache directory in working dir (not $HOME!)
mkdir -p singularity_cache
```

**Note:** When running pipeline, you'll need to set the cache location:
```bash
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache
```

#### 6. Load Nextflow Module

**On SAGA:**
```bash
module load Nextflow/24.04.2
```

**On TSD:**
```bash
module load Nextflow/24.04.2
```

#### 7. Run Pipeline

**Option A: Using test data graph (recommended - reuse downloaded graph):**
```bash
# Set cache location
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache

# Create production samples.csv with your data
cat > samples_production.csv << EOF
sample,read1,read2
patient001,/full/path/to/patient001_R1.fastq.gz,/full/path/to/patient001_R2.fastq.gz
patient002,/full/path/to/patient002_R1.fastq.gz,/full/path/to/patient002_R2.fastq.gz
EOF

# Run with downloaded test data graph
nextflow run main.nf \
    --input_csv samples_production.csv \
    --reference test_data/GRCh38.fa \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --outdir results \
    -profile saga \
    -resume
```

**Option B: With separately downloaded graph:**
```bash
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache

nextflow run main.nf \
    --input_csv samples.csv \
    --reference reference/GRCh38.fa \
    --hprc_graph graphs/hprc-v1.1-mc-grch38.gbz \
    --outdir results \
    -profile saga \
    -resume
```

**Option C: Auto-download graph (slower first run):**
```bash
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache

nextflow run main.nf \
    --input_csv samples.csv \
    --reference reference/GRCh38.fa \
    --outdir results \
    -profile saga \
    -resume
```

**Option D: With graph augmentation (add 1KGP SVs):**
```bash
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache

nextflow run main.nf \
    --input_csv samples.csv \
    --reference test_data/GRCh38.fa \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --kgp_sv_vcf 1KGP_SVs.vcf.gz \
    --augment_graph true \
    --outdir results \
    -profile saga \
    -resume
```

#### 8. For TSD: Create Production Sbatch Script

If running on TSD, modify the test script for production:

```bash
# Copy the test script
cp run_test_saga.sh run_production_tsd.sh

# Edit run_production_tsd.sh to change:
# Line 3: #SBATCH --account=p33_norment
# Line 56: -profile tsd
# Update --input_csv to your production samples.csv

# Then submit
sbatch run_production_tsd.sh
```

**Or run interactively:**
```bash
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache
module load Nextflow/24.04.2

nextflow run main.nf \
    --input_csv samples.csv \
    --reference test_data/GRCh38.fa \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --outdir results \
    -profile tsd \
    -resume
```

## Input Requirements

### Required

1. **Samples CSV** - Tab or comma-separated file with columns:
   - `sample`: Sample identifier (unique)
   - `read1`: Full path to R1 FASTQ file
   - `read2`: Full path to R2 FASTQ file

2. **Reference Genome** - GRCh38/hg38 FASTA file

### Optional

3. **HPRC Graph** - Will auto-download HPRC v1.1 if not provided
4. **1KGP SV VCF** - For graph augmentation (optional)

## Pipeline Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input_csv` | Path to samples CSV file | `samples.csv` |
| `--reference` | Path to GRCh38 FASTA file | `GRCh38.fa` |

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--hprc_graph` | Path to HPRC graph (GBZ) | Auto-download v1.1 |
| `--kgp_sv_vcf` | 1KGP SV VCF for augmentation | `null` |
| `--augment_graph` | Augment graph with 1KGP SVs | `false` |
| `--read_type` | Read type: `short` or `long` | `short` |
| `--output_bam` | Save BAM alignments | `false` |
| `--outdir` | Output directory | `./results` |
| `--cache_dir` | Singularity cache directory | `${launchDir}/singularity_cache` |

## Example Usage

### Basic Multi-Sample Run

```bash
nextflow run main.nf \
    --input_csv samples.csv \
    --reference GRCh38.fa \
    -profile saga
```

### With Local HPRC Graph (Recommended)

```bash
nextflow run main.nf \
    --input_csv samples.csv \
    --reference GRCh38.fa \
    --hprc_graph hprc-v1.1-mc-grch38.gbz \
    -profile saga
```

### With Graph Augmentation

```bash
nextflow run main.nf \
    --input_csv samples.csv \
    --reference GRCh38.fa \
    --hprc_graph hprc-v1.1-mc-grch38.gbz \
    --kgp_sv_vcf 1KGP_SVs.vcf.gz \
    --augment_graph true \
    -profile saga
```

### Long-Read Data (PacBio/ONT)

```bash
nextflow run main.nf \
    --input_csv samples_longread.csv \
    --reference GRCh38.fa \
    --read_type long \
    -profile tsd
```

### Resume Failed Run

```bash
nextflow run main.nf \
    --input_csv samples.csv \
    --reference GRCh38.fa \
    -profile saga \
    -resume
```

## Output Structure

```
results/
├── final/
│   ├── sample1.merged.vcf.gz              # All variants merged
│   ├── sample1.merged.vcf.gz.tbi          # Index
│   ├── sample1.snps_indels.vcf.gz         # SNPs/indels only
│   ├── sample1.snps_indels.vcf.gz.tbi     # Index
│   ├── sample1.svs.vcf.gz                 # SVs only
│   ├── sample1.svs.vcf.gz.tbi             # Index
│   ├── sample2.merged.vcf.gz              # Sample 2...
│   └── ...
│
├── deepvariant/
│   ├── sample1.deepvariant.vcf.gz         # SNPs/indels
│   ├── sample1.deepvariant.vcf.gz.tbi
│   ├── sample1.deepvariant.g.vcf.gz       # gVCF
│   └── sample1.deepvariant.g.vcf.gz.tbi
│
├── pangenie/
│   ├── sample1.pangenie.vcf.gz            # SV calls
│   └── sample1.pangenie.vcf.gz.tbi
│
├── alignments/                             # If --output_bam
│   ├── sample1.bam
│   ├── sample1.bam.bai
│   └── ...
│
├── graph/
│   └── hprc-v1.1-mc-grch38.gbz            # Downloaded graph
│
└── indexes/
    ├── graph.gbz                           # Indexed graph
    ├── graph.dist                          # Distance index
    └── graph.min                           # Minimizer index
```

### Output Files Per Sample

| File | Description |
|------|-------------|
| `sample.merged.vcf.gz` | All variants (SNPs, indels, SVs) |
| `sample.snps_indels.vcf.gz` | SNPs and indels only |
| `sample.svs.vcf.gz` | Structural variants only |
| `sample.bam` | Linear hg38 alignment (optional) |

## Software & Containers

All software runs in Singularity containers (automatically downloaded):

| Tool | Version | Container | Purpose |
|------|---------|-----------|---------|
| **vg toolkit** | v1.51.0 | `quay.io/vgteam/vg:v1.51.0` | Graph alignment, indexing |
| **DeepVariant** | pangenome-aware | `gcr.io/deepvariant-docker/deepvariant:pangenome_aware_deepvariant-head737001992` | SNP/indel calling |
| **PanGenie** | v3.0.0 | `quay.io/biocontainers/pangenie:3.0.0--h7d875b9_0` | SV genotyping |
| **samtools** | v1.18 | `quay.io/biocontainers/samtools:1.18--h50ea8bc_1` | BAM processing |
| **bcftools** | v1.18 | `quay.io/biocontainers/bcftools:1.18--h8b25389_0` | VCF manipulation |

Containers are cached in `${launchDir}/singularity_cache` to avoid filling `$HOME`.

## Resource Requirements

### Default Allocations

| Process | CPUs | Memory | Time | Notes |
|---------|------|--------|------|-------|
| Graph indexing | 16 | 64 GB | 12h | One-time per graph |
| Giraffe align | 16 | 64 GB | 12h | Per sample |
| DeepVariant | 16 | 64 GB | 24h | Per sample |
| PanGenie | 8 | 32 GB | 8h | Per sample |
| Postprocess | 4 | 16 GB | 4h | Per sample |

### Typical Runtime

For 30x WGS data per sample:
- Graph download/indexing: ~2-4 hours (one-time)
- Alignment: ~6-8 hours
- Variant calling: ~12-16 hours
- **Total: ~20-24 hours per sample**

Multiple samples run in parallel (resource dependent).

## HPC Configuration

### SAGA Profile

```groovy
executor: SLURM
account: nn9114k
queue: normal
memory: 32-64 GB per job
cpus: 8-16 per job
```

### TSD Profile

```groovy
executor: SLURM
account: p33_norment
perCpuMemAllocation: true
queue: normal
memory: 32-64 GB per job
cpus: 8-16 per job
```

## Troubleshooting

### Problem: Out of Memory

**Solution:** Edit `nextflow.config` and reduce memory:
```groovy
process {
    withLabel: 'vg' {
        memory = '32 GB'
    }
    withLabel: 'deepvariant' {
        memory = '48 GB'
    }
}
```

### Problem: Job Timeout

**Solution:** Increase time limits in `nextflow.config`:
```groovy
process {
    withLabel: 'deepvariant' {
        time = '48h'
    }
}
```

### Problem: Graph Download Failed

**Solution:** Download manually and provide path:
```bash
cd graphs
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz

# Then use:
nextflow run main.nf --hprc_graph graphs/hprc-v1.1-mc-grch38.gbz ...
```

### Problem: Singularity Cache Filling $HOME

**Solution:** Always set cache directory:
```bash
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache
mkdir -p singularity_cache
```

### Problem: Pipeline Crashes with "No such file"

**Solution:** Use absolute paths in `samples.csv`:
```bash
# Convert relative to absolute paths
realpath data/sample1_R1.fastq.gz
```

### Problem: Resume Not Working

**Solution:** Check work directory exists and use same parameters:
```bash
# Resume with exact same command
nextflow run main.nf --input_csv samples.csv ... -resume
```

## Important Notes

### Sample Names
- Must be unique across all samples
- Used as prefix for all output files
- Avoid special characters (use alphanumeric + underscore)

### File Paths
- Always use **absolute paths** in `samples.csv`
- Relative paths may cause issues with Nextflow work directories

### Singularity Cache
- Default location: `${launchDir}/singularity_cache`
- First run downloads ~5-10 GB of container images
- Subsequent runs reuse cached containers

### Graph Download
- HPRC v1.1 graph is ~50 GB
- First run downloads automatically (can take 1-2 hours)
- Subsequent runs reuse cached graph
- Recommend manual download for production use

### Disk Space
- Graph + indexes: ~60 GB
- Per sample working files: ~50-100 GB (temporary)
- Per sample output: ~2-5 GB (compressed VCFs)
- Plan accordingly for multiple samples!

## SLURM Batch Script Details

The repository includes `run_test_saga.sh` - a ready-to-use SLURM sbatch script for SAGA:

**Script configuration:**
- Job name: `nf_pangenome_test`
- Account: `nn9114k`
- Time limit: 24 hours
- Resources: 2 CPUs, 8 GB RAM (for Nextflow orchestration only)
- Singularity cache: Uses `singularity_cache/` in working directory

**What it does:**
1. Loads Nextflow/24.04.2 module
2. Sets Singularity cache directory
3. Verifies test data exists
4. Runs Nextflow pipeline with `-resume` support
5. Outputs log to `nf_pangenome_test_<jobid>.log`

**To customize for production:**
```bash
# Copy the script
cp run_test_saga.sh run_production.sh

# Edit these variables:
# - Change samples.csv to your data
# - Change reference to full GRCh38 (if not using chr22)
# - Change --outdir as needed
# - Adjust time/memory if needed

sbatch run_production.sh
```

**For TSD users:**
Create a similar script but change:
- Account: `p33_norment`
- Profile: `-profile tsd`

## FAQ

**Q: Can I use GRCh37/hg19 reference?**  
A: No, HPRC graphs are built on GRCh38. You must use GRCh38/hg38.

**Q: What coverage do I need?**  
A: Minimum 20x, recommended 30x for accurate variant calling.

**Q: Can I run single-end reads?**  
A: No, pipeline requires paired-end reads (R1 + R2).

**Q: How long does it take?**  
A: ~20-24 hours per sample (30x WGS) on SAGA/TSD.

**Q: Can I call variants in specific regions only?**  
A: Not currently supported. Pipeline processes whole genome.

**Q: What's the difference from linear reference pipeline?**  
A: Pangenome approach improves accuracy in complex regions, better SV calling, reduces reference bias.

**Q: Can I reuse test data for production?**  
A: Yes! The HPRC graph in `test_data/` can be used for production runs. Just point to it with `--hprc_graph test_data/hprc-v1.1-mc-grch38.gbz`

**Q: Do I need to download data every time?**  
A: No. Once downloaded, the pipeline reuses cached data. The graph, containers, and reference are all cached.

## Citations

If you use this pipeline, please cite the following papers:

### Core Tools

**vg toolkit and vg Giraffe:**
- Garrison, E., Sirén, J., Novak, A.M. et al. Variation graph toolkit improves read mapping by representing genetic variation in the reference. *Nature Biotechnology* **36**, 875–879 (2018). https://doi.org/10.1038/nbt.4227
- Sirén, J., Monlong, J., Chang, X. et al. Pangenomics enables genotyping of known structural variants in 5202 diverse genomes. *Science* **374**(6574), eabg8871 (2021). https://doi.org/10.1126/science.abg8871

**DeepVariant (original):**
- Poplin, R., Chang, P.C., Alexander, D. et al. A universal SNP and small-indel variant caller using deep neural networks. *Nature Biotechnology* **36**, 983–987 (2018). https://doi.org/10.1038/nbt.4235

**Pangenome-Aware DeepVariant:**
- See the latest documentation at: https://github.com/google/deepvariant/blob/r1.9/docs/deepvariant-vg-case-study.md
- Docker container: `gcr.io/deepvariant-docker/deepvariant:pangenome_aware_deepvariant`

**PanGenie:**
- Ebler, J., Ebert, P., Clarke, W.E. et al. Pangenome-based genome inference allows efficient and accurate genotyping across a wide spectrum of variant classes. *Nature Genetics* **54**, 518–525 (2022). https://doi.org/10.1038/s41588-022-01043-w

**HPRC v1.1 Pangenome:**
- Liao, W.W., Asri, M., Ebler, J. et al. A draft human pangenome reference. *Nature* **617**, 312–324 (2023). https://doi.org/10.1038/s41586-023-05896-x

### Data Resources

**HPRC Pangenome Resources:**
- GitHub: https://github.com/human-pangenomics/hpp_pangenome_resources
- Graph files: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html

**Tool Repositories:**
- vg toolkit: https://github.com/vgteam/vg
- DeepVariant: https://github.com/google/deepvariant
- PanGenie: https://github.com/eblerjana/pangenie
- Nextflow: https://www.nextflow.io/

## Helper Scripts in bin/ (Optional)

The `bin/` directory contains standalone helper scripts if you prefer script-based setup:

```bash
# Alternative to Step 2-3: Run helper scripts
bash bin/download_test_data.sh   # Downloads all test data
bash bin/prepare_test.sh          # Creates samples.csv and cache dir
bash bin/run_test_saga.sh         # Runs pipeline interactively
bash bin/run_test_tsd.sh          # Runs pipeline on TSD interactively
```

**Difference between `run_test_saga.sh` (top-level) and `bin/run_test_saga.sh`:**
- **Top-level `run_test_saga.sh`**: SLURM sbatch script for job submission
- **`bin/run_test_saga.sh`**: Interactive shell script (runs in current session)

Use the top-level sbatch script for actual HPC runs!

## Support

For issues and questions:
- Pipeline issues: Open GitHub issue at https://github.com/datngu/nf_pangenome
- SAGA support: support@metacenter.no
- TSD support: tsd-drift@usit.uio.no

## License

MIT License - See LICENSE file for details.

---

## Quick Reference Card

**Test run on SAGA:**
```bash
git clone https://github.com/datngu/nf_pangenome.git && cd nf_pangenome
# Download test data (Step 2), prepare CSV (Step 3)
sbatch run_test_saga.sh
```

**Production run on SAGA:**
```bash
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache
module load Nextflow/24.04.2
nextflow run main.nf --input_csv samples.csv --reference test_data/GRCh38.fa --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz -profile saga -resume
```

**Check results:**
```bash
ls -lh results_test/final/
zcat results_test/final/NA12878_200M.merged.vcf.gz | grep -v "^#" | wc -l
```
