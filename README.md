# Pangenome Variant Calling Pipeline

A streamlined Nextflow pipeline for pangenome-based variant calling using HPRC v1.1 graphs, pangenome-aware DeepVariant, and PanGenie. Designed for SAGA and TSD HPC environments with Singularity containers.

## Overview

This pipeline performs comprehensive variant calling (SNPs, indels, and SVs) by aligning sequencing reads to the Human Pangenome Reference Consortium (HPRC) graph instead of a linear reference genome. This approach improves variant calling accuracy, especially in complex genomic regions and for structural variants.

**Key Features:**
- Multi-sample processing via CSV input
- Pangenome graph alignment with vg giraffe
- Pangenome-aware variant calling with DeepVariant
- SV genotyping with PanGenie
- Automatic HPRC v1.1 graph download
- Singularity containers for reproducibility
- Optimized for SAGA and TSD HPC systems

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

## Quick Start

### Prepare Input CSV

Create a `samples.csv` file with your samples:
```csv
sample,read1,read2
sample1,/path/to/sample1_R1.fq.gz,/path/to/sample1_R2.fq.gz
sample2,/path/to/sample2_R1.fq.gz,/path/to/sample2_R2.fq.gz
sample3,/path/to/sample3_R1.fq.gz,/path/to/sample3_R2.fq.gz
```

### Basic Run on SAGA

```bash
# Set cache directory to avoid filling $HOME
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache
mkdir -p singularity_cache

# Load Nextflow
module load Nextflow/24.04.2

# Run pipeline
nextflow run main.nf \
    --input_csv samples.csv \
    --reference GRCh38.fa \
    --outdir results \
    -profile saga
```

### Basic Run on TSD

```bash
# Set cache directory
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache
mkdir -p singularity_cache

# Load Nextflow
module load Nextflow/24.04.2

# Run pipeline
nextflow run main.nf \
    --input_csv samples.csv \
    --reference GRCh38.fa \
    --outdir results \
    -profile tsd
```

## Production Data Setup

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

# Set environment variable
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity_cache

# Save to file for later use
echo "export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache" > run_env.sh
echo "export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity_cache" >> run_env.sh
chmod +x run_env.sh
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

**Basic run with auto-downloaded graph:**
```bash
source run_env.sh

nextflow run main.nf \
    --input_csv samples.csv \
    --reference reference/GRCh38.fa \
    --outdir results \
    -profile saga \
    -resume
```

**With local graph (faster):**
```bash
source run_env.sh

nextflow run main.nf \
    --input_csv samples.csv \
    --reference reference/GRCh38.fa \
    --hprc_graph graphs/hprc-v1.1-mc-grch38.gbz \
    --outdir results \
    -profile saga \
    -resume
```

**With graph augmentation (1KGP SVs):**
```bash
source run_env.sh

nextflow run main.nf \
    --input_csv samples.csv \
    --reference reference/GRCh38.fa \
    --hprc_graph graphs/hprc-v1.1-mc-grch38.gbz \
    --kgp_sv_vcf 1KGP_SVs.vcf.gz \
    --augment_graph true \
    --outdir results \
    -profile saga \
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

## Citation

If you use this pipeline, please cite:

- **vg toolkit**: Garrison et al. (2018) Nature Biotechnology
- **DeepVariant**: Poplin et al. (2018) Nature Biotechnology  
- **PanGenie**: Ebler et al. (2022) Nature Genetics
- **HPRC**: Liao et al. (2023) Nature

## References & Links

- **HPRC Resources**: https://github.com/human-pangenomics/hpp_pangenome_resources
- **vg toolkit**: https://github.com/vgteam/vg
- **DeepVariant**: https://github.com/google/deepvariant
- **PanGenie**: https://github.com/eblerjana/pangenie
- **Nextflow**: https://www.nextflow.io/

## Support

For issues and questions:
- Pipeline issues: Open GitHub issue
- SAGA support: support@metacenter.no
- TSD support: tsd-drift@usit.uio.no

## License

MIT License - See LICENSE file for details.
