# nf-pangenome

Nextflow pipelines for pangenome-based variant calling using the HPRC pangenome graph.

## Overview

This repository contains **two Nextflow pipelines** for state-of-the-art pangenome-based variant calling:

1. **`main_short_var.nf`** - Short variant calling (SNPs/indels) with vg Giraffe + Pangenome-Aware DeepVariant
2. **`main_sv.nf`** - Structural variant genotyping with PanGenie

Both pipelines use the [HPRC v1.1 pangenome graph](https://github.com/human-pangenomics/hpp_pangenome_resources) (GRCh38-based) for improved variant detection in diverse populations.

## Quick Start

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

# 1. Download pangenie_hprc.vcf (HPRC variants for PanGenie)
echo "Downloading pangenie_hprc.vcf..."
wget --progress=bar:force \
    -O pangenie_hprc.vcf.gz \
    "https://zenodo.org/record/6797328/files/cactus_filtered_ids.vcf.gz?download=1"

zcat pangenie_hprc.vcf.gz > pangenie_hprc.vcf
echo "✓ pangenie_hprc.vcf downloaded"

wget --progress=bar:force \
    -O pangenie_hprc_callset.vcf.gz \
    "https://zenodo.org/record/6797328/files/cactus_filtered_ids_biallelic.vcf.gz?download=1"

zcat pangenie_hprc_callset.vcf.gz > pangenie_hprc_callset.vcf

echo "✓ pangenie_hprc_callset.vcf downloaded"

# 2. Download HPRC v1.1 pangenome graph
echo "Downloading HPRC v1.1 graph..."
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

### Step 4: Configure Nextflow for Your HPC

Edit `nextflow.config` to add your cluster configuration. The file includes templates for SAGA, Orion, and a generic HPC setup.

**Example for adding your HPC:**

```groovy
profiles {
  saga {
    executor.name          = 'slurm'
    executor.account       = 'nn9114k'
    process.clusterOptions = '--time=72:00:00'  
  }

  orion {
    executor.name          = 'slurm'
    process.clusterOptions = '--time=72:00:00'  
  }
  
  // ADD your hpc setup below
  my_cluster {
    executor.name          = 'slurm'
    executor.account       = 'my_account_name'
    process.clusterOptions = '--time=72:00:00 --partition=normal'
    // executor.perCpuMemAllocation = true   
    // Uncomment if your Slurm requires --mem-per-cpu instead of --mem
  }

  singularity {
    singularity.enabled     = true
    singularity.autoMounts  = true
    singularity.pullTimeout = '5h'

    conda.enabled           = false
    docker.enabled          = false
  }
}
```

**Key configuration options:**
- `executor.name`: Job scheduler type (`slurm`, `pbs`, `sge`, etc.)
- `executor.account`: Your HPC project/account name
- `process.clusterOptions`: Additional scheduler flags (time limits, partitions, QoS, etc.)
- `executor.perCpuMemAllocation`: Uncomment if your cluster uses `--mem-per-cpu`

---

## Pipeline 1: Short Variant Calling (`main_short_var.nf`)

Call SNPs and indels using vg Giraffe alignment and Pangenome-Aware DeepVariant.

### Workflow

```
HPRC Graph (GBZ) → Extract Reference → Index Graph
                                            ↓
Sample Reads → Giraffe Align → Project BAM → DeepVariant → VCF (SNPs/indels)
```

### Run the Pipeline

```bash
# Using SAGA profile with Singularity
nextflow run main_short_var.nf \
    -profile saga,singularity \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --meta_csv samples.csv \
    --outdir results_short_variants \
    --output_bam true

# Using your custom profile
nextflow run main_short_var.nf \
    -profile my_cluster,singularity \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --meta_csv samples.csv \
    --outdir results_short_variants
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--hprc_graph` | `test_data/hprc-v1.1-mc-grch38.gbz` | HPRC pangenome graph (GBZ format) |
| `--meta_csv` | `samples.csv` | Sample sheet (sample,read1,read2) |
| `--outdir` | `test_results` | Output directory |
| `--read_type` | `short` | Read type: `short` or `long` |
| `--output_bam` | `false` | Keep BAM files (set `true` to retain) |
| `--ref_prefix` | `GRCh38` | HPRC reference path prefix (without #0#) |

### Output Structure

```
results_short_variants/
├── reference/              # Extracted reference from graph
│   ├── genome_ref.fa       # Reference with HPRC notation (GRCh38#0#chr1)
│   ├── genome_ref.fa.fai
│   └── chrom_rename.txt    # Map for converting HPRC → standard notation
├── indexes/                # Graph indexes (reusable)
│   ├── index.gbz
│   ├── index.dist          # Distance index
│   ├── index.min           # Minimizer index
│   └── index.zipcodes      # Zipcode index for Giraffe
├── alignments/             # GAM files and optional BAMs
│   └── {sample}/
│       ├── {sample}.gam
│       ├── {sample}.bam         # Only if --output_bam true
│       └── {sample}.bam.bai     # Only if --output_bam true
└── variants/               # VCF files (HPRC notation)
    └── {sample}/
        ├── {sample}.hprc.vcf.gz        # VCF with HPRC chromosome names
        ├── {sample}.hprc.vcf.gz.tbi
        ├── {sample}.hprc.g.vcf.gz      # gVCF with HPRC chromosome names
        └── {sample}.hprc.g.vcf.gz.tbi
```

**Note:** The FIX_VCF_CHROMS process is currently commented out. VCF files retain HPRC notation (GRCh38#0#chr1).

---

## Pipeline 2: SV Genotyping (`main_sv.nf`)

Genotype structural variants using PanGenie with HPRC SV catalog.

### Workflow

```
HPRC Graph (GBZ) → Extract Reference
                        ↓
HPRC VCF → PanGenie Index (once)
                ↓
Sample Reads → PanGenie Genotype → VCF (SVs)
```

### Run the Pipeline

```bash
# Using SAGA profile with Singularity
nextflow run main_sv.nf \
    -profile saga,singularity \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --hprc_pangenie_vcf test_data/pangenie_hprc.vcf \
    --meta_csv samples.csv \
    --outdir results_svs

# Using your custom profile
nextflow run main_sv.nf \
    -profile my_cluster,singularity \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --hprc_pangenie_vcf test_data/pangenie_hprc.vcf \
    --meta_csv samples.csv \
    --outdir results_svs
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--hprc_graph` | `test_data/hprc-v1.1-mc-grch38.gbz` | HPRC pangenome graph |
| `--hprc_pangenie_vcf` | `test_data/pangenie_hprc.vcf` | HPRC SV catalog for PanGenie |
| `--hprc_pangenie_callset_vcf` | `test_data/pangenie_hprc_callset.vcf` | Biallelic HPRC SV catalog |
| `--meta_csv` | `samples.csv` | Sample sheet (sample,read1,read2) |
| `--outdir` | `test_results` | Output directory |
| `--ref_prefix` | `GRCh38` | Reference prefix for extraction |

**Alternative SV Catalog:** The `data/` directory contains the 1000 Genomes Project ONT phased SV callset (`1000GP_ONT_shapeit5-phased-callset_final-vcf.phased.vcf.gz`), which can be used instead of the HPRC catalog. Integration is planned (see TODO list).

### Output Structure

```
results_svs/
├── reference/              # Extracted reference (standard chr notation)
│   ├── genome_ref.fa       # Reference without HPRC notation (chr1, chr2...)
│   ├── genome_ref.fa.fai
│   └── chrom_rename.txt
├── pangenie_index/         # PanGenie index (reusable)
│   └── pangenie_index*     # Multiple index files
└── genotypes/              # SV genotype VCFs
    ├── {sample}_genotyping.vcf.gz
    └── {sample}_genotyping.vcf.gz.tbi
```

---

## Profile Usage

Nextflow profiles control where and how your pipeline runs. You can combine multiple profiles:

### Common Profile Combinations

```bash
# SAGA HPC with Singularity containers
nextflow run main_short_var.nf -profile saga,singularity ...

# Your custom HPC with Singularity
nextflow run main_sv.nf -profile my_cluster,singularity ...

# Orion HPC with Singularity
nextflow run main_short_var.nf -profile orion,singularity ...

```

### Profile Types

1. **Executor profiles** (`saga`, `orion`, `my_cluster`): Define HPC-specific settings
2. **Container profiles** (`singularity`): Enable container technology

**Always use both** an executor profile and the `singularity` profile for HPC runs.

---

## Sample CSV Format

Both pipelines require a CSV file with three columns:

```csv
sample,read1,read2
sample1,/path/to/sample1_R1.fq.gz,/path/to/sample1_R2.fq.gz
sample2,/path/to/sample2_R1.fq.gz,/path/to/sample2_R2.fq.gz
```

**Requirements:**
- Use **absolute paths** for read files
- Gzipped FASTQ files (`.fq.gz` or `.fastq.gz`)
- Paired-end reads required

---

## Tools and Containers

All tools are containerized for reproducibility and run via Singularity on HPC systems.

### Short Variant Pipeline (`main_short_var.nf`)

| Process | Container | Tool | Version |
|---------|-----------|------|---------|
| EXTRACT_REFERENCE | `quay.io/vgteam/vg:v1.65.0` | vg paths + samtools | v1.65.0 |
| EXTRACT_REFERENCE | `quay.io/vgteam/vg:v1.65.0` | vg paths + samtools | v1.65.0 |
| INDEX_GRAPH | `quay.io/vgteam/vg:v1.65.0` | vg index + vg minimizer | v1.65.0 |
| GIRAFFE_ALIGN | `quay.io/vgteam/vg:v1.65.0` | vg giraffe | v1.65.0 |
| PROJECT_BAM | `quay.io/vgteam/vg:v1.65.0` | vg surject + samtools | v1.65.0 |
| CALL_SNPS_INDELS | `google/deepvariant:pangenome_aware_deepvariant-1.8.0` | Pangenome-Aware DeepVariant | 1.8.0 |
| FIX_VCF_CHROMS* | `biocontainers/bcftools:1.20--h8b25389_0` | bcftools annotate | 1.20 |

*Currently commented out in workflow
### SV Pipeline (`main_sv.nf`)

| Process | Container | Tool | Version |
|---------|-----------|------|---------|
| EXTRACT_REFERENCE | `quay.io/vgteam/vg:v1.65.0` | vg paths + samtools | v1.65.0 |
| PANGENIE_INDEX | `mgibio/pangenie:v4.2.1-bookworm` | PanGenie-index | v4.2.1 |
| PANGENIE_GENOTYPE | `mgibio/pangenie:v4.2.1-bookworm` | PanGenie | v4.2.1 |

---

## Resource Requirements

### Short Variant Pipeline (`main_short_var.nf`)

| Process | Memory | CPUs | Time (est.) |
|---------|--------|------|-------------|
| EXTRACT_REFERENCE | 32 GB | 8 | ~30 min |
| INDEX_GRAPH | 128 GB | 16 | ~2-4 hours |
| GIRAFFE_ALIGN | 64 GB | 16 | ~4-8 hours per sample (30x), ~6-12 hours (50x) |
| PROJECT_BAM | 128 GB | 16 | ~2-4 hours per sample (30x), ~3-6 hours (50x) |
| CALL_SNPS_INDELS | 64 GB | 16 | ~8-12 hours per sample |
| FIX_VCF_CHROMS | 16 GB | 4 | ~30 min per sample |
| CALL_SNPS_INDELS | 64 GB | 16 | ~8-12 hours per sample (30x), ~16-20 hours (50x) |
| FIX_VCF_CHROMS | 16 GB | 4 | ~30 min per sample (currently disabled) |

| Process | Memory | CPUs | Time (est.) |
|---------|--------|------|-------------|
| EXTRACT_REFERENCE | 32 GB | 8 | ~30 min |
| PANGENIE_INDEX | 128 GB | 16 | ~2-4 hours |
| PANGENIE_GENOTYPE | 64 GB | 16 | ~4-8 hours per sample |

**Note:** Times are estimates for 30x WGS data. Actual times vary based on coverage and cluster load.
**Note:** Times are estimates for 30x WGS data. For 50x WGS, expect ~1.5-2x longer runtimes. Actual times vary based on coverage and cluster load.

---

## TODO List

### High Priority

#### 1. Post-Processing for Variant Calls
- [ ] **PanGenie Output Processing**
  - Add filtering step for low-confidence genotypes
  - Implement coverage-based quality filtering
  - Add process for merging multi-sample VCFs
  - Generate summary statistics (GT counts, quality metrics)
  
- [ ] **DeepVariant Output Processing**
  - Re-enable `FIX_VCF_CHROMS` process to convert HPRC notation to standard (chr1, chr2...)
  - Add BCFtools norm for left-alignment and splitting multiallelic sites
  - Implement VQSR or hard filtering based on QUAL/GQ thresholds
  - Add variant annotation (e.g., with VEP or SnpEff)

#### 2. Resource Optimization for 50x WGS
- [ ] **Test and Optimize Memory Allocation**
  - Test `main_short_var.nf` with 50x WGS samples
  - Adjust `PROJECT_BAM` memory (currently 128 GB, may need 192 GB for 50x)
  - Adjust `CALL_SNPS_INDELS` memory (currently 64 GB, may need 96-128 GB for 50x)
  - Test `PANGENIE_GENOTYPE` with 50x data (currently 64 GB)
  
- [ ] **CPU Optimization**
  - Profile CPU usage for each process
  - Optimize samtools threading in `PROJECT_BAM` (currently 3 threads)
  - Test different `--num_shards` values for DeepVariant (currently = cpus)

#### 3. T2T Pangenome Support
- [ ] **Add T2T-CHM13 Backbone Support**
  - Download T2T-based HPRC graph: `hprc-v1.1-mc-chm13.gbz`
  - Add `--ref_prefix` parameter validation (GRCh38 vs CHM13)
  - Update `EXTRACT_REFERENCE` to handle T2T path names
  - Test chromosome naming conventions for T2T
  - Document T2T-specific parameters in README

### Medium Priority

#### 4. Update to Newer HPRC Release
- [ ] **Short Variant Pipeline (Relatively Easy)**
  - Update default graph parameter to newer HPRC version when available
  - Test compatibility with HPRC v1.2+ graphs
  - Update documentation with new download URLs
  - Validate that Giraffe indexes work with newer graphs
  
- [ ] **SV Pipeline (More Complex)**
  - Check if newer HPRC VCF format requires parser updates
  - Test PanGenie compatibility with updated SV catalogs
  - Update `hprc_pangenie_vcf` parameter to newer version
  - Validate SV genotyping accuracy with updated catalogs
  - May require re-indexing with PanGenie-index

#### 5. Integrate 1000GP ONT SV Catalog
**HPRC SV Catalog (used by PanGenie):**
- Ebert P, Audano PA, Zhu Q, et al. Haplotype-resolved diverse human genomes and integrated analysis of structural variation. *Science*. 2021;372(6537):eabf7117. doi:[10.1126/science.abf7117](https://doi.org/10.1126/science.abf7117)

**1000 Genomes Project ONT SVs (available in `data/` directory):**
- Ebert P, Audano PA, Zhu Q, et al. Haplotype-resolved diverse human genomes and integrated analysis of structural variation. *Nature*. 2025;637:696-706. doi:[10.1038/s41586-025-09290-7](https://doi.org/10.1038/s41586-025-09290-7)
- [ ] **Add 1KGP ONT VCF Support**
  - Create parameter: `--sv_vcf` (default: 1KGP ONT phased VCF from `data/`)
  - Add option to choose between HPRC VCF and 1KGP ONT VCF
  - Implement VCF preprocessing:
    - Filter for high-quality SVs (PASS only)
    - Convert phased genotypes if needed
    - Ensure chromosome naming consistency (chr1 vs 1)
  
- [ ] **Update PanGenie Workflow**
  - Modify `PANGENIE_INDEX` to accept custom VCF input
  - Add validation step for VCF format compatibility
  - Test with 1KGP ONT SV dataset (file exists in `data/1000GP_ONT_shapeit5-phased-callset_final-vcf.phased.vcf.gz`)
  - Compare genotyping results between HPRC and 1KGP catalogs
  
- [ ] **Create Hybrid Workflow**
  - Option to merge HPRC and 1KGP SV catalogs
  - Deduplicate overlapping SVs between catalogs
  - Create unified VCF for PanGenie indexing
  - Document merged catalog usage and performance

### Low Priority

#### 6. Additional Enhancements
- [ ] Add multi-sample joint genotyping for DeepVariant (GLnexus)
- [ ] Implement per-chromosome parallelization for faster processing
- [ ] Add QC metrics dashboard (MultiQC integration)
- [ ] Support for long-read data (PacBio HiFi, ONT)
- [ ] Add benchmark datasets and validation pipeline
- [ ] Create Docker/Singularity recipe for custom container builds

---
## Troubleshooting

### Common Issues

**1. Singularity cache errors**
```bash
# Clear Singularity cache if containers fail to pull
rm -rf $HOME/.singularity/cache
```

**2. Out of memory errors**
- Increase memory in `nextflow.config` or process directives
- Check cluster memory limits

**3. Path not found errors**
- Ensure absolute paths in `samples.csv`
- Verify data files exist before running

**4. Profile not found**
```bash
# List available profiles
grep "profiles {" -A 30 nextflow.config

# Always combine executor + singularity profiles
nextflow run main_short_var.nf -profile saga,singularity ...
```

**5. Container permission errors**
- Ensure `singularity.autoMounts = true` in config
- Check file permissions on input data

---

## Citations

If you use these pipelines, please cite:

### Pangenome Graph and HPRC

**HPRC Pangenome:**
- Liao WW, Asri M, Ebler J, et al. A draft human pangenome reference. *Nature*. 2023;617(7960):312-324. doi:[10.1038/s41586-023-05896-x](https://doi.org/10.1038/s41586-023-05896-x)

### Short Variant Calling Tools

**vg Giraffe (Graph Alignment):**
- Monlong J, Ebler J, Eizenga JM, et al. Pangenome graph construction from genome sequences with structural variants. *bioRxiv*. 2025. doi:[10.1101/2025.09.29.678807](https://doi.org/10.1101/2025.09.29.678807)

**Pangenome-Aware DeepVariant:**
⚠️ **This is NOT standard DeepVariant** - it's specifically trained for pangenome graphs:
- Chang PC, Bazak BR, Nattestad M, Coddington R, Eizenga JM, Carroll A, Paten B. Pangenome-Aware DeepVariant: Accurate SNP and Indel Calling on Whole Genome Sequences Aligned to Pangenome Reference Graphs. *medRxiv*. 2025. URL: [pubmed.ncbi.nlm.nih.gov/40501862](https://pubmed.ncbi.nlm.nih.gov/40501862)

### SV Genotyping Tools

**PanGenie:**
- Ebler J, Ebert P, Clarke WE, et al. Pangenome-based genome inference allows efficient and accurate genotyping across a wide spectrum of variant classes. *Nature Genetics*. 2022;54(4):518-525. doi:[10.1038/s41588-022-01043-w](https://doi.org/10.1038/s41588-022-01043-w)


### Nextflow

- Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. *Nature Biotechnology*. 2017;35(4):316-319. doi:[10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

---

## License

This pipeline is released under the MIT License.

---

## Author

**Dat T. Nguyen**  
Contact: ndat@utexas.edu  
GitHub: [datngu](https://github.com/datngu)

---

## Acknowledgments

- [Human Pangenome Reference Consortium (HPRC)](https://humanpangenome.org/)
- [vg toolkit team](https://github.com/vgteam/vg)
- [PanGenie developers](https://github.com/eblerjana/pangenie)
- [Google DeepVariant team](https://github.com/google/deepvariant)

---

## Additional Resources

- **HPRC Resources**: https://github.com/human-pangenomics/hpp_pangenome_resources
- **vg Documentation**: https://github.com/vgteam/vg/wiki
- **PanGenie Documentation**: https://github.com/eblerjana/pangenie
- **Nextflow Documentation**: https://www.nextflow.io/docs/latest/