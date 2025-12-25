# nf-pangenome

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.0-23aa62.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Version 1.0** - HPRC v1.1 Pangenome Reference Support

Nextflow pipelines for pangenome-based variant calling using the HPRC v1.1 pangenome graph.

## Overview

This repository contains **two production-ready Nextflow pipelines** for state-of-the-art pangenome-based variant calling:

1. **`main_short_var.nf`** - Short variant calling (SNPs/indels) with vg Giraffe + Pangenome-Aware DeepVariant
2. **`main_sv.nf`** - Structural variant genotyping with PanGenie

Both pipelines use the [HPRC v1.1 pangenome graph](https://github.com/human-pangenomics/hpp_pangenome_resources) (GRCh38-based) for improved variant detection in diverse populations.

### Key Features

- ✅ **HPRC v1.1 Support**: Full integration with latest HPRC pangenome graphs
- ✅ **Containerized**: All tools packaged in Docker/Singularity containers
- ✅ **HPC-Ready**: Optimized for SLURM clusters with configurable profiles
- ✅ **Biallelic Conversion**: Automatic post-processing for PanGenie outputs
- ✅ **Scalable**: Parallel processing of multiple samples
- ✅ **Reproducible**: Version-controlled containers and workflows

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
Sample Reads → Giraffe Align (GAM) → Project BAM → DeepVariant → VCF (SNPs/indels)
```

**Processes:**
1. **EXTRACT_REFERENCE** - Extract GRCh38 reference from graph (keeps HPRC notation)
2. **INDEX_GRAPH** - Build Giraffe indexes (distance, minimizer, zipcodes)
3. **GIRAFFE_ALIGN** - Align reads to pangenome with vg Giraffe
4. **PROJECT_BAM** - Project GAM to BAM with vg surject
5. **CALL_SNPS_INDELS** - Call variants with Pangenome-Aware DeepVariant
6. **FIX_VCF_CHROMS** - Convert HPRC notation to standard chromosome names

### Run the Pipeline

```bash
# Using SAGA profile with Singularity
nextflow run main_short_var.nf \
    -profile saga,singularity \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --meta_csv samples.csv \
    --outdir results_short_variants \
    --output_bam true \
    --read_type short

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
└── variants/               # VCF files (standard chromosome notation)
    └── {sample}/
        ├── {sample}.hprc.vcf.gz              # Intermediate VCF with HPRC notation
        ├── {sample}.hprc.vcf.gz.tbi
        ├── {sample}.hprc.g.vcf.gz            # Intermediate gVCF with HPRC notation
        ├── {sample}.hprc.g.vcf.gz.tbi
        ├── {sample}.deepvariant.vcf.gz       # Final VCF with standard chr names
        ├── {sample}.deepvariant.vcf.gz.tbi
        ├── {sample}.deepvariant.g.vcf.gz     # Final gVCF with standard chr names
        └── {sample}.deepvariant.g.vcf.gz.tbi
```

**Note:** The FIX_VCF_CHROMS process converts HPRC notation (GRCh38#0#chr1) to standard notation (chr1).

---

## Pipeline 2: SV Genotyping (`main_sv.nf`)

Genotype structural variants using PanGenie with HPRC SV catalog.

### Workflow

```
HPRC Graph (GBZ) → Extract Reference (standard chr notation)
                        ↓
HPRC Multiallelic VCF → PanGenie Index (once)
                            ↓
Sample Reads → PanGenie Genotype → Convert to Biallelic → Final VCF (SVs)
```

**Processes:**
1. **EXTRACT_REFERENCE** - Extract GRCh38 reference (removes HPRC notation)
2. **PANGENIE_INDEX** - Build PanGenie index from multiallelic VCF
3. **PANGENIE_GENOTYPE** - Genotype SVs per sample
4. **CONVERT_TO_BIALLELIC** - Post-process to biallelic format using `bin/convert-to-biallelic.py`

### Run the Pipeline

```bash
# Using SAGA profile with Singularity
nextflow run main_sv.nf \
    -profile saga,singularity \
    --hprc_graph test_data/hprc-v1.1-mc-grch38.gbz \
    --hprc_pangenie_vcf test_data/pangenie_hprc.vcf \
    --hprc_pangenie_callset_vcf test_data/pangenie_hprc_callset.vcf \
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
├── genotypes/              # Raw multiallelic genotype VCFs
│   ├── {sample}_genotyping.vcf.gz
│   └── {sample}_genotyping.vcf.gz.tbi
└── genotypes_biallelic/    # Biallelic genotype VCFs (final output)
    ├── {sample}_genotyping_biallelic.vcf.gz
    └── {sample}_genotyping_biallelic.vcf.gz.tbi
```

**Note:** The biallelic conversion step uses `bin/convert-to-biallelic.py` which adaptively handles both compressed (.gz) and uncompressed VCF files.

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
| INDEX_GRAPH | `quay.io/vgteam/vg:v1.65.0` | vg index + vg minimizer | v1.65.0 |
| GIRAFFE_ALIGN | `quay.io/vgteam/vg:v1.65.0` | vg giraffe | v1.65.0 |
| PROJECT_BAM | `quay.io/vgteam/vg:v1.65.0` | vg surject + samtools | v1.65.0 |
| CALL_SNPS_INDELS | `google/deepvariant:pangenome_aware_deepvariant-1.8.0` | Pangenome-Aware DeepVariant | 1.8.0 |
| FIX_VCF_CHROMS | `biocontainers/bcftools:1.20--h8b25389_0` | bcftools annotate | 1.20 |

### SV Pipeline (`main_sv.nf`)

| Process | Container | Tool | Version |
|---------|-----------|------|---------|
| EXTRACT_REFERENCE | `quay.io/vgteam/vg:v1.65.0` | vg paths + samtools | v1.65.0 |
| PANGENIE_INDEX | `mgibio/pangenie:v4.2.1-bookworm` | PanGenie-index | v4.2.1 |
| PANGENIE_GENOTYPE | `mgibio/pangenie:v4.2.1-bookworm` | PanGenie | v4.2.1 |
| CONVERT_TO_BIALLELIC | `python:3.9-slim` | Custom Python script | v1.0 |

---

## Resource Requirements

### Short Variant Pipeline (`main_short_var.nf`)

| Process | Memory | CPUs | Time (est.) |
|---------|--------|------|-------------|
| EXTRACT_REFERENCE | 32 GB | 8 | ~30 min |
| INDEX_GRAPH | 128 GB | 16 | ~2-4 hours |
| GIRAFFE_ALIGN | 64 GB | 16 | ~4-8 hours per sample (30x), ~6-12 hours (50x) |
| PROJECT_BAM | 128 GB | 16 | ~2-4 hours per sample (30x), ~3-6 hours (50x) |
| CALL_SNPS_INDELS | 64 GB | 16 | ~8-12 hours per sample (30x), ~16-20 hours (50x) |
| FIX_VCF_CHROMS | 16 GB | 4 | ~30 min per sample |

### SV Pipeline (`main_sv.nf`)

| Process | Memory | CPUs | Time (est.) |
|---------|--------|------|-------------|
| EXTRACT_REFERENCE | 32 GB | 8 | ~30 min |
| PANGENIE_INDEX | 128 GB | 16 | ~2-4 hours |
| PANGENIE_GENOTYPE | 64 GB | 16 | ~4-8 hours per sample |
| CONVERT_TO_BIALLELIC | 16 GB | 4 | ~30 min per sample |

**Note:** Times are estimates for 30x WGS data. For 50x WGS, expect ~1.5-2x longer runtimes. Actual times vary based on coverage and cluster load.

---

## Release Notes

### Version 1.0 (December 2025)

**Initial Release - HPRC v1.1 Support**

#### Features
- ✅ Short variant calling pipeline with vg Giraffe and Pangenome-Aware DeepVariant
- ✅ SV genotyping pipeline with PanGenie
- ✅ Biallelic conversion post-processing for PanGenie outputs
- ✅ HPRC v1.1 pangenome graph support (GRCh38-based)
- ✅ Containerized workflows (Docker/Singularity)
- ✅ HPC cluster support with SLURM profiles
- ✅ Comprehensive documentation and test scripts

#### Known Limitations
- T2T-CHM13 pangenome support not yet implemented
- 1000GP ONT SV catalog integration pending

#### Test Data
- HPRC v1.1 MC GRCh38 graph: 50GB
- NA12878 test reads: 200M and 75M coverage datasets
- HPRC SV catalogs for PanGenie

---

## TODO List

### High Priority - v1.1 Roadmap

#### 1. Post-Processing Enhancements
- [x] **PanGenie Biallelic Conversion** ✅ Implemented in v1.0
  - Automatic conversion of multiallelic to biallelic format
  - Adaptive VCF reader (compressed/uncompressed)
  
- [ ] **DeepVariant VCF Post-Processing**
  - Add BCFtools norm for left-alignment and splitting multiallelic sites
  - Implement VQSR or hard filtering based on QUAL/GQ thresholds
  - Add variant annotation (e.g., with VEP or SnpEff)

- [ ] **Quality Control Metrics**
  - Add filtering step for low-confidence genotypes in PanGenie
  - Coverage-based quality filtering
  - Generate summary statistics (GT counts, quality metrics)
  - Multi-sample VCF merging support

#### 2. Resource Optimization for High-Coverage Data
- [ ] **50x WGS Testing and Optimization**
  - Benchmark `main_short_var.nf` with 50x WGS samples
  - Adjust `PROJECT_BAM` memory (may need 192 GB for 50x)
  - Adjust `CALL_SNPS_INDELS` memory (may need 96-128 GB for 50x)
  - Test `PANGENIE_GENOTYPE` with 50x data
  
- [ ] **CPU and Threading Optimization**
  - Profile CPU usage for each process
  - Optimize samtools threading in `PROJECT_BAM`
  - Test different `--num_shards` values for DeepVariant

#### 3. T2T Pangenome Support
- [ ] **T2T-CHM13 Backbone Implementation**
  - Download and test T2T-based HPRC graph: `hprc-v1.1-mc-chm13.gbz`
  - Add `--ref_prefix` parameter validation (GRCh38 vs CHM13)
  - Update `EXTRACT_REFERENCE` to handle T2T path names
  - Test chromosome naming conventions for T2T
  - Document T2T-specific parameters and usage

### Medium Priority - v1.2+ Roadmap

#### 4. HPRC Version Updates
- [ ] **Short Variant Pipeline** (Relatively Straightforward)
  - Update to HPRC v1.2+ graphs when available
  - Test compatibility with newer graph formats
  - Update documentation with new download URLs
  - Validate Giraffe index compatibility
  
- [ ] **SV Pipeline** (More Complex)
  - Test PanGenie with updated HPRC SV catalogs
  - Update `hprc_pangenie_vcf` to newer versions
  - Validate SV genotyping accuracy
  - May require PanGenie re-indexing

#### 5. 1000 Genomes Project ONT SV Integration
**Note:** Most complex task - requires significant workflow changes

- [ ] **1KGP ONT VCF Support**
  - Add `--sv_vcf` parameter (default: 1KGP ONT phased VCF from `data/`)
  - Option to choose between HPRC VCF and 1KGP ONT VCF
  - Implement VCF preprocessing:
    - Filter for high-quality SVs (PASS only)
    - Convert phased genotypes if needed
    - Ensure chromosome naming consistency
  
- [ ] **PanGenie Workflow Updates**
  - Modify `PANGENIE_INDEX` to accept custom VCF input
  - Add VCF format validation step
  - Test with 1KGP ONT dataset (`data/1000GP_ONT_shapeit5-phased-callset_final-vcf.phased.vcf.gz`)
  - Compare results: HPRC vs 1KGP catalogs
  
- [ ] **Hybrid Catalog Support**
  - Merge HPRC and 1KGP SV catalogs
  - Deduplicate overlapping SVs
  - Create unified VCF for PanGenie
  - Document performance and accuracy

### Low Priority - Future Enhancements

#### 6. Advanced Features
- [ ] Multi-sample joint genotyping for DeepVariant (GLnexus integration)
- [ ] Per-chromosome parallelization for faster processing
- [ ] MultiQC integration for QC dashboard
- [ ] Native long-read support (PacBio HiFi, ONT)
- [ ] Benchmark datasets and validation pipeline
- [ ] Custom container builds (Docker/Singularity recipes)
- [ ] Support for additional pangenome backbones (e.g., custom populations)

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

- *Nature* doi:[10.1038/s41586-023-05896-x](https://doi.org/10.1038/s41586-023-05896-x)

### Short Variant Calling Tools

**vg Giraffe (Graph Alignment):**

- *bioRxiv* doi:[10.1101/2025.09.29.678807](https://doi.org/10.1101/2025.09.29.678807)

**Pangenome-Aware DeepVariant:**
⚠️ **This is NOT standard DeepVariant** - it's specifically trained for pangenome graphs:
- *bioRxiv* doi:[10.1101/2025.06.05.657102](https://doi.org/10.1101/2025.06.05.657102)

### SV Genotyping Tools

**PanGenie:**
- *Nature Genetics*. doi:[10.1038/s41588-022-01043-w](https://doi.org/10.1038/s41588-022-01043-w)


**1000 Genomes Project ONT SVs (available in `data/` directory):**
- *Nature* doi:[10.1038/s41586-025-09290-7](https://doi.org/10.1038/s41586-025-09290-7)

### Nextflow

- *Nature Biotechnology* doi:[10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### Development Guidelines
- Follow existing code style and conventions
- Add tests for new features
- Update documentation accordingly
- Use semantic versioning for releases

---

## Support and Contact

**Issues:** Please report bugs and feature requests via [GitHub Issues](https://github.com/datngu/nf_pangenome/issues)

**Discussions:** For questions and discussions, use [GitHub Discussions](https://github.com/datngu/nf_pangenome/discussions)

---

## License

This pipeline is released under the MIT License. See [LICENSE](LICENSE) file for details.

---

## Author

**Dat T. Nguyen**  
Contact: ndat@utexas.edu  
GitHub: [datngu](https://github.com/datngu)  

---

## Acknowledgments

- [Human Pangenome Reference Consortium (HPRC)](https://humanpangenome.org/) for the pangenome graphs and resources
- [vg toolkit team](https://github.com/vgteam/vg) for graph genome tools
- [PanGenie developers](https://github.com/eblerjana/pangenie) for SV genotyping software
- [Google DeepVariant team](https://github.com/google/deepvariant) for variant calling tools
- SAGA HPC cluster (University of Oslo) and TSD secure computing platform
- Nextflow community for workflow management framework

---

## Additional Resources

- **HPRC Resources**: https://github.com/human-pangenomics/hpp_pangenome_resources
- **vg Documentation**: https://github.com/vgteam/vg/wiki
- **PanGenie Documentation**: https://github.com/eblerjana/pangenie
- **Nextflow Documentation**: https://www.nextflow.io/docs/latest/
- **Pangenome-Aware DeepVariant**: https://github.com/google/deepvariant/blob/r1.8/docs/pangenome-aware-variant-calling.md

---

**Citation:** If you use this pipeline in your research, please cite the relevant tools and the HPRC pangenome as described in the Citations section above.
