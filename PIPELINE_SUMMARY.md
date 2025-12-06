# Pangenome Variant Calling Pipeline - Technical Summary

## Pipeline Version
- **Version**: 1.0
- **Last Updated**: December 6, 2025
- **Nextflow DSL**: 2
- **Minimum Nextflow**: 22.10.0

## Overview

This pipeline implements state-of-the-art pangenome-based variant calling using the Human Pangenome Reference Consortium (HPRC) graph. It combines vg toolkit for graph alignment, pangenome-aware DeepVariant for SNP/indel calling, and PanGenie for structural variant genotyping.

## Key Features

### Multi-Sample Processing
- **CSV input format**: Process multiple samples in parallel
- **Automatic sample tracking**: Each sample tracked through entire pipeline
- **Independent processing**: Samples processed independently after graph preparation

### Container-Based Reproducibility
- **Singularity containers**: All software in pre-built containers
- **Version locked**: Specific versions for reproducibility
- **Cache management**: Configurable cache directory to preserve disk space

### HPC Optimization
- **SAGA/TSD profiles**: Pre-configured for Norwegian HPC systems
- **SLURM integration**: Automatic job submission and management
- **Resource tuning**: Adjustable CPU/memory per process

## Architecture

### Pipeline Structure

```
nf_pangenome/
├── main.nf                      # Main workflow orchestration
├── nextflow.config              # Configuration and profiles
├── modules/                     # Process definitions
│   ├── download_graph.nf        # Download HPRC graph
│   ├── augment_graph.nf         # Add population SVs (optional)
│   ├── index_graph.nf           # Build graph indexes
│   ├── giraffe_align.nf         # Align reads to graph
│   ├── project_bam.nf           # Convert GAM to BAM
│   ├── call_snps_indels.nf      # DeepVariant variant calling
│   ├── call_svs.nf              # PanGenie SV genotyping
│   └── postprocess.nf           # Merge and filter variants
└── bin/                         # Helper scripts
    ├── download_test_data.sh    # Download test data
    ├── prepare_test.sh          # Setup test environment
    ├── run_test_saga.sh         # Test on SAGA
    └── run_test_tsd.sh          # Test on TSD
```

### Workflow Details

#### Phase 1: Graph Preparation (One-time)
```
INPUT: HPRC v1.1 URL or local file
  ↓
DOWNLOAD_GRAPH (if not provided)
  ↓ hprc-v1.1-mc-grch38.gbz (~50GB)
  ↓
AUGMENT_GRAPH (optional)
  ↓ Add 1KGP SVs to graph
  ↓
INDEX_GRAPH
  ↓ Build distance & minimizer indexes
  ↓
OUTPUT: Indexed graph ready for alignment
```

#### Phase 2: Per-Sample Processing (Parallel)
```
INPUT: sample_name, R1.fq.gz, R2.fq.gz
  ↓
GIRAFFE_ALIGN
  ├─ Align reads to pangenome graph
  ├─ Uses graph-aware alignment algorithm
  └─ Output: sample.gam (graph alignment)
  ↓
PROJECT_BAM
  ├─ Project graph alignment to linear reference
  ├─ Convert GAM → BAM (hg38 coordinates)
  ├─ Sort and index BAM
  └─ Output: sample.bam + sample.bam.bai
  ↓
┌─────────────┴─────────────┐
│                           │
CALL_SNPS_INDELS       CALL_SVS
│                           │
├─ DeepVariant         ├─ PanGenie
├─ Uses graph          ├─ Genotype known
├─  context               SVs from graph
├─ Calls SNPs          └─ Output:
│   and indels            sample.pg.vcf.gz
└─ Output:
   sample.dv.vcf.gz
   sample.dv.g.vcf.gz
│                           │
└─────────────┬─────────────┘
              ↓
       POSTPROCESS
         ├─ Merge SNPs/indels + SVs
         ├─ Filter and normalize
         ├─ Sort and index
         └─ Output:
            sample.merged.vcf.gz
            sample.snps_indels.vcf.gz
            sample.svs.vcf.gz
```

## Software Components

### Core Tools

| Tool | Version | Container | Function |
|------|---------|-----------|----------|
| vg | v1.51.0 | quay.io/vgteam/vg:v1.51.0 | Graph construction, indexing, alignment |
| DeepVariant | pangenome-aware | gcr.io/deepvariant-docker/deepvariant:pangenome_aware_deepvariant-head737001992 | SNP/indel calling with graph context |
| PanGenie | v3.0.0 | quay.io/biocontainers/pangenie:3.0.0--h7d875b9_0 | SV genotyping from graph |
| samtools | v1.18 | quay.io/biocontainers/samtools:1.18--h50ea8bc_1 | BAM manipulation |
| bcftools | v1.18 | quay.io/biocontainers/bcftools:1.18--h8b25389_0 | VCF manipulation |

### Tool Descriptions

**vg (variation graph toolkit)**
- Graph-based reference genome representation
- Includes Giraffe aligner for fast graph alignment
- Handles complex genomic variation natively

**Pangenome-aware DeepVariant**
- Extension of Google's DeepVariant
- Uses graph context to improve variant calling
- Reduces reference bias in complex regions

**PanGenie**
- Genotypes structural variants from pangenome graphs
- Uses k-mer counting for efficient genotyping
- Handles large SVs (>50bp) effectively

## Data Requirements

### Input Data

#### Samples CSV Format
```csv
sample,read1,read2
patient001,/full/path/to/patient001_R1.fastq.gz,/full/path/to/patient001_R2.fastq.gz
patient002,/full/path/to/patient002_R1.fastq.gz,/full/path/to/patient002_R2.fastq.gz
```

**Requirements:**
- Header row required: `sample,read1,read2`
- Sample names must be unique
- Use absolute paths for FASTQ files
- FASTQ files must be gzip compressed

#### Reference Genome
- **Required**: GRCh38/hg38 (no ALT contigs)
- **Format**: FASTA (uncompressed or gzip)
- **Recommended**: NCBI analysis set

#### HPRC Pangenome Graph
- **Version**: HPRC v1.1 (freeze1)
- **Format**: GBZ (compressed graph format)
- **Size**: ~50 GB
- **Auto-download**: Yes (if not provided)
- **URL**: https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz

### Sequencing Requirements

| Parameter | Minimum | Recommended | Maximum |
|-----------|---------|-------------|---------|
| Coverage | 20x | 30x | 100x |
| Read Length | 100bp | 150bp | 250bp |
| Read Type | Paired-end | Paired-end | Paired-end |
| Insert Size | 300bp | 350-500bp | 1000bp |
| Quality | Q20 | Q30 | Q40 |

## Resource Requirements

### Compute Resources

#### Graph Preparation (One-time per graph)
| Process | CPUs | Memory | Time | Disk |
|---------|------|--------|------|------|
| Download | 1 | 4 GB | 1-2h | 50 GB |
| Augment | 8 | 32 GB | 2-4h | 60 GB |
| Index | 16 | 64 GB | 2-4h | 70 GB |

#### Per-Sample Processing
| Process | CPUs | Memory | Time | Disk |
|---------|------|--------|------|------|
| Alignment | 16 | 64 GB | 6-8h | 50 GB |
| Project BAM | 8 | 32 GB | 2-4h | 30 GB |
| DeepVariant | 16 | 64 GB | 12-16h | 20 GB |
| PanGenie | 8 | 32 GB | 4-6h | 15 GB |
| Postprocess | 4 | 16 GB | 1-2h | 5 GB |

**Total per sample (30x WGS):** ~20-24 hours, ~64 GB peak memory, ~50-100 GB temporary disk

### Storage Requirements

| Component | Size | Notes |
|-----------|------|-------|
| HPRC graph | 50 GB | One-time download |
| Graph indexes | 10-20 GB | Built once per graph |
| Container images | 5-10 GB | Cached in singularity_cache/ |
| Per-sample temp | 50-100 GB | Deleted after completion |
| Per-sample output | 2-5 GB | Final VCF files |

**Planning Example:**
- 10 samples, 30x coverage each
- Graph + indexes: 70 GB
- Containers: 10 GB
- Temporary (parallel): 200-500 GB (4-5 samples running simultaneously)
- Final outputs: 20-50 GB
- **Total needed: ~300-600 GB**

## HPC Configuration

### SAGA Profile
```groovy
profiles {
    saga {
        executor {
            name = 'slurm'
            account = 'nn9114k'
            queueSize = 50
        }
        process {
            cpus = 16
            memory = 32.GB
            time = 12.h
            queue = 'normal'
        }
    }
}
```

### TSD Profile
```groovy
profiles {
    tsd {
        executor {
            name = 'slurm'
            account = 'p33_norment'
            perCpuMemAllocation = true
            queueSize = 50
        }
        process {
            cpus = 16
            memory = 32.GB
            time = 12.h
            queue = 'normal'
        }
    }
}
```

## Performance Characteristics

### Scalability
- **Sample parallelization**: Unlimited (resource-dependent)
- **Graph preparation**: One-time, not parallelized across samples
- **Alignment**: Per-sample, highly parallel within sample
- **Variant calling**: Per-sample, can run all samples simultaneously

### Bottlenecks
1. **Graph download**: Network speed (1-2 hours)
2. **Graph indexing**: Memory-intensive (64 GB)
3. **DeepVariant**: Most time-consuming step (12-16h per sample)

### Optimization Strategies
1. **Pre-download graph**: Save 1-2 hours on first run
2. **Use local graph**: Avoid repeated downloads
3. **Increase parallelization**: Run multiple samples simultaneously
4. **Resource tuning**: Adjust CPUs/memory based on HPC availability

## Output Details

### VCF Files Per Sample

#### merged.vcf.gz
- All variant types combined
- SNPs, indels, and SVs
- Sorted and indexed
- Ready for downstream analysis

#### snps_indels.vcf.gz
- Only SNPs and small indels (<50bp)
- From DeepVariant
- High-quality calls
- Suitable for GWAS, QC

#### svs.vcf.gz
- Only structural variants (≥50bp)
- From PanGenie
- Graph-based genotyping
- Includes deletions, insertions, inversions

### VCF Format Details

**INFO fields:**
- `SVTYPE`: Variant type (SNP, INS, DEL, etc.)
- `SVLEN`: SV length (for SVs)
- `AF`: Allele frequency (if available)

**FORMAT fields:**
- `GT`: Genotype (0/0, 0/1, 1/1)
- `GQ`: Genotype quality
- `DP`: Read depth
- `AD`: Allelic depths
- `PL`: Phred-scaled likelihoods

## Testing

### Test Data
Pipeline includes test data download script that fetches:
- GRCh38 chr22 reference (~50 MB)
- HPRC v1.1 graph (~50 GB)
- NA12878 reads: 200M and 75M coverage from Zenodo
- Total download: ~55 GB

### Running Tests

```bash
# Step 1: Download test data
bash bin/download_test_data.sh

# Step 2: Prepare environment
bash bin/prepare_test.sh

# Step 3: Run test
bash bin/run_test_saga.sh  # or run_test_tsd.sh
```

### Expected Test Results
- Runtime: ~6-8 hours (chr22 only)
- Output VCFs: ~100-200 MB per sample
- Variant counts: ~50,000-100,000 variants per sample

## Troubleshooting Guide

### Common Issues

#### 1. "No such file or directory" errors
**Cause:** Relative paths in samples.csv
**Solution:** Use absolute paths
```bash
# Convert to absolute
realpath data/sample1_R1.fastq.gz
```

#### 2. Out of memory errors
**Cause:** Insufficient memory allocation
**Solution:** Increase memory in nextflow.config
```groovy
process {
    withLabel: 'deepvariant' {
        memory = '128 GB'
    }
}
```

#### 3. Job timeout
**Cause:** Time limit too short
**Solution:** Increase time in nextflow.config
```groovy
process {
    withLabel: 'deepvariant' {
        time = '48h'
    }
}
```

#### 4. Graph download failures
**Cause:** Network issues
**Solution:** Download manually
```bash
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz
```

#### 5. Singularity cache filling $HOME
**Cause:** Default cache location
**Solution:** Set custom cache directory
```bash
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache
mkdir -p singularity_cache
```

## Best Practices

### Data Preparation
1. Use high-quality sequencing data (Q30+)
2. Ensure adequate coverage (30x recommended)
3. Use absolute paths in sample sheets
4. Validate FASTQ files before running

### Resource Management
1. Pre-download HPRC graph for production use
2. Set Singularity cache to large filesystem
3. Monitor disk space during runs
4. Use `-resume` for interrupted runs

### Quality Control
1. Check variant counts per sample
2. Verify genotype quality distributions
3. Compare with linear reference results
4. Validate known variants in controls

### Production Deployment
1. Test with small dataset first
2. Adjust resources based on HPC availability
3. Plan disk space for multiple samples
4. Set up monitoring and logging
5. Document run parameters

## Version History

### v1.0 (December 2025)
- Initial release
- Multi-sample CSV input support
- HPRC v1.1 graph integration
- SAGA and TSD profiles
- Singularity container support
- Automated testing scripts

## Future Enhancements

### Planned Features
- [ ] Region-specific variant calling
- [ ] Multi-threading optimization
- [ ] Joint genotyping across samples
- [ ] Phasing support
- [ ] Long-read support (PacBio/ONT)
- [ ] Quality control metrics
- [ ] HTML report generation

### Under Consideration
- [ ] Cloud deployment support
- [ ] Docker container option
- [ ] Custom graph construction
- [ ] Population-specific graphs
- [ ] Integration with annotation tools

## Support and Contribution

### Getting Help
- GitHub Issues: Report bugs and request features
- Documentation: Refer to README.md
- SAGA Support: support@metacenter.no
- TSD Support: tsd-drift@usit.uio.no

### Contributing
Contributions welcome! Please:
1. Fork repository
2. Create feature branch
3. Test changes thoroughly
4. Submit pull request with description

## License
MIT License - See LICENSE file for details

## Acknowledgments
- HPRC Consortium for pangenome resources
- vg team for toolkit and support
- Google DeepVariant team
- PanGenie developers
- Norwegian HPC infrastructure (SAGA/TSD)
