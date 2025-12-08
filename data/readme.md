# Data Directory

This directory contains reference data for the PanGenie genotyping pipeline.

## 1000 Genomes Project ONT Phased VCF

### File
`1000GP_ONT_shapeit5-phased-callset_final-vcf.phased.vcf.gz/`

### Description
Phased VCF from the 1000 Genomes Project long-read (Oxford Nanopore) dataset, phased with SHAPEIT5.

### Citation
Ebert, P., Audano, P.A., Zhu, Q. et al. Haplotype-resolved diverse human genomes and integrated analysis of structural variation. *Nature* **637**, 696–706 (2025). https://doi.org/10.1038/s41586-025-09290-7

## HPRC Pangenome Graph

### Download
If you prefer to extract variants from the HPRC pangenome graph:

```bash
wget -O data/hprc-v1.1-mc-grch38.gbz \
    https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz
```

### Citation
Liao, W.W., Asri, M., Ebler, J. et al. A draft human pangenome reference. *Nature* **617**, 312–324 (2023). https://doi.org/10.1038/s41586-023-05896-x

## Usage

### Option 1: Use 1000GP ONT Phased VCF
```bash
nextflow run main.nf \
    --input_csv samples.csv \
    --reference reference.fa \
    --phased_vcf data/1000GP_ONT_shapeit5-phased-callset_final-vcf.phased.vcf.gz \
    -profile saga
```

### Option 2: Extract VCF from HPRC Graph
```bash
nextflow run main.nf \
    --input_csv samples.csv \
    --reference reference.fa \
    --hprc_graph data/hprc-v1.1-mc-grch38.gbz \
    --use_graph_vcf true \
    -profile saga
```
