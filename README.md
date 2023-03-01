# Caerus reanalysis pipeline
Reanalysis pipeline for unsolved 100K EGLH rare disease cases

## Background



## Aims



## Data sources
GEL CIP-API: JSON payload and small variant VCF files for each case
CellBase API: Variant annotation
PanelApp API: Panel information

## Static files

Chromosome notation standardisation
- Chromosome scaffold notation: NCBI website GRCh38.p14 assembly page
- https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt

GRCh37 VCF file liftover
- Genome FASTA file
- b37 to b38 chain file

## Data structures

Code storage structure in GitHub

File storage structure in DNAnexus

## Data flow in program

What is a normal input
How does the program deal with this
What is the output


## Software requirements
PyPI packages described in requirements.txt

Also requires: <VERSIONS>
- bedtools
- bcftools
- vcftools
- CrossMap.py
