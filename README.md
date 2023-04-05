# Caerus reanalysis pipeline
Reanalysis pipeline for unsolved 100K EGLH rare disease cases

## Background

Rare disease cases have a low diagnostic rate in clinical genomics. This can be due to several reasons including lack of evidence for causal variants in shared genomic knowledgebases or low disease penetrance. However, if we reanalyse an unsolved case at a later time point we may be able to leverage more recent information to reach a diagnosis.

This program takes as input one or more IDs of cases in the 100,000 Genomes Project (100k) recruited via EGLH. It then retrieves the JSON file containing clinical data for that case, as well as any available unfiltered SNV VCFs. The clinical data is combined with various filtering parameters to return an output list of variants which meet those parameters.

Currently, JSON and VCF files are only available for a limited number of cases which were previously solved. The cases which can currently be analysed are listed in available_cases.txt.

## Requirements

Required PyPI packages are described in the requirements.txt file.

The software also requires:

- bedtools v2.30.0
- bcftools v1.15
- vcftools v0.1.17
- CrossMap.py v0.6.4

In order to be able to download the JSON and VCF files to process a case, which are stored in DNAnexus, you must be a member of the EGLH bioinformatics team with DNAnexus access credentials.

## Usage

The program downloads one or more VCF files for each case, which may be relatively large. It also generates numerous intermediate files with a high combined total file size. If you are running the program for multiple cases, it may be advisable to do this via a DNAnexus cloud workstation.

The program can be called using various inputs:

- A single case ID
- A space-separated list of case IDs
- The path to a .txt file containing a case ID on each line
- The argument 'all' which will process all cases for which files are available

Examples of program execution for these are:

```
src/caerus.py <case ID>
src/caerus.py <case ID 1> <case ID 2> ... <case ID N>
src/caerus.py <file listing case IDs>.txt
src/caerus.py all
```

## Outputs

The key outputs, which are created in the root directory, are the Excel output and the sorted versions of each available VCF.

The Excel file has three worksheets:

1. Case info: Case ID, family members with VCFs available, filtering parameter thresholds applied, whether the case was originally solved, whether family filtering was performed in reanalysis, whether liftover to GRCh38 was required
2. Variants: Details of the variants which passed filtering parameters
3. Panels: Details of the PanelApp panels used to create a BED file - their original and current versions, and the number of genes and regions in the current version

The sorted VCFs have been filtered on a BED file and variant call quality metrics (QUAL score, DP/read depth, and GQ/genotype quality).

Additional outputs are created in the intermediate_files folder, which will be created by the program if it does not already exist. These include intermediate VCF files, the BED file used to filter variants on genotype-phenotype associations, and temporary files created during the annotation process.

## Data sources

GEL interpretation portal: Original source of JSON/VCF files. Not accessed via this software - requires access credentials and an internal CUH network connection.

PanelApp API: Used to retrieve data for current panel versions for specified phenotypes.

Ensembl BioMart API (GRCh38): Used to obtain chromosome, start and end positions for each gene in the case's panels in order to create a BED file.

CellBase API: Used to retrieve variant annotations.

CADD tool API: Used to retrieve PHRED-scaled CADD scores for SNVs where CellBase is unable to provide this.

## Static files

header.txt : Defines the header line to add to VCF files if they need to be lifted over from GRCh37 to GRCh38. This line defines a new INFO field tag, 'PCV' (for putative causal variant), which is used to tag variants reported in the original analysis before liftover.

cb_config.json : Defines the configuration to use when accessing the CellBase API.

chrom_map.txt : Defines the mapping between standard and UCSC chromosome notation.

bed_template.txt : Contains the template query for requesting BED file output from BioMart.

hg19ToHg38.over.chain.gz : Chain file for lifting over GRCh37 VCFs to GRCh38.

GCF_000001405.40_GRCh38.p14_genomic.fna : GRCh38 reference genome FASTA file.

GCF_000001405.40_GRCh38.p14_genomic.fna.fai : GRCh38 reference genome FASTA file index.
