# CNV detection pipeline from genotyping data
Digital karyotyping pipeline to detect de novo copy number abnormalities arising in cultured cell lines i.e. to determine differences between cell lines and the starting material from which they were derived in term of copy number variation. The genomic screening for chromosomal abnormalities can be used as quality control to establish and maintain stem cell lines.

## About the project
The analysis is divided in two steps:
    1. Sample hybridization using Illumina Global Screening Array-24 v1.0 [[1]](#1) and Genotyping (GT) module of GenomeStudio software to call the genotypes
        * For each probe, GT module estimates Log R ratio (LRR) and B-allele frequency (BAF) using a clustering module applied to the distribution of signal intensities
    1. Pairwise comparison between cell line (e.g. iPSC) and starting material (e.g. fibroblast) using GT module estimates of LRR and BAF to detect copy number variations (CNVs) through two-sample Hidden Markov Model [[2]](#2) 

## Built with
 * GenomeStudio for Genotyping with PLINK Input Report Plug-in v2.1.4 (only on Windows)
 * R (3.5.1) with packages: 
     * argparse 
     * doParallel 
     * ggplot2 
     * RColorBrewer 
     * ggrepel
     * ComplexHeatmap
     * circlize
     * gridExtra
  * Plink1.9
  * BCFtools
  * bedtools
  * VCFtools

## Usage
For details please refer to [documentation](https://gitlab.mpcdf.mpg.de/luciat/cnv_detection/-/blob/master/pipeline_description/CNV_analysis_pipeline.pdf)

## References
<a id="1">[1]</a> Infinium Global Screening Array-24 v1.0, A powerful, high-quality, cost-effective array for population-scale genetic studies [Technical Data Sheet](http://www.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/infinium-commercial-gsa-data-sheet-370-2016-016.pdf)  

<a id="2">[2]</a> Danecek P., McCarthy S.A., HipSci C. and Durbin  R., A Method for Checking Genomic Integrity in Cultured Cell Lines from SNP Genotyping Data. _PLoS One_ 11 (2016). 

