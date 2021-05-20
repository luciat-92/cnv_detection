# CNV detection pipeline from genotyping data
Digital karyotyping pipeline to detect de novo copy number abnormalities arising in cultured cell lines i.e. to determine differences between cell lines and the starting material from which they were derived in term of copy number variation. The genomic screening for chromosomal abnormalities can be used as quality control to establish and maintain stem cell lines.

### Built with
    - GenomeStudio for Genotyping with PLINK Input Report Plug-in v2.1.4
    - R (3.5.1) with packages: 
        - argparse 
        - doParallel 
        - ggplot2 
        - RColorBrewer 
        - ggrepel
        - ComplexHeatmap
        - circlize
        - grid
        - gridExtra
    - Plink1.9
    - BCFtools
    - bedtools
    - VCFtools


