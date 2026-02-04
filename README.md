# Obesity-irAEs
Codes for analyzing and visualizing data on how obesity fuels irAEs induced by anti-PD-1 therapy
# This repository contains analysis codes used for the publication: Obesity fuels immune-related adverse events of anti-PD-1 therapy
Codes provided here are to perform DEGs analysis, Pathway Enrichment, Transcription Factor Analysis, and Cell Communication analysis for single cell RNA-seq of adipose and colon tissues,and Public PBMC and BALF dataset, as well as code for DEGs analysis and Pathway Enrichment for BMDM.
1. Single_Cell_RNA_Seq_Data_Analysis_Pipeline: The code in R to perform single cell data analysis, including cell annotation, DEGs analysis, and Pathway Enrichment result visualization for single cell RNA-seq of adipose and colon tissues.
2. Cell_Communication_Analysis: The code in R to perform Cell-Cell communication analysis.
3. Pyscenic_Analysis: The code in Linux shell to perform Pyscenic transcription factor analysis.
4. Public_Data_Analysis: The code in R to perform single cell data analysis of Public PBMC and BALF dataset.
5. Pyscenic_results_visualization: The code in R to find differentially regulated TFs and visualize results of transcription factor analysis.
6. Bulk_RNA_Seq_Data_Analysis: The code in R to perform DEGs analysis and pathway enrichment result visualization of BMDM Bulk RNA seq data.
# Downloading the data
- The raw data for single-cell RNA-seq has been deposited in GSA (https://ngdc.cncb.ac.cn/gsa/) under CRA030409.
- The processed data for single-cell RNA-seq has been deposited in OMIX (https://ngdc.cncb.ac.cn/omix/release/OMIX012013) under OMIX012013.
- The raw data for Bulk RNA-seq of BMDM has been deposited in GSA (https://ngdc.cncb.ac.cn/gsa/) under CRA030408.
- The processed data for Bulk RNA-seq of BMDM has been deposited in OMIX (https://ngdc.cncb.ac.cn/omix/release/OMIX012014) under OMIX012014.
- Public Human PBMC dataset can be downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/) under accession GSE285888.
- Public Human BALF dataset can be downloaded from EGA (https://ega-archive.org/) under accession EGAS00001006762.
# Operation systems
- Linux version 3.10.0-1062.el7.x86_64
- Windows 10 64-bit
# Hardware requirements
The pipeline requires only a standard computer with enough RAM to support the operations. For optimal performance, we recommended and used a computer with the following specs for testing:
Intel® Core™ Ultra 7 Processor 265K and RAM 128 GB.
# Installation guide
R packages required for the pipeline can be installed from CRAN (https://cran.r-project.org/) using the install.packages() function, or from Bioconductor (https://bioconductor.org/) using the BiocManager::install() function.
pySCENIC used for transcription factor analysis can be installed from https://github.com/aertslab/pySCENIC.
# Packages
- python (version 3.8.20)
- pyscenic (version 0.12.1)
- R packages:
- Seurat (version 5.3.0)
- DropletUtils (version 1.26.0)
- DoubletFinder (version 2.0.6)
- ggsci (version 3.2.0)
- RColorBrewer (version 1.1.3)
- ggplot2 (version 3.5.2)
- dplyr (version 1.1.4)
- viridis (version 0.6.5)
- BiocParallel (version 1.40.2)
- harmony (version 1.2.3)
- patchwork (version 1.3.0)
- clustertree (version 0.5.1)
- tidyverse (version 2.0.0)
- paletteer (version 1.6.0)
- hdf5r (version 1.3.12)
- future (version 1.49.0)
- presto (version 1.0.0)
- SCP (version 0.5.6)
- ggpubr (version 0.6.0)
- SeuratDisk (version 0.0.0.9021)
- ggrepel (version 0.9.6)
- AUCell (version 1.28.0)
- GSEABase (version 1.68.0)
- msigdbr (version 10.0.2)
- KEGGREST (version 1.46.0)
- ggSCvis (version 0.0.3)
- CellChat (version 2.1.2)
- edgeR (version 4.4.2)
- 
# Demo 
1. Single_Cell_RNA_Seq_Data_Analysis_Pipeline:
- Input: required input data files, which can be downloaded from OMIX012013
- Output: expected output results
2. Bulk_RNA_Seq_Data_Analysis:
- Input: required input data files, which can be downloaded from OMIX012014
- Output: expected output results
3. Users can also directly provide their own input files with the format required, and perform analysis by modifying our codes
