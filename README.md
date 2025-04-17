Overview: 
This repository contains the code and data necessary to run a Gene Set Enrichment Analysis (GSEA) on a curated dataset. 

What this code does:
This code uses a pre-generated set of lists of genes that have been identified as up or downregulated in breast cancer cells (as compared to immune cells penetrating these tumors), and analyzes if these genes are enriched in a different dataset containing tumor and normal samples. This provides a method to validate that these gene lists correspond to "genetic signatures" that can be used to differentiate between the three breast cancer subtypes and thereby inform future treatment strategies.

Source of Data:
We utilized the Gene Expression Omnibus (GEO) Datasets GSE114747 (Dataset A, comprising single-cell RNA sequencing data of tumor and immune cells in the breast cancer microenvironment), and GSE75688 (Dataset B, comprising single-cell RNA sequencing data of primary breast cancer cells and lymph node metastases, segregated by breast cancer subtype). We ran a differential gene expression (DEG) analysis on Dataset B, extracting the top 100 genes that were upregulated and downregulated in tumor cells, compared to immune cells. We did this for each of the three breast cancer subtypes (ER+, HER2+, TNBC) to obtain six CSV files of the top 100 differentially expressed genes in each case. We also merged all patient samples for the "tumor" and "normal" cases for Dataset A to provide a single averaged expression CSV file. 

What this repo contains:
- The six CSV files containing the top 100 differentially expressed genes for each of the cases listed above, named accordingly (from Dataset B).
- A tumor vs normal merged expression CSV file (from Dataset A). 
- A python script to run the GSEA analysis, using the six CSV files (from Dataset B) as the gene sets and the tumor vs normal CSV file (from Dataset A) as the dataset to analyze for enrichment.

Installation:
Simply download this repo's contents and open the python script in any IDE. Then right-click on the downloaded folder's path and copy it to the data_folder variable in line 6. Then you can run the script and it will generate all twelve GSEA plots for these variables. You can modify font size and layout parameters to improve readability of this figure as well.
