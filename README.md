# Pan Cancer Gene Reports (PCGR)

## About 
PCGR is an automated informatic pipeline designed by Taylor Harding at the University of Minnesota that creates a pan-cancer analytical report for a user-specified gene of interest (GOI). PCGR is primarily written in the R statistical programming language.  

PCGR sequesters GOI-specific data from cBioPortal (https://www.cbioportal.org/). These data include mutation, copy number, expression and clinical information from >200 studies representing >47,000 unique cases of cancer accross 32 broad tissue types. 

PCGR obtains all natural variants for a GOI observed within the >60,000 exomes analized by the Exome Aggregation Consortium project (http://exac.broadinstitute.org/) 

PCGR additionally obtains information from several database APIs (i.e. Ensembl, UniProtKB) that are used to structurally annotate GOI transcripts and peptides.

See 'example_output.pdf' for an example report.  

## How to use PCGR
PCGR is currently designed to be run in a Linux environment with the R statistical programming language installed. PCGR has not been tested on versions of R <3.6.0  

cBioPortal is continually uloading new studies and the associated processed data. For this reason, PCGR runs in two steps:  
  1. Initial cBioPOrtal query that obtains metasata of all studies and cases to be used downstream. 
      - An initialized R workspace is created that can be used for later GOI queries.
      - Tissue type annotation must be manually specified by the user by editing the cancer_type_manual_annotation_in.tab file
        - If this annotation file is not up to date a template output file will be created and the user will be prompted to update the tissue type annotation
  2. GOI query that sequesters GOI-specific data, performs analysis and generates a PDF report.  

The first step is optional and should only be run if the user wishes to include new data uploaded to cBioPortal that is not included in the initialized workspace provided in this repository.  

### Running PCGR for the first time

PCGR has been developed to run on Debian-style linux distributions (i.e. 64-bit Ubuntu > 18.0). Functionality on Windows/OSX has not been fully tested.

Ensure that you have the following installed/set on your system:
* latest version of R (tested on versions > 3.6.x) 
* ensure that .libPaths are set to a writable directory prior to running PCGR if using R for the first time.
* May require installation of certain parser libraries: XML, Curl, openssl. 
  - On Ubuntu platforms: use "sudo apt-get install" to install the following libraries:
      - libxml2-dev
      - libcurl4-openssl-dev
      - openssl-dev
    
PCGR will install required R package dependencies (this may take a significant amount of time upon first run). Major packages/dependancies include:
- RCurl
- XML
- cgdsr
- BiocManager
- rtracklayer
- biomaRt
- trackviewer
- BSgenome.Hsapiens.NCBI.GRCh38
- rmarkdown
- knittr
- rjson
- Biostrings

#### Running from bash terminal

1. Open a bash terminal and clone the PCGR repository:

```{bash eval=FALSE}
  git clone https://github.com/tituslabumn/Pan-Cancer-Gene-Reports.git
```

2. To run an PCGR:
    - Set your working directory to the cloned Pan-Cancer-Gene-Reports directory
    - PCGR is executed by running the Single_gene_cBioPortal_query.sh bash script
      - The first argument (y/n) specifies whether or not to create a new initialized workspace (study/case data)
      - The second argument is the gene symbol for the GOI
      
```{bash eval=FALSE}
  cd Pan-Cancer-Gene-Reports
  sudo ./Single_gene_cBioPortal_query.sh "y" "NRAS"
```

3. If the 'cancer_type_manual_annotation_in.tab' file in the cloned repository is not up to date with the current database:
    - PCGR will halt and export a cancer_type_manual_annotation_out.tab file for manual annotation
    - Manually annotate each NA populated row with the diesired tissue type under the final_tissue column. 
    - Save the updated file as 'cancer_type_manual_annotation_in.tab' in the same directory
    - Re-run PCGR
    
4. PCGR will create an output directory with sub directories for each GOI querried
    - When finished PCGR will output the PDF report into its respective output directory

If initialization fails to execute to completion use the following instuctons to run PCGR with the latest working initialized workspace provided in this repository:
    
### Running additional queries without re-initializing

1. Run PCGR with an "n" passed to the first argument of the run script

```{bash eval=FALSE}
  ./Single_gene_cBioPortal_query.sh "n" "KRAS"
```




