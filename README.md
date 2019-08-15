# Pan-Cancer Gene Reports (PCGR) v1.0

## About 
PCGR is an automated informatic pipeline written by Taylor Harding at the University of Minnesota that creates a pan-cancer analytical report for a user-specified gene of interest (GOI). PCGR is primarily written in the R statistical programming language.  

PCGR sequesters GOI-specific data from cBioPortal (https://www.cbioportal.org/). These data include mutation, copy number, expression and clinical information from >200 studies representing >47,000 unique cases of cancer accross 32 broad tissue types. 

PCGR obtains all natural variants (excluding those associated with cancer cases) for the supplied GOI contained within the genomes and exomes collected by the gnomAD project (https://gnomad.broadinstitute.org/). 

PCGR additionally obtains information from several database APIs (i.e. Ensembl, UniProtKB) that are used to structurally annotate GOI transcripts and peptides.

See 'example_output.pdf' for an example report (GOI = NRAS).  

PCGR source code is currently designed to be run in a Debian Linux environment with aditional requirments outlined below. PCGR has not been tested on versions of R <3.6.0  

## How to run PCGR from a Docker container
Due to the large number of specific software requirements it is HIGHLY RECCOMENDED that the user run PCGR from within a the PCGR Docker container that has been built to provide an environment that has all the required packages pre-installed. Because components of the PCGR pipline use RSelenium, the PCGR container must be run alongside a selenium-standalone-chrome container. To streamline execution of these environments 

Ensure that you have:
* Docker (and docker-compose) installed on your system
* 5GB free space available on your system

Clone this repository which includes shell scripts to manage container composition:

```{bash eval=FALSE}
  git clone https://github.com/tituslabumn/Pan-Cancer-Gene-Reports.git
```

Set working directory and execute the container startup shell script with the path to your output directory as an argument. 


```{bash eval=FALSE}
  cd Pan-Cancer-Gene-Reports/RUN_PCGR_CONTAINER
  ./RUN_PCGR_CONTAINER.sh [/path/to/your/output/directory]
```

This will:
* Automatically pull images from dockerhub if needed.
* Start both the PCGR and selenium server containers.
* Clone the most recent version of the PCGR source code to the PCGR container at entrypoint.
* Initiate an interactive bash shell from which PCGR querys are run


To run a PCGR query in the bash shell of the container execute the run script with the GOI gene symbol as an argument:

```{bash eval=FALSE}
  ./run_PCGR_query.sh NRAS
```

After the query is complete the output PDF and final R workspace image will be saved to the user-provided output directory.  

To exit and delete the containers (data in output directory persists):

```{bash eval=FALSE}
  exit
```

### Updating the initialized workspace (optional)
PCGR quries run using an R workspace image (provided in this repository) that contains the gene-agnostic data regarding cancer studies, samples, etc. from cBioPortals API. This prevents having to re-sequester and process this data upon each query. cBioPortal is, however, frequently updated. To re-build the initialized workspace image run the associated shell script from within the PCGR container:

```{bash eval=FALSE}
  ./initialize_PCGR.sh
```

This command will need to be run twice. 
The first run will produce a sorted 'cancer_type_manual_annotation_out.tab' file. Manually inspect and nnotate/complete this file's final_tissue column and save as 'cancer_type_manual_annotation_in.tab' in the same directory. Run the initialization again.

### Requirements
If running PCGR in your own environment without using the provided PCGR Docker container ensure that you have installed these requirements. (the code in this repository was written for use in container-specific directories and may therefore require adjustment to source without error)

NOTE: RSelenium still requires a running selenium-standalone-chrome container.

```{bash eval=FALSE}
  docker run -d -rm -p 4444:4444 -v [/path/to/your/output/directory]:/home/seluser/Downloads selenium/standalone-chrome
```

PCGR has been developed to run on Debian-style linux distributions.

Ensure that you have the following installed/set on your system:
* latest version of R (tested on versions > 3.6.x) 
* ensure that .libPaths are set to a writable directory prior to running PCGR if using R for the first time.
* Requires installation of certain parser libraries (XML, Curl, openssl), latex (TinyTex), etc.:
  - On debian linux platforms: use "apt-get install" to install the following libraries:
      - libxml2-dev
      - libcurl4-openssl-dev
      - libssl-dev
      - pandoc
      - default-jdk
      - git
      - Docker
    
PCGR will install required R package dependencies upon running initialization (this may take a significant amount of time upon first run). 
Major packages/dependancies include:
- BiocManager (version='devel')
- RCurl
- XML
- xml2
- cgdsr
- rtracklayer
- biomaRt
- trackviewer
- BSgenome.Hsapiens.NCBI.GRCh38
- rmarkdown
- knittr
- rjson
- Biostrings
- tinytex
- RSelenium

### Known issues
Occasionally one of the several web connections/APIs (e.g. UniProt) will fail to connect for unknown reasons (ssl error). In these cases you may need to retry your query at another time.



