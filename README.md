# Pan-Cancer Gene Reports (PCGR) v1.0

## About 
PCGR is an automated informatic pipeline designed by Taylor Harding at the University of Minnesota that creates a pan-cancer analytical report for a user-specified gene of interest (GOI). PCGR is primarily written in the R statistical programming language.  

PCGR sequesters GOI-specific data from cBioPortal (https://www.cbioportal.org/). These data include mutation, copy number, expression and clinical information from >200 studies representing >47,000 unique cases of cancer accross 32 broad tissue types. 

PCGR obtains all natural variants (excluding those associated with cancer cases) for the supplied GOI contained within the genomes and exomes collected by the gnomAD project (https://gnomad.broadinstitute.org/). 

PCGR additionally obtains information from several database APIs (i.e. Ensembl, UniProtKB) that are used to structurally annotate GOI transcripts and peptides.

See 'example_output.pdf' for an example report (GOI = NRAS).  

PCGR is currently designed to be run in a Debian Linux environment with aditional requirments outlined below. PCGR has not been tested on versions of R <3.6.0  

To clone this source code:

```{bash eval=FALSE}
  git clone https://github.com/tituslabumn/Pan-Cancer-Gene-Reports.git
```

## How to use PCGR
Due to the large number of specific software requirements it is HIGHLY RECCOMENDED that the user run PCGR from within a the PCGR Docker container that has been built to provide an environment that has all the required packages pre-installed. 

Ensure that you have Docker installed on your system and pull the PCGR image from DockerHub (>2GB):

```{bash eval=FALSE}
  docker pull tsharding/pcgr_v1.0
```

Create an active container from the Docker image. 
Upon running, the container will:
* Pull the latest version of PCGR from this repository
* Set the container's working directory to the cloned directory
* Initiate a bash shell through which the user can run a PCGR query for their GOI

```{bash eval=FALSE}
  docker run -ti --rm  -v <path_to_your_host_output_directory>:/OUTPUT tsharding/pcgr_v1.0
```

The above command will run the container interactivly (-ti), delete the container after it closes (--rm) and mount (-v) the container's output to the to the directory chosen by the user (output files will persist after container is closed/removed).

To run a PCGR query in the bash shell of the container execute the run script with the GOI gene symbol as an argument:

```{bash eval=FALSE}
  ./run_PCGR_query.sh NRAS
```

After the query is complete the output PDF and final R workspace image will be saved to the output directory (output will only persist if directory is properly mounted).  

To exit and delete the container:

```{bash eval=FALSE}
  exit
```

### Updating the initialized workspace
PCGR quries run using an R workspace image (provided in this repository) that contains the gene-agnostic data regarding cancer studies, samples, etc. from cBioPortals API. This prevents having to re-sequester and process this data upon each query. cBioPortal is, however, frequently updated. To re-build the initialized workspace image run the associated shell script from within the PCGR container:

```{bash eval=FALSE}
  ./initialize_PCGR.sh
```

This command will need to be run twice. The first run will produce a
manual annotation input and output files do not have equal number of rows:\n\tmanually annotate 'cancer_type_manual_annotation_out.tab' file \n\tadd final_tissue column \n\tsave as 'cancer_type_manual_annotation_in.tab'\n\trun again\n\n\n

### Requirements
If running PCGR in your own environment without using the provided PCGR Docker container ensure that you have installed these requirements.

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

### Known issues
Occasionally one of the several web connections/APIs (e.g. UniProt) will fail to connect for unknown reasons. In these cases you may need to retry your query at another time.



