# Pan Cancer Gene Reports (PCGR)

## About 
PCGR is an automated informatic pipeline designed by Taylor Harding at the University of Minnesota that creates a pan-cancer analytical report for a user-specified gene of interest (GOI).  
PCGR sequesters GOI-specific data from cBioPortal (https://www.cbioportal.org/). These data include mutation, copy number, expression and clinical information from >200 studies representing >47,000 unique cases of cancer accross 32 broad tissue types. 
PCGR obtains all natural variants for a GOI observed within the >60,000 exomes analized by the Exome Aggregation Consortium project (http://exac.broadinstitute.org/) 
PCGR additionally obtains information from several database APIs (i.e. Ensembl, UniProtKB) that are used to structurally annotate GOI transcripts and peptides.

## How to use PCGR
PCGR is currently designed to be run in a Linux environment with the R statistical programming language installed. PCGR has not been tested on versions of R <3.6.0  

### Running PCGR for the first time
cBioPortal is continually uloading new studies and the associated processed data. For this reason, PCGR runs in two steps:  

        a. Initial cBioPOrtal query that obtains metasata of all studies and cases to be used downstream. An initialized R workspace is created that can be used for later GOI queries.
        b. GOI query that sequesters GOI-specific data, performs analysis and generates a PDF report.  

The first step is optional and should only be run if the user wishes to include new data uploaded to cBioPortal that is not included in the initialized workspace provided in this repository.  

1. Seven utilizes Fiji/ImageJ application. If you already don't have the latest version, Please download it from [here](https://fiji.sc/).

2. Install the updated version of seven from Github repository and move all the seven files to the Fiji.app root folder.

3. Move the tiff files from an experiment to a folder with the name of the file string beginning i.e. 95.1 - 95.25 into a folder 95 or so.

4.	From Fiji/Image J Select :

        a.  Plugins>Macros>Edit>seven_multi
        b.  Run seven_multi file.
        c.  Select the folder with tiff files which you want to analyze. 

5.	Seven Runs through the folder and creates a ZIP mask for each image.

6. The seven will ask you to Manually edit the image file incase the tiff image is not clear.

7. Manually erase any erroneous cells by comparing with the original tiff file, if image is not clearly inverted.

8. For each image in the experiment folder, seven will create a sliced folder ( for ex. sliced_1-1 or so). These folder will contained all the stats exported from an image which will be helpful in conducting analysis.



