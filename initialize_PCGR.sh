#!/bin/bash

#Generate initialized workspace for single-gene queries

#arg is GOI gene symbol
echo "Started at: $(date +"%T")"
Rscript initialize_PCGR.R
echo "Finished at: $(date +"%T")"

