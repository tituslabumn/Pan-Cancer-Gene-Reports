#!/bin/bash

#Run full PCGR pipeline for gene of interest

#arg is GOI gene symbol
echo "Started at: $(date +"%T")"
Rscript run_PCGR_query.R $1
echo "Finished at: $(date +"%T")"

