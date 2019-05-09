#!/bin/bash

#Run full cBioPortal analysis pipeline

#first arg is whether to initialize ("y" or "n")
#second arg is GOI
echo "Started at: $(date +"%T")"
Rscript cBP_single_gene_query.R $1 $2

#cd ./output/$2
#chmod +x VEP_run.sh
#./VEP_run.sh

#cd ..
#cd ..

#Rscript cBP_post_VEP.R $2

echo "Finished at: $(date +"%T")"

