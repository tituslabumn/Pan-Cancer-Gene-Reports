#!/bin/bash

#this script feeds .csv files generated in R into ensemblVEP's PERL servery query API
#run this script in the wd within which the output is to be written(name specified as second arg (char string))

#first argument is the .csv file with the following columns and example content
#identifyer can either be rs ID (dbSNP data) or sample_AAchange (cBioPortal data)
	#	chromosome	start	end	allele	strand identifyer
	#	5		16010	16010	A/T	+	ID####

#third argument is the path of the directory containing the VEP hs_GRCh38 cache and installed plugins

echo "#########################################################"
echo "               Preparing ensemblVEP query"
echo "#########################################################"
echo ""

echo "Started at: $(date +"%T")"
echo ""

echo "Input file:"
echo $1
echo ""

echo "cache/plugin directory"
echo $3
echo ""

#set the path to the vep script below accordingly
#many of the below arguments require --cache and downloaded caches when running INSTAL.pl
#currently querying at USA east coast host server (should be faster)
#$2 is char string of output file name assigned in current wd
#use the following link to see documentation for vep script arguments
	# https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html
~/ensembl-vep/./vep -i $1 --cache --dir $3 --host useastdb.ensembl.org --tab --output_file $2 --hgvs --canonical --domains --check_existing --max_af --af --af_esp --af_gnomad --gencode_basic --force_overwrite --plugin Downstream

#when output as tab info and metadata lines are prefixed with "##" and colnames row is prefexed with "#"
#change for R's default read.table() comment.char argument of "#" while preserving colnames row
echo "Output written to "$2""
echo ""
echo "Changing comment.chars for parsing into R"

sed -i 's/##/~~/' $2
sed -i 's/#//' $2
sed -i 's/~~/#/' $2

echo ""

echo "Finished at: $(date +"%T")"

#example run command for MYO10 directory
#~/myosin_pancancer_project/code/./ensembl_VEP_query.sh ./MYO10_biomaRt_SNPs_filtered_VEP_IN.txt "MYO10_biomaRt_SNPs_filtered_VEP_OUT" /mnt/DATA/.vep
