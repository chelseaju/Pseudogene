#!/bin/bash

## bash analysis_v4.sh DIR COVERAGE READLEN ABUNDANCE
## ex bash analysis_v4.sh select_one_pseudogene_110 1 100 10

DIR=$1
COVERAGE=$2
READLEN=$3
ABUNDANCE=$4
TOPOUT="tophat_out"
SUBDIR=$COVERAGE"X_"$READLEN"L_"$ABUNDANCE"A"
BAM="accepted_hits.bam"
CHROMOSOME=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8'
	'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16'
	'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' 'chrY' 'chrM')

#ENST2ESG="/home/chelseaju/Database/EnsemblTranscriptome/Genes/ENST2ENSG.py"
#ENSG2PGOHUM="/home/chelseaju/Database/EnsemblTranscriptome/Genes/ENSG2PGOHUM.py"

ENST2ENSG="/u/scratch/c/chelseaj/database/EnsemblTranscriptome/Genes/ENST2ENSG.py"
ENSG2PGOHUM="/u/scratch/c/chelseaj/database/EnsemblTranscriptome/Genes/ENSG2PGOHUM.py"

#CHROMOSOME=('chr1' 'chr2')

echo "Data Analysis Version 4:"

echo "Step 0: Folders and Files Preparation"
CMD0="mkdir -p $DIR/$SUBDIR/$TOPOUT/mapping"
CMD1="python bamfile_sorter.py -i $BAM -d $DIR/$SUBDIR/$TOPOUT"
echo ">$CMD0"
echo ">$CMD1"
$CMD0
$CMD1

# STEP 1: Read Counts - a series of scripts
## Compute the expected number of read for each gene
## Map each read to known gene or identify the mapped region
## Compute the observed number of read for each region
## Compute the distribution of read generated from parent genes

echo ""
echo "Step 1: Read Counts for Genes"

# expect read count
python expected_counter.py -d $DIR/$SUBDIR/$TOPOUT
python $ENST2ENSG -i expected_read_count.txt -o ENSG_expected_read_count.txt -d $DIR/$SUBDIR/$TOPOUT

for chr in "${CHROMOSOME[@]}"
do
  echo ""
  echo "CHROMOSOME ${chr} =========================================="
  
  python exon_identifier_v2.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr}
  python gene_identifier_v2.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr}
  python observed_separator.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr} -t genes
  python $ENSG2PGOHUM -d $DIR/$SUBDIR/$TOPOUT/mapping -i ${chr}_genes_distribution.txt -o ${chr}_PGOHUM_distribution.txt
  python $ENST2ENSG -d $DIR/$SUBDIR/$TOPOUT/mapping -i ${chr}_PGOHUM_distribution.txt -o ${chr}_ENSG_distribution.txt

done

# observed read distribution
python distribution_equation.py -d $DIR/$SUBDIR/$TOPOUT -t ENSG
python distribution_matrix.py -d $DIR/$SUBDIR/$TOPOUT -t ENSG



echo ""
echo "COMPLETE !! "

