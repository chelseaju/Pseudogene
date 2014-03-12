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

## MAC
#ENST_ENSG_ENSP="/home/chelseaju/Database/Ensembl/ENST_ENSG_ENSP_74.txt"
#ENST2ENSG="/home/chelseaju/Database/Ensembl/script/ENST2ENSG.py"

## HOFFMAN
ENST_ENSG_ENSP="/u/home/c/chelseaj/project/database/Ensembl/ENST_ENSG_ENSP_74.txt"
ENST2ENSG="/u/home/c/chelseaj/project/database/Ensembl/script/ENST2ENSG.py"

## LAB
#ENST_ENSG_ENSP="/home/chelseaju/Database/Ensembl/ENST_ENSG_ENSP_74.txt"
#ENST2ENSG="/home/chelseaju/Database/Ensembl/script/ENST2ENSG.py"

#CHROMOSOME=('chr1' 'chr2')

echo "Data Analysis Version 4:"

echo ""
echo "Step 1 - Creating Folder"
mkdir -p $DIR/$SUBDIR/$TOPOUT/mapping

echo ""
echo "Step 2 - Converting BAM to SAM"
python bamfile_sorter.py -i $BAM -d $DIR/$SUBDIR/$TOPOUT


# STEP 1: Read Counts - a series of scripts
## Compute the expected number of read for each gene
## Map each read to known gene or identify the mapped region
## Compute the observed number of read for each region
## Compute the distribution of read generated from parent genes

echo ""
echo "Step 1: Read Counts for Genes"

# expect read count
python expected_counter.py -d $DIR/$SUBDIR/$TOPOUT
python $ENST2ENSG -i $DIR/$SUBDIR/$TOPOUT/expected_read_count.txt -o $DIR/$SUBDIR/$TOPOUT/ENSG_expected_read_count.txt -d $ENST_ENSG_ENSP

# merge the same ENSG
cp $DIR/$SUBDIR/$TOPOUT/ENSG_expected_read_count.txt $DIR/$SUBDIR/$TOPOUT/ENSG_expected_read_count.backup
awk '{A[$1]+=$2; next} END{for (i in A) {print i,"\t",A[i]}}' $DIR/$SUBDIR/$TOPOUT/ENSG_expected_read_count.backup $DIR/$SUBDIR/$TOPOUT/ENSG_expected_read_count.txt

for chr in "${CHROMOSOME[@]}"
do
  echo ""
  echo "CHROMOSOME ${chr} =========================================="
  
  python exon_identifier_v2.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr}
  python gene_identifier_v2.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr}
  python observed_separator.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr} -t genes
  python $ENST2ENSG -d $ENST_ENSG_ENSP -i $DIR/$SUBDIR/$TOPOUT/mapping/${chr}_genes_distribution.txt -o $DIR/$SUBDIR/$TOPOUT/mapping/${chr}_ENSG_distribution.txt

done

# observed read distribution
python distribution_equation.py -d $DIR/$SUBDIR/$TOPOUT -t ENSG
python distribution_matrix.py -d $DIR/$SUBDIR/$TOPOUT -t ENSG



echo ""
echo "COMPLETE !! "

