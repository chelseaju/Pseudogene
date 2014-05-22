#!/bin/bash
## sh evaluation.sh DIR A B BAM SUFFIX

DIR=$1
SUBDIR=$2"X_"$3"L_"$4"A"
LASSODIR=$5
TOPOUT="tophat_out"
LASSO_X="$DIR/$LASSODIR/$SUBDIR""_reconstructed_x_matrix.txt"
LASSO_Y="$DIR/$LASSODIR/$SUBDIR""_validationY_v2.txt"
CUFFLINKS="cufflinks_after"


## MAC
#PSEUDOGENE="/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/Pseudogene_74.bed"
#PARENT="/Users/Chelsea/Bioinformatics/CJDatabase/Ensembl/ENSG_74.bed"
#BED="/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/ParentENSG_Pseudogene_74.bed"
#ENST_ENSG_ENSP="/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/ENST_ENSG_ENSP_74.txt"
#ENST2ENSG="/Users/Chelsea/Bioinformatics/CJDatabase/Ensembl/script/ENST2ENSG.py"

## HOFFMAN
#PSEUDOGENE="/u/home/c/chelseaj/project/database/Pseudogene/Pseudogene_74.bed"
#PARENT="/u/home/c/chelseaj/project/database/Ensembl/ENSG_74.bed"
#BED="/u/home/c/chelseaj/project/database/Pseudogene/ParentENST_Pseudogene_74.bed"
#ENST_ENSG_ENSP="/u/home/c/chelseaj/project/database/Ensembl/ENST_ENSG_ENSP_74.txt"
#ENST2ENSG="/u/home/c/chelseaj/project/database/Ensembl/script/ENST2ENSG.py"

## LAB
PARENT="/home/chelseaju/Database/Ensembl/ENSG_74.bed"
PSEUDOGENE="/home/chelseaju/Database/Pseudogene/Pseudogene_74.bed"
BED="/home/chelseaju/Database/Pseudogene/ParentENSG_Pseudogene_74.bed"
ENST_ENSG_ENSP="/home/chelseaju/Database/Ensembl/ENST_ENSG_ENSP_74.txt"
ENST2ENSG="/home/chelseaju/Database/Ensembl/script/ENST2ENSG.py"

echo "Read Correction Version 1:"

echo ""
echo "Step 0 - Loading Modules"
#module load python/2.7
#module load bedtools


echo ""
echo "Step 1 - Creating Folder"
mkdir -p $DIR/$SUBDIR/$TOPOUT/correction


echo ""
echo "Step 2 - Preparing Files"
## remove gene prefix
perl -pi -e 's/'$SUBDIR'_//g' $LASSO_Y
perl -pi -e 's/'$SUBDIR'_//g' $LASSO_X

## format X matrix, remove the "^M sign"
tr '\015' '\012' < $LASSO_X > $DIR/$SUBDIR/$TOPOUT/correction/distribution_matrix.txt
cp $DIR/$SUBDIR/$TOPOUT/correction/distribution_matrix.txt $LASSO_X

## merge gene ID with location
join -t $'\t' -1 1 -2 4 -a 1 -o 0,2.1,2.2,2.3,2.6 <(head -1 $LASSO_X | tr '\t' '\n' | sort) <(cat $PSEUDOGENE $PARENT | sort -k 4) > $DIR/$SUBDIR/$TOPOUT/correction/locus.txt
join -t $'\t' -1 1 -2 1 -a 1 <(sort -k 1 $DIR/$SUBDIR/$TOPOUT/correction/locus.txt) <(sort -k 1 $LASSO_Y) > $DIR/$SUBDIR/$TOPOUT/correction/expected_count.txt

echo ""
echo "Step 3 - Marking Multireads"

## mark multireads
python multiread_marker_v2.py -d $DIR/$SUBDIR/$TOPOUT/ -p T

## count unique reads for each region first
python uniqueread_assigner.py -d $DIR/$SUBDIR/$TOPOUT/ -p T

## assign multireads to diferent regions
python multiread_assigner_v2.py -d $DIR/$SUBDIR/$TOPOUT/ -p T

## re-assign unique reads using dynamic programming
python uniqueread_reassigner.py -d $DIR/$SUBDIR/$TOPOUT/ -p T

## update bam file
python alignment_updator.py -d $DIR/$SUBDIR/$TOPOUT/

## sort new bam file by positions and name
python bamfile_sorter.py -i corrected_hits.bam -d $DIR/$SUBDIR/$TOPOUT

# step 4: run cufflinks to compute FPKM
echo ""
echo "Step 4 - Assembling Transcripts"
cufflinks -o $DIR/$SUBDIR/$CUFFLINKS $DIR/$SUBDIR/$TOPOUT/corrected_hits_sorted.bam


