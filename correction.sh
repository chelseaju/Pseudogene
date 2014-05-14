#!/bin/bash
## sh evaluation.sh DIR A B BAM SUFFIX

DIR=$1
PREFIX=$2"X_"$3"L_"$4
LASSODIR=$5
TOPOUT="tophat_out"
LASSO_X="$DIR/$LASSODIR/$PREFIX""_reconstructed_x_matrix.txt"
LASSO_Y="$DIR/$LASSODIR/$PREFIX""_validationY_v2.txt"
LASSO_LOCUS="$DIR/$LASSODIR/$PREFIX""_validationY_locus.txt"


## MAC
#BED="/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/ParentENSG_Pseudogene_74.bed"
#ENST_ENSG_ENSP="/home/chelseaju/Database/Ensembl/ENST_ENSG_ENSP_74.txt"
#ENST2ENSG="/home/chelseaju/Database/Ensembl/script/ENST2ENSG.py"

## HOFFMAN
#BED="/u/home/c/chelseaj/project/database/Pseudogene/ParentENST_Pseudogene_74.bed"
#ENST_ENSG_ENSP="/u/home/c/chelseaj/project/database/Ensembl/ENST_ENSG_ENSP_74.txt"
#ENST2ENSG="/u/home/c/chelseaj/project/database/Ensembl/script/ENST2ENSG.py"

## LAB
PARENT="/home/chelseaju/Database/Ensembl/ENSG_74.bed"
PSEUDOGENE="/home/chelseaju/Database/Pseudogene/Pseudogene_74.bed"
#ENST_ENSG_ENSP="/home/chelseaju/Database/Ensembl/ENST_ENSG_ENSP_74.txt"
#ENST2ENSG="/home/chelseaju/Database/Ensembl/script/ENST2ENSG.py"

echo "Read Correction Version 1:"

echo ""
echo "Step 1 - Creating Folder"
mkdir -p $DIR/$PREFIX/$TOPOUT/correction


echo ""
echo "Step 2 - Preparing Files"
## remove gene prefix
perl -pi -e 's/'$PREFIX'_//g' $LASSO_Y
perl -pi -e 's/'$PREFIX'_//g' $LASSO_X

## format X matrix, remove the "^M sign"
tr '\015' '\012' < $LASSO_X > $DIR/$PREFIX/$TOPOUT/correction/distribution_matrix.txt
cp $DIR/$PREFIX/$TOPOUT/correction/distribution_matrix.txt $LASSO_X

## merge gene ID with location
join -t $'\t' -1 1 -2 4 -a 1 -o 0,2.1,2.2,2.3 <(head -1 $LASSO_X | tr '\t' '\n' | sort) <(cat $PSEUDOGENE $PARENT | sort -k 4) > $DIR/$PREFIX/$TOPOUT/correction/locus.txt
join -t $'\t' -1 1 -2 1 -a 1 <(sort -k 1 $DIR/$PREFIX/$TOPOUT/correction/locus.txt) <(sort -k 1 $LASSO_Y) > $DIR/$PREFIX/$TOPOUT/correction/expected_count.txt


echo ""
echo "Step 3 - Marking Multireads"

## mark multireads
python multiread_marker.py -d $DIR/$PREFIX/$TOPOUT/ -p T

## count unique reads for each region first
python uniqueread_assigner.py -d $DIR/$PREFIX/$TOPOUT/ -p T

## assign multireads to diferent regions
python multiread_assigner.py -d $DIR/$PREFIX/$TOPOUT/ -p T

## re-assign unique reads using dynamic programming
python uniqueread_reassigner.py -d $DIR/$PREFIX/$TOPOUT/ -p T