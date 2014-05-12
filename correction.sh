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
BED="/home/chelseaju/Database/Pseudogene/ParentENSG_Pseudogene_74.bed"
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

## merge gene ID with location
join -t $'\t' -1 1 -2 4 <(sort -k1 $LASSO_Y) <(sort -k4 $BED) | cut -f1-5 > $LASSO_LOCUS
mv $LASSO_LOCUS $DIR/$PREFIX/$TOPOUT/correction/expected_count.txt

echo ""
echo "Step 3 - Marking Multireads"

## mark multireads
python multiread_marker.py -d $DIR/$PREFIX/$TOPOUT/ -p T

## count unique reads for each region first
python uniqueread_assigner.py -d $DIR/$PREFIX/$TOPOUT/ -p T