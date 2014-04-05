#!/bin/bash
## sh evaluation.sh DIR A B BAM SUFFIX

DIR=$1
COVERAGE=$2
READLEN=$3
ABUNDANCE=$4


LASSO_COEFF="lasso_coefficient.txt"
LASSO_COEFF_LOCUS="lasso_coefficient_locus.txt"

## MAC
BED="/Users/Chelsea/Bioinformatics/CJDatabase/Pseudogene/ParentENSG_Pseudogene_74.bed"
#ENST_ENSG_ENSP="/home/chelseaju/Database/Ensembl/ENST_ENSG_ENSP_74.txt"
#ENST2ENSG="/home/chelseaju/Database/Ensembl/script/ENST2ENSG.py"

## HOFFMAN
#BED="/u/home/c/chelseaj/project/database/Pseudogene/ParentENST_Pseudogene_74.bed"
#ENST_ENSG_ENSP="/u/home/c/chelseaj/project/database/Ensembl/ENST_ENSG_ENSP_74.txt"
#ENST2ENSG="/u/home/c/chelseaj/project/database/Ensembl/script/ENST2ENSG.py"

## LAB
#ENST_ENSG_ENSP="/home/chelseaju/Database/Ensembl/ENST_ENSG_ENSP_74.txt"
#ENST2ENSG="/home/chelseaju/Database/Ensembl/script/ENST2ENSG.py"

## merge gene ID with location
join -1 1 -2 4 <(sort $DIR/$LASSO_COEFF) <(sort -k4 $BED) > $DIR/LASSO_COEFF_LOCUS