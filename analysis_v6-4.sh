#!/bin/bash

## this script follows analysis_v4.sh:
## 1. prepare for lasso analysis by combining the regions
## 2. run lasso analysis, and validate results
## 3. compare the lasso coefficient with distribution equations

## bash analysis_v5-2.sh DIR 


DIR=$1
TOPOUT="tophat_out"
EXP=('5X_101L_4A' '5X_101L_6A' '5X_101L_8A' '5X_101L_R1A' '5X_101L_R2A' '5X_101L_R3A'
	'7X_101L_4A' '7X_101L_6A' '7X_101L_8A' '7X_101L_R1A' '7X_101L_R2A' '7X_101L_R3A'
	'10X_101L_4A' '10X_101L_6A' '10X_101L_8A' '10X_101L_R1A' '10X_101L_R2A' '10X_101L_R3A'
	'13X_101L_4A' '13X_101L_6A' '13X_101L_8A' '13X_101L_R1A' '13X_101L_R2A' '13X_101L_R3A'
	'17X_101L_4A' '17X_101L_6A' '17X_101L_8A' '17X_101L_R1A' '17X_101L_R2A' '17X_101L_R3A'
	'20X_101L_4A' '20X_101L_6A' '20X_101L_8A' '20X_101L_R1A' '20X_101L_R2A' '20X_101L_R3A'
	'23X_101L_4A' '23X_101L_6A' '23X_101L_8A' '23X_101L_R1A' '23X_101L_R2A' '23X_101L_R3A'
	'27X_101L_4A' '27X_101L_6A' '27X_101L_8A' '27X_101L_R1A' '27X_101L_R2A' '27X_101L_R3A'
	'30X_101L_4A' '30X_101L_6A' '30X_101L_8A' '30X_101L_R1A' '30X_101L_R2A' '30X_101L_R3A')

#EXP=('5X_101L_4A' '5X_101L_6A' '5X_101L_8A')


EXP=( "${EXP[@]/#/$DIR/}" )
EXP=( "${EXP[@]/%//$TOPOUT}" )

echo ""
echo "Step 1: Unify the Unknowns"
echo "python unknown_unifier_v2.py -o $DIR -t genes -d ${EXP[@]}"
python unknown_unifier_v2.py -o $DIR -t genes -d ${EXP[@]}

echo ""
echo "Step 2: Update the Unknonws"

for subdir in "${EXP[@]}"
do
	echo ""
	echo "Sub-directory ${subdir} ===================="
	python unknown_updator_v2.py -d $DIR -i $subdir/genes_distribution.matrix -t genes
	cp $subdir/genes_distribution.matrix $subdir/genes_distribution.backup
	R --no-save --slave < unknown_merger_v1.R --args $subdir genes
	rm $subdir/genes_distribution.backup
done



echo ""
echo "Step 2: Data examination, check for variance of the distribution matrix"
R --no-save --slave < data_examination_v2.R --args $DIR genes all all 

echo ""
echo "Step 3: Prepare Data for Lasso"
#R --no-save --slave < lasso_preparation_v3.R --args $DIR genes 54

echo ""
echo "Step 3: Training Coefficient with Lasso"
#R --no-save --slave < lasso_trainer.R --args $DIR genes

echo ""
echo "Step 3: Prepare Data for Validation"
#mkdir -p $DIR/LassoValidation
#cp $DIR/$X30A1/$TOPOUT/genes_distribution.matrix $DIR/LassoValidation/genes_30X1A_matrix_X.txt
#cp $DIR/$X30AR2/$TOPOUT/genes_distribution.matrix $DIR/LassoValidation/genes_30XR2A_matrix_X.txt
#cp $DIR/$X30AR3/$TOPOUT/genes_distribution.matrix $DIR/LassoValidation/genes_30XR3A_matrix_X.txt
#cp $DIR/$X30A1/$TOPOUT/genes_expected_read_count.txt $DIR/LassoValidation/genes_30X1A_expected_Y.txt
#cp $DIR/$X30AR2/$TOPOUT/genes_expected_read_count.txt $DIR/LassoValidation/genes_30XR2A_expected_Y.txt
#cp $DIR/$X30AR3/$TOPOUT/genes_expected_read_count.txt $DIR/LassoValidation/genes_30XR3A_expected_Y.txt

echo ""
echo "Step 4: Lasso Prediction"
#R --no-save --slave < lasso_predictor_v4.R --args $DIR 30X1A genes
#R --no-save --slave < lasso_predictor_v4.R --args $DIR 30XR2A genes
#R --no-save --slave < lasso_predictor_v4.R --args $DIR 30XR3A genes

echo ""
echo "Step 5: Lasso Validation"
#R --no-save --slave < lasso_validator_v4.R --args $DIR 30X1A genes
#R --no-save --slave < lasso_validator_v4.R --args $DIR 30XR2A genes
#R --no-save --slave < lasso_validator_v4.R --args $DIR 30XR3A genes
