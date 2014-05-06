#!/bin/bash

## this script follows analysis_v4.sh:
## 1. prepare for lasso analysis by combining the regions
## 2. run lasso analysis, and validate results
## 3. compare the lasso coefficient with distribution equations

## bash analysis_v5-2.sh DIR 


DIR=$1
TOPOUT="tophat_out"

X3A4="3X_101L_4A"
X3A6="3X_101L_6A"
X3A8="3X_101L_8A"
X3AR1="3X_101L_R1A"
X3AR2="3X_101L_R2A"
X3AR3="3X_101L_R3A"

X5A4="5X_101L_4A"
X5A6="5X_101L_6A"
X5A8="5X_101L_8A"
X5AR1="5X_101L_R1A"
X5AR2="5X_101L_R2A"
X5AR3="5X_101L_R3A"

X7A4="7X_101L_4A"
X7A6="7X_101L_6A"
X7A8="7X_101L_8A"
X7AR1="7X_101L_R1A"
X7AR2="7X_101L_R2A"
X7AR3="7X_101L_R3A"

X10A4="10X_101L_4A"
X10A6="10X_101L_6A"
X10A8="10X_101L_8A"
X10AR1="10X_101L_R1A"
X10AR2="10X_101L_R2A"
X10AR3="10X_101L_R3A"

X13A4="13X_101L_4A"
X13A6="13X_101L_6A"
X13A8="13X_101L_8A"
X13AR1="13X_101L_R1A"
X13AR2="13X_101L_R2A"
X13AR3="13X_101L_R3A"

X17A4="17X_101L_4A"
X17A6="17X_101L_6A"
X17A8="17X_101L_8A"
X17AR1="17X_101L_R1A"
X17AR2="17X_101L_R2A"
X17AR3="17X_101L_R3A"

X20A4="20X_101L_4A"
X20A6="20X_101L_6A"
X20A8="20X_101L_8A"
X20AR1="20X_101L_R1A"
X20AR2="20X_101L_R2A"
X20AR3="20X_101L_R3A"

X23A4="23X_101L_4A"
X23A6="23X_101L_6A"
X23A8="23X_101L_8A"
X23AR1="23X_101L_R1A"
X23AR2="23X_101L_R2A"
X23AR3="23X_101L_R3A"

X27A4="27X_101L_4A"
X27A6="27X_101L_6A"
X27A8="27X_101L_8A"
X27AR1="27X_101L_R1A"
X27AR2="27X_101L_R2A"
X27AR3="27X_101L_R3A"

X30A4="30X_101L_4A"
X30A6="30X_101L_6A"
X30A8="30X_101L_8A"
X30AR1="30X_101L_R1A"
X30AR2="30X_101L_R2A"
X30AR3="30X_101L_R3A"



echo ""
echo "Step 1: Unify the Unknowns"
python unknown_unifier.py -o $DIR -t genes -d $DIR/$X3A4/$TOPOUT/ $DIR/$X3A6/$TOPOUT/ $DIR/$X3A8/$TOPOUT/ $DIR/$X3AR1/$TOPOUT/ $DIR/$X3AR2/$TOPOUT/ $DIR/$X3AR3/$TOPOUT/ $DIR/$X5A4/$TOPOUT/ $DIR/$X5A6/$TOPOUT/ $DIR/$X5A8/$TOPOUT/ $DIR/$X5AR1/$TOPOUT/ $DIR/$X5AR2/$TOPOUT/ $DIR/$X5AR3/$TOPOUT/ $DIR/$X7A4/$TOPOUT/ $DIR/$X7A6/$TOPOUT/ $DIR/$X7A8/$TOPOUT/ $DIR/$X7AR1/$TOPOUT/ $DIR/$X7AR2/$TOPOUT/ $DIR/$X7AR3/$TOPOUT/ $DIR/$X10A4/$TOPOUT/ $DIR/$X10A6/$TOPOUT/ $DIR/$X10A8/$TOPOUT/ $DIR/$X10AR1/$TOPOUT/ $DIR/$X10AR2/$TOPOUT/ $DIR/$X10AR3/$TOPOUT/ $DIR/$X13A4/$TOPOUT/ $DIR/$X13A6/$TOPOUT/ $DIR/$X13A8/$TOPOUT/ $DIR/$X13AR1/$TOPOUT/ $DIR/$X13AR2/$TOPOUT/ $DIR/$X13AR3/$TOPOUT/ $DIR/$X17A4/$TOPOUT/ $DIR/$X17A6/$TOPOUT/ $DIR/$X17A8/$TOPOUT/ $DIR/$X17AR1/$TOPOUT/ $DIR/$X17AR2/$TOPOUT/ $DIR/$X17AR3/$TOPOUT/ $DIR/$X20A4/$TOPOUT/ $DIR/$X20A6/$TOPOUT/ $DIR/$X20A8/$TOPOUT/ $DIR/$X20AR1/$TOPOUT/ $DIR/$X20AR2/$TOPOUT/ $DIR/$X20AR3/$TOPOUT/ $DIR/$X23A4/$TOPOUT/ $DIR/$X23A6/$TOPOUT/ $DIR/$X23A8/$TOPOUT/ $DIR/$X23AR1/$TOPOUT/ $DIR/$X23AR2/$TOPOUT/ $DIR/$X23AR3/$TOPOUT/ $DIR/$X30A4/$TOPOUT/ $DIR/$X30A6/$TOPOUT/ $DIR/$X30A8/$TOPOUT/ $DIR/$X30AR1/$TOPOUT/ $DIR/$X30AR2/$TOPOUT/ $DIR/$X30AR3/$TOPOUT/

echo ""
echo "Step 2: Data examination, check for variance of the distribution matrix"
R --no-save --slave < data_examiniation_v2.R --args $DIR genes all all 

echo ""
echo "Step 3: Prepare Data for Lasso"
R --no-save --slave < lasso_preparation_v3.R --args $DIR genes 54

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
