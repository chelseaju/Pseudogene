#!/bin/bash

## this script follows analysis_v4.sh:
## 1. prepare for lasso analysis by combining the regions
## 2. run lasso analysis, and validate results
## 3. compare the lasso coefficient with distribution equations

## bash analysis_v5-2.sh DIR 


DIR=$1
TOPOUT="tophat_out"

X3A1="3X_101L_1A"
X3A5="3X_101L_5A"
X3A10="3X_101L_10A"
X3AR1="3X_101L_R1A"
X3AR2="3X_101L_R2A"
X3AR3="3X_101L_R3A"

X5A1="5X_101L_1A"
X5A5="5X_101L_5A"
X5A10="5X_101L_10A"
X5AR1="5X_101L_R1A"
X5AR2="5X_101L_R2A"
X5AR3="5X_101L_R3A"

X7A1="7X_101L_1A"
X7A5="7X_101L_5A"
X7A10="7X_101L_10A"
X7AR1="7X_101L_R1A"
X7AR2="7X_101L_R2A"
X7AR3="7X_101L_R3A"

X10A1="10X_101L_1A"
X10A5="10X_101L_5A"
X10A10="10X_101L_10A"
X10AR1="10X_101L_R1A"
X10AR2="10X_101L_R2A"
X10AR3="10X_101L_R3A"

X13A1="13X_101L_1A"
X13A5="13X_101L_5A"
X13A10="13X_101L_10A"
X13AR1="13X_101L_R1A"
X13AR2="13X_101L_R2A"
X13AR3="13X_101L_R3A"

X17A1="17X_101L_1A"
X17A5="17X_101L_5A"
X17A10="17X_101L_10A"
X17AR1="17X_101L_R1A"
X17AR2="17X_101L_R2A"
X17AR3="17X_101L_R3A"

X20A1="20X_101L_1A"
X20A5="20X_101L_5A"
X20A10="20X_101L_10A"
X20AR1="20X_101L_R1A"
X20AR2="20X_101L_R2A"
X20AR3="20X_101L_R3A"

X23A1="23X_101L_1A"
X23A5="23X_101L_5A"
X23A10="23X_101L_10A"
X23AR1="23X_101L_R1A"
X23AR2="23X_101L_R2A"
X23AR3="23X_101L_R3A"

X27A1="27X_101L_1A"
X27A5="27X_101L_5A"
X27A10="27X_101L_10A"
X27AR1="27X_101L_R1A"
X27AR2="27X_101L_R2A"
X27AR3="27X_101L_R3A"

X30A1="30X_101L_1A"
X30A5="30X_101L_5A"
X30A10="30X_101L_10A"
X30AR1="30X_101L_R1A"
X30AR2="30X_101L_R2A"
X30AR3="30X_101L_R3A"



echo ""
echo "Step 1: Unify the Unknowns"
python unknown_unifier.py -o $DIR -t genes -d $DIR/$X3A1/$TOPOUT/ $DIR/$X3A5/$TOPOUT/ $DIR/$X3A10/$TOPOUT/ $DIR/$X3AR1/$TOPOUT/ $DIR/$X3AR2/$TOPOUT/ $DIR/$X3AR3/$TOPOUT/ $DIR/$X5A1/$TOPOUT/ $DIR/$X5A5/$TOPOUT/ $DIR/$X5A10/$TOPOUT/ $DIR/$X5AR1/$TOPOUT/ $DIR/$X5AR2/$TOPOUT/ $DIR/$X5AR3/$TOPOUT/ $DIR/$X7A1/$TOPOUT/ $DIR/$X7A5/$TOPOUT/ $DIR/$X7A10/$TOPOUT/ $DIR/$X7AR1/$TOPOUT/ $DIR/$X7AR2/$TOPOUT/ $DIR/$X7AR3/$TOPOUT/ $DIR/$X10A1/$TOPOUT/ $DIR/$X10A5/$TOPOUT/ $DIR/$X10A10/$TOPOUT/ $DIR/$X10AR1/$TOPOUT/ $DIR/$X10AR2/$TOPOUT/ $DIR/$X10AR3/$TOPOUT/ $DIR/$X13A1/$TOPOUT/ $DIR/$X13A5/$TOPOUT/ $DIR/$X13A10/$TOPOUT/ $DIR/$X13AR1/$TOPOUT/ $DIR/$X13AR2/$TOPOUT/ $DIR/$X13AR3/$TOPOUT/ $DIR/$X17A1/$TOPOUT/ $DIR/$X17A5/$TOPOUT/ $DIR/$X17A10/$TOPOUT/ $DIR/$X17AR1/$TOPOUT/ $DIR/$X17AR2/$TOPOUT/ $DIR/$X17AR3/$TOPOUT/ $DIR/$X20A1/$TOPOUT/ $DIR/$X20A5/$TOPOUT/ $DIR/$X20A10/$TOPOUT/ $DIR/$X20AR1/$TOPOUT/ $DIR/$X20AR2/$TOPOUT/ $DIR/$X20AR3/$TOPOUT/ $DIR/$X23A1/$TOPOUT/ $DIR/$X23A5/$TOPOUT/ $DIR/$X23A10/$TOPOUT/ $DIR/$X23AR1/$TOPOUT/ $DIR/$X23AR2/$TOPOUT/ $DIR/$X23AR3/$TOPOUT/ $DIR/$X30A1/$TOPOUT/ $DIR/$X30A5/$TOPOUT/ $DIR/$X30A10/$TOPOUT/ $DIR/$X30AR1/$TOPOUT/ $DIR/$X30AR2/$TOPOUT/ $DIR/$X30AR3/$TOPOUT/

echo ""
echo "Step 2: Prepare Data for Lasso"
#mkdir -p $DIR/LassoTraining
#R --no-save --slave < lasso_preparation_v2.R --args $DIR genes

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
