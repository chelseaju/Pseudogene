#!/bin/bash

## this script follows analysis_v4.sh:
## 1. prepare for lasso analysis by combining the regions
## 2. run lasso analysis, and validate results
## 3. compare the lasso coefficient with distribution equations

## bash analysis_v4-2.sh DIR 


DIR=$1
TOPOUT="tophat_out"

X10A1="10X_101L_1A"
X10AR1="10X_101L_R1A"
X10AR2="10X_101L_R2A"
X10AR3="10X_101L_R3A"

X20A1="20X_101L_1A"
X20AR1="20X_101L_R1A"
X20AR2="20X_101L_R2A"
X20AR3="20X_101L_R3A"

X30A1="30X_101L_1A"
X30AR1="30X_101L_R1A"
X30AR2="30X_101L_R2A"
X30AR3="30X_101L_R3A"

X3AR1="3X_101L_R1A"
X5AR1="5X_101L_R1A"
X7AR1="7X_101L_R1A"
X13AR1="13X_101L_R1A"
X17AR1="17X_101L_R1A"
X23AR1="23X_101L_R1A"
X27AR1="27X_101L_R1A"


echo ""
echo "Step 1: Unify the Unknowns"
python unknown_unifier.py -o $DIR -t ENSG -d $DIR/$X10A1/$TOPOUT/ $DIR/$X10AR1/$TOPOUT/ $DIR/$X10AR2/$TOPOUT/ $DIR/$X10AR3/$TOPOUT/ $DIR/$X20A1/$TOPOUT/ $DIR/$X20AR1/$TOPOUT/ $DIR/$X20AR2/$TOPOUT/ $DIR/$X20AR3/$TOPOUT/ $DIR/$X30A1/$TOPOUT/ $DIR/$X30AR1/$TOPOUT/ $DIR/$X30AR2/$TOPOUT/ $DIR/$X30AR3/$TOPOUT/ $DIR/$X3AR1/$TOPOUT/ $DIR/$X5AR1/$TOPOUT/ $DIR/$X7AR1/$TOPOUT/ $DIR/$X13AR1/$TOPOUT/ $DIR/$X17AR1/$TOPOUT/ $DIR/$X23AR1/$TOPOUT/ $DIR/$X27AR1/$TOPOUT/

echo ""
echo "Step 2: Prepare Data for Lasso"
mkdir -p $DIR/LassoTraining
R --no-save --slave < lasso_preparation_v2.R --args $DIR ENSG

echo ""
echo "Step 3: Training Coefficient with Lasso"
#R --no-save --slave < lasso_trainer.R --args $DIR ENSG

echo ""
echo "Step 3: Prepare Data for Validation"
#mkdir -p $DIR/LassoValidation
#cp $DIR/$X30A1/$TOPOUT/ENSG_distribution.matrix $DIR/LassoValidation/ENSG_30X1A_matrix_X.txt
#cp $DIR/$X30AR2/$TOPOUT/ENSG_distribution.matrix $DIR/LassoValidation/ENSG_30XR2A_matrix_X.txt
#cp $DIR/$X30AR3/$TOPOUT/ENSG_distribution.matrix $DIR/LassoValidation/ENSG_30XR3A_matrix_X.txt
#cp $DIR/$X30A1/$TOPOUT/ENSG_expected_read_count.txt $DIR/LassoValidation/ENSG_30X1A_expected_Y.txt
#cp $DIR/$X30AR2/$TOPOUT/ENSG_expected_read_count.txt $DIR/LassoValidation/ENSG_30XR2A_expected_Y.txt
#cp $DIR/$X30AR3/$TOPOUT/ENSG_expected_read_count.txt $DIR/LassoValidation/ENSG_30XR3A_expected_Y.txt

echo ""
echo "Step 4: Lasso Prediction"
#R --no-save --slave < lasso_predictor_v4.R --args $DIR 30X1A ENSG
#R --no-save --slave < lasso_predictor_v4.R --args $DIR 30XR2A ENSG
#R --no-save --slave < lasso_predictor_v4.R --args $DIR 30XR3A ENSG

echo ""
echo "Step 5: Lasso Validation"
#R --no-save --slave < lasso_validator_v4.R --args $DIR 30X1A ENSG
#R --no-save --slave < lasso_validator_v4.R --args $DIR 30XR2A ENSG
#R --no-save --slave < lasso_validator_v4.R --args $DIR 30XR3A ENSG
