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

echo "Step 1: Unify the Unknowns"
python unknown_unifier.py -o $DIR -t genes -d $DIR/$X10A1/$TOPOUT/ $DIR/$X10AR1/$TOPOUT/ $DIR/$X10AR2/$TOPOUT/ $DIR/$X10AR3/$TOPOUT/ $DIR/$X20A1/$TOPOUT/ $DIR/$X20AR1/$TOPOUT/ $DIR/$X20AR2/$TOPOUT/ $DIR/$X20AR3/$TOPOUT/ $DIR/$X30A1/$TOPOUT/ $DIR/$X30AR1/$TOPOUT/ $DIR/$X30AR2/$TOPOUT/ $DIR/$X30AR3/$TOPOUT/ 

echo ""
echo "Step 2: Prepare Data for Lasso"
mkdir -p $DIR/LassoTraining
R --no-save --slave < lasso_preparation.R --args $DIR genes

echo ""
echo "Step 3: Training Coefficient with Lasso"
#R --no-save --slave < lasso_trainer.R --args $DIR ENSG

echo ""
echo "Step 3: Making Lasso Folder"
mkdir -p $DIR/$X30A1/lasso_out
mkdir -p $DIR/$X30AR1/lasso_out
mkdir -p $DIR/$X30AR2/lasso_out
mkdir -p $DIR/$X30AR3/lasso_out


echo ""
echo "Step 4: Lasso Prediction"
#R --no-save --slave < lasso_predictor.R --args $DIR $X30A1 ENSG glmnet
#R --no-save --slave < lasso_predictor.R --args $DIR $X30AR1 ENSG glmnet
#R --no-save --slave < lasso_predictor.R --args $DIR $X30AR2 ENSG glmnet
#R --no-save --slave < lasso_predictor.R --args $DIR $X30AR3 ENSG glmnet
#R --no-save --slave < lasso_predictor.R --args $DIR $X30A1 ENSG genlasso
#R --no-save --slave < lasso_predictor.R --args $DIR $X30AR1 ENSG genlasso
#R --no-save --slave < lasso_predictor.R --args $DIR $X30AR2 ENSG genlasso
#R --no-save --slave < lasso_predictor.R --args $DIR $X30AR3 ENSG genlasso


echo ""
echo "Step 5: Lasso Validation"
#R --no-save --slave < lasso_validator.R --args $DIR $X30A1 ENSG glmnet
#R --no-save --slave < lasso_validator.R --args $DIR $X30AR1 ENSG glmnet
#R --no-save --slave < lasso_validator.R --args $DIR $X30AR2 ENSG glmnet
#R --no-save --slave < lasso_validator.R --args $DIR $X30AR3 ENSG glmnet
#R --no-save --slave < lasso_validator.R --args $DIR $X30A1 ENSG genlasso
#R --no-save --slave < lasso_validator.R --args $DIR $X30AR1 ENSG genlasso
#R --no-save --slave < lasso_validator.R --args $DIR $X30AR2 ENSG genlasso
#R --no-save --slave < lasso_validator.R --args $DIR $X30AR3 ENSG genlasso
