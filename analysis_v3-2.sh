#!/bin/bash

## this script follows analysis_v3.sh:
## 1. prepare for lasso analysis by combining the regions
## 2. run lasso analysis, and validate results
## 3. compare the lasso coefficient with distribution equations

## bash analysis_v3-2.sh DIR 

DIR=$1
TOPOUT="tophat_out"

X10A10="10X_100L_10A"
X10A30="10X_100L_30A"
X10A50="10X_100L_50A"
X20A10="20X_100L_10A"
X20A30="20X_100L_30A"
X20A50="20X_100L_50A"
X30A10="30X_100L_10A"
X30A30="30X_100L_30A"
X30A50="30X_100L_50A"


LOOKUP=$DIR/"lookup_table.txt"

echo ""
echo "Step 1: Prepare for Lasso Analysis"
python unknown_unifier.py -o $DIR -t genes -d $DIR/$X10A10/$TOPOUT/ $DIR/$X10A30/$TOPOUT/ $DIR/$X10A50/$TOPOUT/ $DIR/$X20A10/$TOPOUT/ $DIR/$X20A30/$TOPOUT/ $DIR/$X20A50/$TOPOUT/ $DIR/$X30A10/$TOPOUT/ $DIR/$X30A30/$TOPOUT/ $DIR/$X30A50/$TOPOUT/

python unknown_updator.py -d $DIR -f $DIR/$X10A10/$TOPOUT/genes_distribution.matrix -t genes
python unknown_updator.py -d $DIR -f $DIR/$X10A10/$TOPOUT/genes_distribution.eqn -t genes
python unknown_updator.py -d $DIR -f $DIR/$X10A30/$TOPOUT/genes_distribution.matrix -t genes
python unknown_updator.py -d $DIR -f $DIR/$X10A30/$TOPOUT/genes_distribution.eqn -t genes
python unknown_updator.py -d $DIR -f $DIR/$X10A50/$TOPOUT/genes_distribution.matrix -t genes
python unknown_updator.py -d $DIR -f $DIR/$X10A50/$TOPOUT/genes_distribution.eqn -t genes

python unknown_updator.py -d $DIR -f $DIR/$X20A10/$TOPOUT/genes_distribution.matrix -t genes
python unknown_updator.py -d $DIR -f $DIR/$X20A10/$TOPOUT/genes_distribution.eqn -t genes
python unknown_updator.py -d $DIR -f $DIR/$X20A30/$TOPOUT/genes_distribution.matrix -t genes
python unknown_updator.py -d $DIR -f $DIR/$X20A30/$TOPOUT/genes_distribution.eqn -t genes
python unknown_updator.py -d $DIR -f $DIR/$X20A50/$TOPOUT/genes_distribution.matrix -t genes
python unknown_updator.py -d $DIR -f $DIR/$X20A50/$TOPOUT/genes_distribution.eqn -t genes

python unknown_updator.py -d $DIR -f $DIR/$X30A10/$TOPOUT/genes_distribution.matrix -t genes
python unknown_updator.py -d $DIR -f $DIR/$X30A10/$TOPOUT/genes_distribution.eqn -t genes
python unknown_updator.py -d $DIR -f $DIR/$X30A30/$TOPOUT/genes_distribution.matrix -t genes
python unknown_updator.py -d $DIR -f $DIR/$X30A30/$TOPOUT/genes_distribution.eqn -t genes
python unknown_updator.py -d $DIR -f $DIR/$X30A50/$TOPOUT/genes_distribution.matrix -t genes
python unknown_updator.py -d $DIR -f $DIR/$X30A50/$TOPOUT/genes_distribution.eqn -t genes


python unknown_unifier.py -o $DIR -t transcripts -d $DIR/$X10A10/$TOPOUT/ $DIR/$X10A30/$TOPOUT/ $DIR/$X10A50/$TOPOUT/ $DIR/$X20A10/$TOPOUT/ $DIR/$X20A30/$TOPOUT/ $DIR/$X20A50/$TOPOUT/ $DIR/$X30A10/$TOPOUT/ $DIR/$X30A30/$TOPOUT/ $DIR/$X30A50/$TOPOUT/

python unknown_updator.py -d $DIR -f $DIR/$X10A10/$TOPOUT/transcripts_distribution.matrix -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X10A10/$TOPOUT/transcripts_distribution.eqn -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X10A30/$TOPOUT/transcripts_distribution.matrix -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X10A30/$TOPOUT/transcripts_distribution.eqn -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X10A50/$TOPOUT/transcripts_distribution.matrix -t transcripts 
python unknown_updator.py -d $DIR -f $DIR/$X10A50/$TOPOUT/transcripts_distribution.eqn -t transcripts

python unknown_updator.py -d $DIR -f $DIR/$X20A10/$TOPOUT/transcripts_distribution.matrix -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X20A10/$TOPOUT/transcripts_distribution.eqn -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X20A30/$TOPOUT/transcripts_distribution.matrix -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X20A30/$TOPOUT/transcripts_distribution.eqn -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X20A50/$TOPOUT/transcripts_distribution.matrix -t transcripts 
python unknown_updator.py -d $DIR -f $DIR/$X20A50/$TOPOUT/transcripts_distribution.eqn -t transcripts

python unknown_updator.py -d $DIR -f $DIR/$X30A10/$TOPOUT/transcripts_distribution.matrix -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X30A10/$TOPOUT/transcripts_distribution.eqn -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X30A30/$TOPOUT/transcripts_distribution.matrix -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X30A30/$TOPOUT/transcripts_distribution.eqn -t transcripts
python unknown_updator.py -d $DIR -f $DIR/$X30A50/$TOPOUT/transcripts_distribution.matrix -t transcripts 
python unknown_updator.py -d $DIR -f $DIR/$X30A50/$TOPOUT/transcripts_distribution.eqn -t transcripts


echo ""
echo "Step 2: Training Coefficient with Lasso"
#R --no-save --slave < lasso_trainer.R --args $DIR genes


echo ""
echo "Step 3: Lasso Prediction"
#R --no-save --slave < lasso_predictor.R --args $DIR $X30A10/$TOPOUT genes glmnet
#R --no-save --slave < lasso_predictor.R --args $DIR $X30A30/$TOPOUT genes glmnet
#R --no-save --slave < lasso_predictor.R --args $DIR $X30A50/$TOPOUT genes glmnet
#R --no-save --slave < lasso_predictor.R --args $DIR $X30A10/$TOPOUT genes genlasso
#R --no-save --slave < lasso_predictor.R --args $DIR $X30A30/$TOPOUT genes genlasso
#R --no-save --slave < lasso_predictor.R --args $DIR $X30A50/$TOPOUT genes genlasso


echo ""
echo "Step 4: Lasso Validation"
#R --no-save --slave < lasso_validator.R --args $DIR $X30A10/$TOPOUT genes glmnet
#R --no-save --slave < lasso_validator.R --args $DIR $X30A30/$TOPOUT genes glmnet
#R --no-save --slave < lasso_validator.R --args $DIR $X30A50/$TOPOUT genes glmnet
#R --no-save --slave < lasso_validator.R --args $DIR $X30A10/$TOPOUT genes genlasso
#R --no-save --slave < lasso_validator.R --args $DIR $X30A30/$TOPOUT genes genlasso
#R --no-save --slave < lasso_validator.R --args $DIR $X30A50/$TOPOUT genes genlasso


















