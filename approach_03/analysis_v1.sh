## sh analysis.sh reference human_72.pseudogene human72 100 100
DIR=$1
PSEUDO=$2
PREFIX=$3
GENOME="/u/scratch/c/chelseaj/database/genome.fa"
GENOMEINDEX="/u/scratch/c/chelseaj/database/Bowtie2Index/genome"
COVERAGE=$4
READLEN=$5
INSERT=200
INSERTSTD=20
CUFFLINKS="cufflinks_out"
TOPOUT="tophat_out"
SUBDIR=$COVERAGE"X_"$READLEN"L"

echo "Data Analysis Version 1:"

# step 1: parse out the expression
CMD4="python expression_parser_from_bam.py -p $DIR/$PREFIX""_parents.bed -s $DIR/$PREFIX""_pseudo.bed -r $PSEUDO -e $DIR/$SUBDIR/$TOPOUT/accepted_hits.bam -o $DIR/$SUBDIR/$PREFIX"

echo "\n\nStep 1: Parsing Expression"
echo ">$CMD4"
$CMD4


# step 2: run linear regression analysis
CMD5="R --no-save --slave < model_expression_value.R --args $DIR/$SUBDIR/$PREFIX"".expression $DIR/$SUBDIR/$PREFIX"
echo "\n\nStep 2: Model Pseudogene and Parent Genes Expression"
echo ">$CMD5"

R --no-save --slave < model_expression_value.R --args $DIR/$SUBDIR/$PREFIX"".expression $DIR/$SUBDIR/$PREFIX



