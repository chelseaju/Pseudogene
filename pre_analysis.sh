#!/bin/bash

## sh pre_analysis.sh reference human_72.pseudogene human72

GENOME="/u/scratch/c/chelseaj/database/genome.fa"
MUTATION=0.001
DIR=$1
PSEUDO=$2
PREFIX=$3

## Step 0: Create necessary directory
CMD0="mkdir -p $DIR"
$CMD0

## Step 1: Generate personal genome
#CMD1 = "python personalized_genome.py -g $REFERENCE -m $MUTATION -o $DIR/personal_genome"
#echo "\n\nStep 1: Generate Personal Genome\n"
#echo ">$CMD1"
#$CMD1

## Step 2: Extract pseudogene information 
CMD2="python pseudogene_pseudo_parser.py -i $PSEUDO -o $DIR/$PREFIX"
echo "\n\nStep 2: Extract Pseudogene Information:"
echo ">$CMD2"
$CMD2

## Step 3: retrieve parent information: consensus mRNA sequence
CMD3="R --no-save --slave < pseudogene_parent_retriever.R --args $DIR/$PREFIX""_parents.peptide $DIR/$PREFIX"
echo "\n\nStep 3: Retrieve Parent Information"
echo ">$CMD3"
R --no-save --slave < pseudogene_parent_retriever.R --args $DIR/$PREFIX""_parents.peptide $DIR/$PREFIX

# Step 4: construct parent sequence with personalized or given genome
CMD4="python pseudogene_parent_bed_retriever.py -i $DIR/$PREFIX""_parents.fa -o $DIR/$PREFIX"
CMD4a="bedtools getfasta -name -s -split -fi $GENOME -bed $DIR/$PREFIX""_parents.bed -fo $DIR/$PREFIX""_parents.fa"

echo "\n\nStep 4: Construct Parent Sequence with Personal (or given) Genome"
echo ">$CMD4"
$CMD4

echo ">$CMD4a"
$CMD4a


### Execute this step if smaller set of data is used
CMD5="head 500 $DIR/$PREFIX""_parents.bed > $DIR/$PREFIX""_top500_parents.bed"
CMD5a="bedtools getfasta -name -s -split -fi $GENOME -bed $DIR/$PREFIX""_top500_parents.bed -fo $DIR/$PREFIX""_top500_parents.fa"

echo "\n\nStep 5: Select a subset of parent data"
echo ">$CMD5"
#$CMD5
echo ">$CMD5a"
#$CMD5a

### Execute this step if the filter set is used (remove sequence from non-conventional chromosome)
CMD6="grep 'chr[0-9XY]' $DIR/$PREFIX""_parents.bed > $DIR/$PREFIX""_filter_parents.bed"
CMD6a="bedtools getfasta -name -s -split -fi $GENOME -bed $DIR/$PREFIX""_filter_parents.bed -fo $DIR/$PREFIX""_filter_parents.fa"

echo "\n\nStep 6: Filter out non-conventional chromosome"
echo ">$CMD6"
#$CMD6
echo ">$CMD6a"
#$CMD6a

### Execute this step to randomly assign expession level to different genes
CMD7="python assign_random_expression.py -i $DIR/$PREFIX""_parents.bed -o $DIR/$PREFIX""_replicates_parents -n 10 -x 100"
CMD7a="bedtools getfasta -name -s -split -fi $GENOME -bed $DIR/$PREFIX""_replicates_parents.bed -fo $DIR/$PREFIX""_replicates_parents.fa"

echo "\n\nStep 7: Randomly Assign Expression Value"
echo ">$CMD7"
$CMD7
echo ">$CMD7a"
$CMD7a


