#!/bin/bash
## sh evaluation.sh DIR A B BAM SUFFIX

DIR=$1
SAMPLEA=$2
SAMPLEB=$3
BAM=$4
SUFFIX=$5

TOPHAT="tophat_out"
CUFFLINK=cufflinks_$SUFFIX
CUFFMERGE=$DIR/$SAMPLEA"_"$SAMPLEB/cuffdiff_$SUFFIX
ASSEMBLY=$CUFFMERGE/assemblies.txt


#GTF="/u/home/c/chelseaj/project/database/Ensembl/Homo_sapiens.GRCh37.74.gtf"
#GENOME="/u/home/c/chelseaj/project/database/Ensembl/Genome/Homo_sapiens_GRCh37_74_genome.fa"

#GTF="/home/chelseaju/Database/Ensembl/Homo_sapiens.GRCh37.74.gtf"
GTF="$DIR/SELECT_GENES.gtf"
GENOME="/home/chelseaju/Database/Ensembl/Genome/Homo_sapiens_GRCh37_74_genome.fa"

## Step 0: Create directory
echo ""
echo "Step 0 - Creating Directories"
echo "mkdir -p $CUFFMERGE"
mkdir -p $CUFFMERGE

## Step 1: Create Assemblies.txt
echo ""
echo "Step 1 - Creating Assemblies.txt"
echo "$DIR/$SAMPLEA/$CUFFLINK/transcripts.gtf" > $ASSEMBLY
#echo "$DIR/$SAMPLEA-2/$CUFFLINK/transcripts.gtf" >> $ASSEMBLY
#echo "$DIR/$SAMPLEA-3/$CUFFLINK/transcripts.gtf" >> $ASSEMBLY
echo "$DIR/$SAMPLEB/$CUFFLINK/transcripts.gtf" >> $ASSEMBLY
#echo "$DIR/$SAMPLEB-2/$CUFFLINK/transcripts.gtf" >> $ASSEMBLY
#echo "$DIR/$SAMPLEB-3/$CUFFLINK/transcripts.gtf" >> $ASSEMBLY
echo "Write to $ASSEMBLY"

## Step 2: Running cuffmerge
echo ""
echo "Step 2 - Running Cuffmerge"
echo "cuffmerge -o $CUFFMERGE -g $GTF -s $GENOME $ASSEMBLY"
#cuffmerge -p 8 -o $CUFFMERGE -g $GTF -s $GENOME $ASSEMBLY
cuffmerge -p 8 -o $CUFFMERGE  -g $GTF -s $GENOME $ASSEMBLY 

## Step 3: Runing cuffdiff
echo ""
echo "Step 3 - Running Cuffdiff"
echo "cuffdiff -o $CUFFMERGE -b $GENOME -L $SAMPLEA,$SAMPLEB -u $CUFFMERGE/merged.gtf $DIR/$SAMPLEA/$TOPHAT/$BAM $DIR/$SAMPLEB/$TOPHAT/$BAM"

cuffdiff -p 8 -o $CUFFMERGE -b $GENOME -L $SAMPLEA,$SAMPLEB -u $CUFFMERGE/merged.gtf $DIR/$SAMPLEA/$TOPHAT/$BAM $DIR/$SAMPLEB/$TOPHAT/$BAM



