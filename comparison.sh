#!/bin/bash
## sh evaluation.sh DIR A B BAM SUFFIX

DIR=$1
A=$2
B=$3
BEFORE="cufflinks_before"
AFTER="cufflinks_after"


echo ""
echo "Step 1 - Creating Folder"
mkdir -p $DIR"/"$A"_"$B"/"

echo ""
echo "Step 2 - Link Cufflink ID"
#python cufflinks_identifier.py -d 
#python cufflink_counter.py -d $DIR"/"$A"/"$BEFORE

python cufflink_identifier.py -d $DIR"/"$A"/"$BEFORE
python cufflink_identifier.py -d $DIR"/"$A"/"$AFTER
#python cufflink_identifier.py -d $DIR"/"$B"/"$BEFORE
#python cufflink_identifier.py -d $DIR"/"$B"/"$AFTER

python cufflink_counter.py -d $DIR"/"$A"/"$BEFORE
python cufflink_counter.py -d $DIR"/"$A"/"$AFTER
