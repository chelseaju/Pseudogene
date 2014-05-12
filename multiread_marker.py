"""
Usage: python multiread_marker.py -d directory -p paired
Input: -d the directory and prefix of output files
       -p paired end file [True|False] 
Output: A list of read that mapped to multiple places
Function: Go through the bam file, and count the read appearance
    output the read name is the count is > 2 for paired end or > 1 for single end
Author: Chelsea Ju
Date: 2014-04-05
"""

import pysam, argparse, datetime


"""
    Function : helper function to output message with time
"""
def echo(msg):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(msg))


def read_counter(bamfile):

    counter = {}
    bamFH = pysam.Samfile(bamfile, 'rb')
    for read in bamFH:
        name = read.qname

        if(counter.has_key(name)):
            counter[name] += 1
        else:
            counter[name] = 1
    bamFH.close()

    return counter

def export_data(counts, outfile, paired):

    threshould = 2
    if(not paired):
        threshould = 1

    fh = open(outfile, 'w')
    for k in counts.keys():
        if(counts[k] > threshould):
            fh.write("%s\n" %(k))
    fh.close()
    echo("Writing Multireads to File : %s" %(outfile))

def main(parser):
    
    options = parser.parse_args()
    dir = options.dir
    paired = options.pair
    
    ## check dir
    if(dir[-1] != "/"):
        dir += "/"

    ## start counting the readd
    bam = dir + "accepted_hits.bam"
    count_hash = read_counter(bam)
      
    ## output data
    outfile = dir + "correction/multireads.txt"
    export_data(count_hash, outfile, paired)


if __name__ == "__main__":   
   
    parser = argparse.ArgumentParser(prog='multiread_marker.py')
    parser.add_argument("-d", "--directory", dest="dir", type=str, help="directory of input files", required = True)
    parser.add_argument("-p", "--paired", dest="pair", type=bool, help="paired end read [True|False]", required = True)

    main(parser)
