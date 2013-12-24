""" personalized_genome.py 
Usage: python personalized_genome -g default_genome -m mutation_rate -o outprefix
Input: mutation rate and the consensus genome sequence
Output: personalized genome, and the base different from consensus
Author: Chelsea Ju
Function: Read in the consensus genome, and randomly change the base with a given mutation rate.
Date: 2013-08-05
"""


import sys, getopt, re, random, copy
 
#global variables
GENOME = dict()
GENOME_SIZE=[]
NT = ['A', 'C', 'G', 'T']

def usage():
    options = [ ("-g genome.fa, --genome = genome in fasta format", ""),
                ("-m mutation rate, --mutation = mutation rate, default = 0.001", ""),
                ("-o outprefix, --outprefix = prefix of output file name", "")
    ]
    
    print "Usage: python personalized_genome -g default_genome -m mutation_rate -o outprefix"
    
    for (o,d) in options:
        print "\t" + o + "\t" + d
    print "\n"


def extract_fasta(reference_fh):
    seq = []
    title = ""
    start = 1
    for line in reference_fh:
        line = line.rstrip()
        m = re.match('^>(.*)', line)
        if(m):
            if(len(seq) != 0):
                s = "".join(seq)
                GENOME[title] = s
                GENOME_SIZE.append((title, start, start + len(s) - 1))
                start = start + len(s)
                seq = []
            title = m.group(1)

        else:
            seq.append(line)
    
    s = "".join(seq)
    GENOME[title] = s
    GENOME_SIZE.append((title, start, start + len(s) - 1))

 
def mutate_genome(rate, snp_fh):

    (title, start, genome_size) = GENOME_SIZE[-1]
    mutation_position = random.sample(xrange(1,genome_size), int(rate * genome_size))

    snp_fh.write("Mutation Rate = %d" %(rate))
    for x in mutation_position:
        location = filter(lambda chromosome: x >= chromosome[1] and x <= chromosome[2], GENOME_SIZE)
        (chr, s, e) = location[0]
        position = x - s
        before = GENOME[chr][position].upper()
        
        if(before != 'N'):
            tmp_NT = copy.deepcopy(NT)
            tmp_NT.remove(before)
            after = random.choice(tmp_NT)        
            GENOME[chr] = GENOME[chr][0:position] + after + GENOME[chr][position+1:]
        
            snp_fh.write('%s position %d %s --> %s\n' %(chr, position, before, after))

def export_fasta(fasta_fh):
    for key in sorted(GENOME.keys()):
        fasta_fh.write(">%s\n%s\n" %(key, GENOME[key]))


    
def main():
    
    # read in user's input
    try:
        opts, args = getopt.getopt(sys.argv[1:], "g:m:o:", ["genome=", "mutation=", "outprefix=" ])
    
    except getopt.GetoptError, err:
    
    # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    for o, a in opts:
        if o in ("-g", "--genome"):
            reference_file = str(a)
        elif o in ("-o", "--outprefix"):
            outfile = str(a)
        elif o in ("-m", "--mutation"):
            mutation_rate = float(a)
        else:
            assert False, "unhandled option"
    
    try:
        reference_file != None
        outfile != None
    
    except:
        print "Missing arguments"
        usage()
        sys.exit(2)

    if(mutation_rate == None):
        mutation_rate = 0.01

    personal_genome = outfile+ ".fa"
    personal_snp = outfile+".snp"

    reference_fh = open(reference_file, 'rb')
    genome_fh = open(personal_genome, 'w')
    snp_fh = open(personal_snp, 'w')
    
    extract_fasta(reference_fh)
    mutate_genome(mutation_rate, snp_fh)
    export_fasta(genome_fh)

    reference_fh.close()
    genome_fh.close()
    snp_fh.close()
        
    print "Write personal genome to file %s" %(personal_genome)
    print "Write personal SNP to file %s" %(personal_snp)

    print "DONE\n\n"



if __name__ == "__main__":
    main()






