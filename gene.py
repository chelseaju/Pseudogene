"""
Function: Construct a gene object with its location
Created: 2013-07-31
Author: Chelsea Ju
"""
class Gene:
    
    def __init__(self, name, chr, strand, start, end):
        self.name = name
        self.chr = chr
        self.strand = strand
        self.start = int(start)
        self.end = int(end)
        self.primary = 0
        self.secondary = 0
        self.size = self.end - self.start + 1
        
    def __repr__(self):
        return repr((self.name, self.chr, self.strand, self.start, self.end, self.size))

    def set_name(self, name):
        self.name = name
        
    def set_start(self, position):
        self.start = position
    
    def set_end(self, position):
        self.end = position

    def set_size(self, length):
        exon_length = length.split(",")
        self.size = sum(int(x) for x in exon_length) 

    def add_primary_fragment(self, count):
        self.primary += count

    def add_secondary_fragment(self, count):
        self.secondary += count

    def set_primary_fragment(self, count):
        self.primary = count

    def set_secondary_fragment(self, count):
        self.secondary = count

    
    def get_fragment(self, type = 0):
        # type = 0 : combine primary and secondary
        # type = 1 : primary fragments
        # type = 2 : secondary fragments

        if(type == 0):
            return float(self.primary) + float(self.secondary)
        elif(type == 1):
            return float(self.primary)
        elif(type == 2):
            return float(self.secondary)
        else:
            print ("Invalid Fragments Type: %d" %(type))
            sys.exit(2)

        
    def get_fpkm(self,total_fragment, type = 0):
        # type = 0 : combine primary and secondary
        # type = 1 : primary fragments
        # type = 2 : secondary fragments
        try:
            self.size != 0
        except:
            print("Gene size can not be zero")
            sys.exit(2)
        
        try:
            total_fragment < 1
        except:
            print ("Invalid total fragment count")
            sys.exit(2)

        if(type == 0):
            return (float(self.primary) + float(self.secondary))* float(10**9) / (float(self.size) * float(total_fragment))
        elif(type == 1):
            return (float(self.primary) * float(10**9)) / (float(self.size) * float(total_fragment))
        elif(type == 2):
            return (float(self.secondary) * float(10**9))/ (float(self.size) * float(total_fragment))
        else:
            print ("Invalid fragments type: %d" %(type))
            sys.exit(2)