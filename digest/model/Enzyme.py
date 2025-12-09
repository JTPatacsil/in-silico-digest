import sys
from .Seq import *
from .Fragment import *

# Enzyme class - contains information relevant to protein digestion
# Proteins are written N -> C
#   0 is N, i is C
# cleavage_terminus - where is the splice site. N or C
#   Indicates whether or not the cleavage peptide is the last or first of a new sequence
# Cleavage_sites are where the enzyme cleaves
#   ex: after R,K
# exception sites: if any 
    # ex: after P
# example:
# trypsin = Enzyme("C", ["K","R"], [("K","P"),("R","P")])


class Enzyme:
    def __init__(self, cleavage_terminus = "", cleavage_sites = [], exception_sites = [(None,None)]):
        
        # Will allow for N, C, or an empty string
        if cleavage_terminus != "N" and cleavage_terminus != "C" and cleavage_terminus != "":
            print("Invalid cleavage terminus given for enzyme")
            sys.exit(1)
            
        
        self.cleave_terminus = str(cleavage_terminus)
        self.cleave_sites = cleavage_sites
        self.exception_sites = exception_sites
        
        # for use later, indexes of where the enzyme will cleave
        self.cleave_indexes = []

        
        
    # This will find the cleave indexes for a given sequence
    # The cleave indexes are adjusted for whether the enzyme cleaves at the N or the C terminus
    # C terminus have indexes of i+1, to cleave after the specified amino acid
    # N terminus have indexes of i to cleave before the spcified amino acid
    # Note: no error checking to make sure the sequence given is of a Seq class
    def find_cleave_indexes(self, Seq):

        if self.cleave_terminus == "N":
            exc_sites = self.exception_sites
            exc_sites = [x[0] for x in exc_sites]
        
        elif self.cleave_terminus == "C":
            exc_sites = self.exception_sites
            exc_sites = [x[1] for x in exc_sites]
        
        sites = [0]
        
        
        for i,aa in enumerate(Seq.seq):
            if aa in self.cleave_sites:
                if self.cleave_terminus == "N" and i != 0: # Beginning of sequence
                    if Seq.seq[i-1] in exc_sites:
                        continue
                    
                if self.cleave_terminus == "C" and i !=len(Seq.seq)-1: # end of sequence
                    if Seq.seq[i+1] in exc_sites:
                        continue
                
                sites.append(i)
            
        if self.cleave_terminus == "C":
            sites = [0] + [x+1 for x in sites[1:]]
        
        sites.append(len(Seq.seq)-1)
        
        # removes duplicates if present and arranges them from N -> C
        sites = sorted(list(set(sites)))
        
        return sites
    
    
    # Will cleave a sequence, based on an enzyme. Generates fragments
    # Fragments are stored in the given Seq object
    # Note: no error checking to make sure the sequence given is of a Seq class
    #   could be fixed by checking in the find.cleave_indexes function
    
    def cleave(self, Seq):
        # Will find the indexes of a cleavage sites of a sequence
        # if none is found, expect [0,len(Seq.seq)]
        cleave_indexes = self.find_cleave_indexes(Seq)
        fragments = []
        
        for i in range(0,len(cleave_indexes)-1):
            first = cleave_indexes[i]; last = cleave_indexes[i+1]
            seq = Seq.seq[first:last]
            
            prev_aa = Seq.seq[last-2] # grabs the aa before the site
            
            # grabs the aa after the site
            try:
                next_aa = Seq.seq[last]
            except IndexError: # if at the end of the sequence
                next_aa = None
                
            frag = Fragment(seq, prev_aa, 
                            next_aa)
            
            frag.pos = cleave_indexes[i]
            
            fragments.append(frag)
            
        Seq.fragments = fragments
    
### Child classes of Enzyme
### Just different types of enzymes that can be used for the digest
class Trypsin(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus="C",
            cleavage_sites=["K", "R"],
            exception_sites= [("K","P"),("R","P")]
        )

class Lys_c(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus="C",
            cleavage_sites = ["K"],
            exception_sites= [("K","P"),("R","P")]
        )

class Lys_n(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus="N",
            cleavage_sites = ["K"],
            exception_sites= [("K","P"),("R","P")]
        )

class CNBr(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus= "C",
            cleavage_sites = ["M"]
        )

class AspN(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus= "N",
            cleavage_sites= ["D"]
            )

class ArgC(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus = "C",
            cleavage_site = ["R"],
            exception_sites = [("R","P")]
            )

class Pepsin_1_3(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus= "C",
            cleavage_site = ["F","L"]
            )

class Pepsin_gt2(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus = "C",
            cleavage_sites = ["F","L","W","Y","A","E","Q"]
            )

class PtKinase_K(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus = "C",
            cleavage_sites = ["A","F","Y","W","I","L"]
            )
