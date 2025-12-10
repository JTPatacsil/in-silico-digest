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

        # for representation:
        self.name = "Base Enzyme Class"

    def __str__(self):
        return self.name

    def disp_clv_sites(self):
        return ", ".join(self.cleave_sites)

    def disp_exc_sites(self):
        if self.cleave_terminus == "N":
            e_sites = [x[0] for x in self.exception_sites]

        elif self.cleave_terminus == "C":
            e_sites = [x[1] for x in self.exception_sites]

        e_sites = set(e_sites)
        e_sites = list(e_sites)
        if e_sites[0] == None:
            return "None"
        else:
            return "If " + ", ".join(e_sites) + " is " + self.cleave_terminus + "-term of cleave site"

        
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
        
        sites = []
        
        
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
            sites = [x+1 for x in sites]

        # doing this so that the cleave_function can find all fragments
        sites.append(0) # first aa
        sites.append(len(Seq.seq)) # last aa
        
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

        # cleave_indexes start at 0 and then go to the last cleave site
        for i in range(0,len(cleave_indexes)-1):
            first = cleave_indexes[i]
            last = cleave_indexes[i+1]

            seq = Seq.seq[first:last]

            try: # if the aa is out of range
                prev_aa = Seq.seq[last-1] # grabs the aa before the site
            except:
                prev_aa = None
            
            # grabs the aa after the site
            try:
                next_aa = Seq.seq[last]
            except IndexError: # if at the end of the sequence
                next_aa = None
                prev_aa = None
                
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
        self.name = "Trypsin"


class Lys_c(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus="C",
            cleavage_sites = ["K"]
        )
        self.name = "Lys C"


class Lys_n(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus="N",
            cleavage_sites = ["K"]
        )
        self.name = "Lys N"


class CNBr(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus= "C",
            cleavage_sites = ["M"]
        )
        self.name = "CNBr"


class AspN(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus= "N",
            cleavage_sites= ["D"]
            )
        self.name = "Asp N"

class ArgC(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus = "C",
            cleavage_sites = ["R"],
            exception_sites = [("R","P")]
            )
        self.name = "Arg C"

class Pepsin_1_3(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus= "C",
            cleavage_sites = ["F","L"]
            )
        self.name = "Pepsin (pH = 1.3)"

class Pepsin_gt2(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus = "C",
            cleavage_sites = ["F","L","W","Y","A","E","Q"]
            )
        self.name = "Pepsin (pH > 2)"

class Proteinase_K(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus = "C",
            cleavage_sites = ["A","F","Y","W","I","L"]
            )
        self.name = "Proteinase K"

class Thermolysin(Enzyme):
    def __init__(self):
        super().__init__(
            cleavage_terminus = "N",
            cleavage_sites = ['A','F','I','L','M','V'],
            exception_sites = []
            )

        e_sites = ["D","E"]
        for x in e_sites:
            for y in self.cleave_sites:
                self.exception_sites.append((x,y))
        self.name = "Thermolysin"



all_enzymes = [Trypsin(),ArgC(),AspN(),Lys_n(),Lys_c(),CNBr(),
               Proteinase_K(),Pepsin_1_3(),Pepsin_gt2(), Thermolysin()]

