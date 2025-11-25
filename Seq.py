import sys
import custom_io as io
from Fragment import *

# Contains information about a protein sequence, base is name and sequence
# Can be constructed with a file name or a UniProt accession number

class Seq:
    def __init__(self,seq = "", name = "",file = None,UniProt_acc = None):
        
        # Access UniProt, if accession number is provided
        if UniProt_acc != None:
            self.access_uniprot(UniProt_acc)
            
        # Read file, if provided
        elif file !=None:
            self.read_file(file)
            
        # if no accession given, constructs class based on arguments
        else:     
            self.seq = str(seq)
            self.name = name
            
        # For use later, once the sequence is cleaved by an enzyme
        self.fragments = []
        
    # Used to construct class basedd on a file sequence
    def read_file(self, file):
        l = io.read_file(file)
        self.seq = l[0]
        self.name = l[1]
    
    # Used to construct class basedd on a file sequence
    def access_uniprot(self, acc):
        up = io.access_uniprot(acc)
        self.seq = up[0]
        self.name = up[1]
        
    # String represention of the Seq class, given in FASTA format
    def __str__(self):
        s = ">" + self.name + "\n" + self.seq
        return s
    
    def __repr__(self):
        return f"Seq({self.name},{self.seq})"
    
    # Will return the length of the sequence
    def __len__(self):
        return len(self.seq)
    

    
    # For error checking. Gets the length of all the fragments in the seq object
    def get_all_fragment_lengths(self):
        n = 0
        for frag in self.fragments:
            n = n + len(frag.seq)
        
        return n
    