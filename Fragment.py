# Sequence fragment that would be generated from an enzymatic digest
# Initialized based on a sequence, and the aa around a splice site
#   Will calculate the length and weight of a fragment
class Fragment:
    def __init__(self, seq, prev_aa, next_aa, pos = 0, name = "0"):
        self.seq = str(seq)
        self.surrounding_aa = (str(prev_aa),str(next_aa))
        self.pos = int(pos)
        self.name = str(name)
        
    def __str__(self):
        s = "Seq: " + str(self.seq) + "\nSurrounding amino acids: " + str(self.surrounding_aa)
        return s
    
    def __repr__(self):
        return f"Fragment({self.name},{self.seq},{self.surrounding_aa})"
    
    # specifies how to add two fragments together
    def __add__(self,other):
        n_seq = self.seq + other.seq
        return Fragment(n_seq, self.surrounding_aa[0], other.surrounding_aa[1],
                        pos = self.pos)
    
    def __len__(self):
        return len(self.seq)
    
    # calculates the molecular weight of a sequence
    # If an unknown amino acid is given, the function will provide the 
    # average weight of an amino acid, 110.0 Da. 
    
    def seq_weight(self):
        seq = self.seq
        
        d = {"A":89.09 , "R":174.20,"N":123.12,"D":133.10,
             "C":121.16, "E":147.13,"Q":146.14,"G":75.07,
             "H":155.15, "I":131.17,"L":131.17,"K":146.19,
             "M":149.21, "F":165.19,"P":115.13,"S":105.09,
             "T":119.12, "W":204.23,"Y":181.19,"V":117.15
             }
        
        w = 0
        for aa in seq:
            try:
                w = w + d[aa]
            except KeyError:
                # Should Unknown weights be weighted 0 or the aveage weight?? 
                #   Or the weight of just the backbone
                # For now, average weight of an amino acid
                # the weight of the backbone would be conservative, underestimating 
                w = w + 110.0
                # pass
            
        return w + 19 # takes into account the H2O and H at the ends of the chain
    
    def seq_length(self):
        return len(self.seq)
    
    # POSITION. NOT INDEX
    def frag_position(self):
        start = self.pos + 1
        end = self.pos + self.seq_length()
        
        return (start, end)
    
    # determines if a fragment is valid, based on the user's specification
    # DOES NOT check exception sites of enzymes
    # Only checks the fragment length and weight
    
    def isValidFragment(self,min_l,max_l,min_w,max_w):            
        valid = False
        
        if self.isValidLength(min_l,max_l) and self.isValidWeight(min_w, max_w):
            valid = True
            
        return valid
    
    # helper function, determines if a sequence is of valid weight
    def isValidWeight(self,min_w,max_w):
        valid = False
        if self.seq_weight() >= min_w and self.seq_weight() <= max_w:
            valid = True
            
        return valid
    
    # helper function, determines if a sequence is of valid length
    def isValidLength(self,min_l,max_l):
        valid = False
        if len(self.seq) >= min_l and len(self.seq) <= max_l:
            valid = True
            
        return valid
    
    def export_as_tsv(self):
        l = [self.name, self.seq, self.pos, self.seq_length(),round(self.seq_weight(),3),
                      self.surrounding_aa[0],self.surrounding_aa[1]]
        l = map(str, l)
        return "\t".join(l)