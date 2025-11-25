# get the molecular weight dictionary of amino avids
# https://www.researchgate.net/figure/Molecular-Formula-Molecular-Mass-N-Content-and-Calculated-CP-of-20-Amino-Acids-1_tbl1_331719969

# Should look like:
    #https://web.expasy.org/peptide_mass/

import sys
import custom_io as io
from Seq import *
from Enzyme import *
from Fragment import *
    
# digest should take into account
#   sequence and the enzyme
#   the minimum and maximum fragment lengths
#   the min and max fragment weight
#   the number of missed cleavages

def get_missed_cleavages(Seq, num_missed):
    frags = Seq.fragments
    res = []
    
    if num_missed == 0:
        for f in frags:
            f.name = "0"
            res.append(f)
        return res
    
    for i in range(0,len(frags) - num_missed):
        
        f = frags[i]        
        
        if num_missed > 0:
            for j in range(1,num_missed + 1):
                f = f + frags[i+j]
        
        f.name = str(num_missed)
        res.append(f)
        
    return res


def digest(Seq, Enzyme, min_l, max_l, min_w, max_w, missed_cleavages):
    # Generates Fragments
    
    Enzyme.cleave(Seq)
    
    # Psuedocode
    # loop through 0:missed_cleavages
    #   generate missed cleavages fragments
    #   check fragments to see if they are valid
    
    # generate all possible fragments
    all_frags = []
    for i in range(0,missed_cleavages+1):
        mc = get_missed_cleavages(Seq,i)
        all_frags.extend(mc)
    
    # valid fragments
    valid_frags = []
    
    # loop through all the fragments generated
    for f in all_frags:
        # if the fragment is valid, append to the list
        if f.isValidFragment(min_l, max_l, min_w, max_w):
            valid_frags.append(f)
        
    return valid_frags

