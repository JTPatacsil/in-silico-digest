
# Should look like:
    #https://web.expasy.org/peptide_mass/

import sys
from .Seq import *
from .Enzyme import *
from .Fragment import *
    
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
        
        # fragment name is the number of missed cleavages it has
        f.name = str(num_missed)
        res.append(f)
        
    return res


def enzyme_digest(Seq, Enzyme, min_l, max_l, min_w, max_w, missed_cleavages):
    # Generates Fragments
    
    Enzyme.cleave(Seq)
    
    # generate all possible fragments
    # we do this because we do not know which fragments are going to be generated
    # so generate all and hten filter for valid
    
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

