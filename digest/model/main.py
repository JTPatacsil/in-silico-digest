# Protein Digest
# Write a simple web-server application using TurboGears to carry out 
#   an in silico enzymatic digest of a user-provided protein sequence.

# user should specify:
    # min and max length
    # min and max molec weight
    # Number of missed cleavages
    # specific enzyme
    
# Outout a table of peptides with:
    # peptide length
    # molecular weight
    # Number of missed cleavages
    # amino acids left and right of each peptide in the proteins seq

import sys, os
# Add parent directory to path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

from . import custom_io as io
from digest import *


#### Importing data into program #####
try:
    # Sequence
    seq_input = sys.argv[1]

    # Length
    min_length = int(sys.argv[2])
    max_length = int(sys.argv[3])

    # Molecular weight
    min_weight = float(sys.argv[4])
    max_weight = float(sys.argv[5])

    #number of missed cleavages
    missed_cleavages = int(sys.argv[6])

    # Specific enzyme
    enzyme = sys.argv[7]

except IndexError:
    print("Insufficient arguments given")
    print("Please give: sequence, min length, max length, min weight, max weight, number of missed cleavages, enzyme")
        
    # sys.exit(1)

    enzyme = "Trypsin"
    seq_input = "P12345"
    min_l = 0
    max_l = 999
    min_w = 0
    max_w = 1500000
    missed_cleavages = 0


#### Logical error checking
if min_l > max_l:
    print("Minimum peptide length can't be larger than maximum peptide length",
          file = sys.stderr)
    sys.exit(1)
    
if min_w > max_w:
    print("Minimum molecular length can't be larger than maximum molecular length",
          file = sys.stderr)
    sys.exit(1)
    

#### Setting up program
if seq_input.find(".") != -1: # If a "." is in a file name:
    try:
        seq = io.read_file(seq_input)
    except FileNotFoundError:
        print("File not found on OS")
else:
    seq = io.access_uniprot(seq_input)

# setting enzyme to enzyme class
if enzyme == None or enzyme == "Trypsin":
    enzyme = Trypsin()
elif enzyme == "Arg C":
    enzyme = ArgC()
elif enzyme == "Asp N":
    enzyme = AspN()
elif enzyme == "Lys N":
    enzyme = Lys_n()
elif enzyme == "Lys C":
    enzyme = Lys_c()
elif enzyme == "CNBr":
    enzyme = CNBr()
elif enzyme == "Protein Kinase K":
    enzyme = PtKinase_K()
elif enzyme == "Pepsin (pH 1.3)":
    enzyme = Pepsin_1_3()
elif enzyme == "Pepsin (pH > 2)":
    enzyme = Pepsin_gt2()

## Running program
seq = Seq(UniProt_acc = seq_input)
enzyme_digest(seq, enzyme, min_l,max_l,min_w,max_w,missed_cleavages)

io.make_frag_table(seq.valid_fragments)