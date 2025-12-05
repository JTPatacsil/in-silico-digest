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


import xml.etree.ElementTree as ET
import urllib
import sys, os
import gzip

import urllib.request
from urllib.error import HTTPError

def get_seq(cont):
    # access uniprot
    # shouldn't be a dot in the name
    if cont.find(".") == -1:
        s = access_uniprot(cont)        
    else:
        s = read_file(cont)
        
    return s



def access_uniprot(uni_protID):
    # url = "https://rest.uniprot.org/uniprotkb/P12345?format=xml"
    url = "https://rest.uniprot.org/uniprotkb/" + uni_protID + "?format=xml"
    try:
        xml = urllib.request.urlopen(url)
    except HTTPError:
        print("HTTP Error, bad accession number")
        sys.exit(1)

    name = uni_protID
    
    ns = "{http://uniprot.org/uniprot}"
    tree = ET.parse(xml)
    root = tree.getroot()

    for seq in root.iter(ns + 'sequence'):
        if seq != None:
            seq = seq.text
            
        else:
            print("No sequence element found in UniProt file", file = sys.stderr)
            sys.exit(1)
    
    return [seq,name]

# Should be able to handle a fasta, gzip, zip file
def read_file(f):
    # checking to make sure the file exists
    if os.path.exists(f):
        pass
    else:
        print("File not found")
        sys.exit(1)

    # handeling file types
    if f.endswith("gzip") or f.endswith("gz"):
        f = gzip.open(f)
    else:
        f = open(f)


    # Reading file
    name     = ""
    contents = []
    for line in f:
        # print(line)
        if line.startswith(">"):
            if len(contents) == 0:  #first entry of the file
                name = line.split(">")[1]
            elif len(contents) >= 1:
                print("Reading first FASTA entry")
                break
        
        # runs regardless of FASTA or txt file
        else:
            contents.append(line.strip())
        
    return ["".join(contents),name]    
    

def make_frag_table(fragment_list):
    header = "\t".join(["num_missed_cleavages","sequence","position_start","sequence_length",
                        "sequence_weight","left_splice_aa","right_splice_aa"])
    
    print(header)
    for f in fragment_list:
        print(f.export_as_tsv())

