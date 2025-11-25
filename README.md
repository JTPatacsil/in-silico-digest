# in-silico-digest
In silico digest project for BINF5240 (Bioinformatics Computing at Georgetown University)


The goal of this project is to perform an enzymatic digest of a user-provided protein sequence. The user will interact with the program via a TurboGears web-interface and will be able to specify different parameters (specified below) for the digest. The output will be a table of peptides with their legnth, molecular weight, # of missed cleavages, and the amino-acids to left and right of each peptide in the protein sequence

### Parameters
* min and max length of the fragment
* min and max weight of the fragment
* number of missed cleavages allowed
* enzyme to perform the digest. Implemented enzymes are listed below:
    * Trypsin
    * Lys C
    * Lys N
    * CNBr (Cyanogen bromide)
    * AspN
    * ArgC
    * Pepsin (pH = 1.3)
    * Pepsin (pH > 2)
    * Proteinase K