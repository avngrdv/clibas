"""
Standard config parameters to default to in case an incomplete config file is 
passed.
"""

translation_table = {
                     'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                     'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                     'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                     'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                     'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                     'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                     'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                     'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                     'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                     'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                     'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                     'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                     'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                     'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                     'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
                     'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
                    }

aas = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
       'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')        
           
aa_SMILES = {
    "A": "N[C@@H](C)C(=O)",          
    "C": "N[C@@H](CS)C(=O)",
    "D": "N[C@@H](CC(=O)O)C(=O)",
    "E": "N[C@@H](CCC(=O)O)C(=O)",
    "F": "N[C@@H](Cc1ccccc1)C(=O)",
    "G": "NCC(=O)",
    "H": "N[C@@H](Cc1c[nH]cn1)C(=O)",
    "I": "N[C@@H]([C@H](CC)C)C(=O)",
    "K": "N[C@@H](CCCCN)C(=O)",    
    "L": "N[C@@H](CC(C)C)C(=O)",
    "M": "N[C@@H](CCSC)C=O",
    "N": "N[C@@H](CC(=O)N)C(=O)",
    "P": "O=C[C@@H]1CCCN1",
    "Q": "N[C@@H](CCC(=O)N)C(=O)",
    "R": "N[C@@H](CCCNC(=N)N)C(=O)",
    "S": "N[C@@H](CO)C(=O)",
    "T": "N[C@@H]([C@H](O)C)C(=O)",
    "V": "N[C@@H](C(C)C)C(=O)",
    "W": "N[C@@H](Cc1c[nH]c2c1cccc2)C(=O)",
    "Y": "N[C@@H](Cc1ccc(O)cc1)C(=O)",   
}

bases = ('A', 'C', 'G', 'T')

base_SMILES = {
    "A": "NC1=NC=NC2=C1N=CN2[C@H]3C[C@H](O)[C@@H](COP(O)(O)=O)O3",
    "C": "NC(C=CN1[C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)=NC1=O",
    "G": "O=C1C(N=CN2[C@H]3C[C@H](O)[C@@H](COP(O)(O)=O)O3)=C2N=C(N)N1",
    "T": "O=C(C(C)=CN1[C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)NC1=O",
}

complement_table = ['ACTGNYKBD', 'TGACNRMVH']













