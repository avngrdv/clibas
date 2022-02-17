# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 17:59:16 2020
@author: Alex Vinogradov
"""

import numpy as np
def dense_morgan(r, monomer_SMILES, w=None):
    '''
    Create a feature matrix dims = (number of aas, number of features)
    For each amino acid looked up in constants, create a list of Morgan
    fingerprints. Number of features = number of unique fingerprints.
    
    bit_features and info output values can be ignored in most cases
    (used primarily for mapping integrated gradient attributions to 
     the chemical structure of the underlying peptide)
    
    Parameters
    ----------
    r : int; maximum fingerprint radius
    w : flag True, if the resulting matrix should written to an .npy file

    Returns
    -------
    F : feature matrix
    bit_features : fingerprint accesion number (internal RDkit repr)
    info : list of dicts; fingerprint information (internal RDkit repr)
    '''
    from rdkit.Chem import AllChem
    from rdkit import Chem
    aas = [Chem.MolFromSmiles(x) for x in monomer_SMILES]
    
    #construct a list of all bit features
    bit_features = []
    for aa in aas:
        fingerprints = AllChem.GetMorganFingerprint(aa, r)
        keys = list(fingerprints.GetNonzeroElements().keys())
        for k in keys:
            bit_features.append(k)
            
    bit_features = list(set(bit_features))
        
    #assemble the F matrix, encoding fingerprints as a dense bit string
    F = np.zeros((len(monomer_SMILES), len(bit_features)))
    info = []
    for i, aa in enumerate(aas):
        fp_info = {}
        fingerprints = AllChem.GetMorganFingerprint(aa, r, bitInfo=fp_info).GetNonzeroElements()
        for f in fingerprints:
            F[i,bit_features.index(f)] = 1

        info.append(fp_info)

    if w is not None:
        np.save(w, F)

    return F, bit_features, info

def bit_morgan(r, n, monomer_SMILES, w=None):
    '''
    Deprecated. Unnecessarily bloats the data and has a high
    chance of feature collision unless is n is unreasonably large.
    
    Parameters
    ----------
    r : int; atomic radius of the fingerprint.
    n : int; bit length.
    w : a flag to write down the matrix as an .npy file; bool

    Returns
    -------
    Full Morgan fingerprint F matrix; ready for featurization.

    '''
    from rdkit.Chem import AllChem
    from rdkit import Chem
    
    #load amino acid SMILES
    aas = [Chem.MolFromSmiles(x) for x in monomer_SMILES]
    
    #generate an empty F matrix and fill it
    F = np.zeros((len(aas), n))
    for i,aa in enumerate(aas):
        bit = AllChem.GetMorganFingerprintAsBitVect(aa, r, nBits=n)
        F[i] = list(bit.ToBitString())
    
    #write the file if need be
    if w is not None:
        np.save(w, F)

    return F

def _get_SMILES_chars(monomer_SMILES):
    '''
    Generate a dictionary enumerating all unicode characters
    found in the SMILES tuple.
    
    return:     a list containing each unique character in
                SMILES tuple. sorted.
                
                a dictionary, which has unique characters as
                keys and their value + 1 as values.
    '''
    ch = ''.join(x for x in monomer_SMILES)
    ch = set(list(ch))
    ch = list(ch)
    ch.sort()
    return ch, {c: i+1 for i,c in enumerate(ch)}


def SMILES_repr_v1(monomer_SMILES, w=None):
    '''
    Deprecated. In this representation, 
    SMILES characters are represented as one-hot vectors. 
    The representation results in humongous dimensions.
    '''
    
    SMILES_ch, _ = _get_SMILES_chars(monomer_SMILES)
    dim = (len(SMILES_ch),)
    one_hot_SMILES_ch = np.diagflat(np.ones(dim, dtype=int))
    
    x_dim = len(monomer_SMILES)
    y_dim = max([len(x) for x in monomer_SMILES]) * dim[0]

    F = np.zeros((x_dim, y_dim)) - 1
    for i,ch in enumerate(monomer_SMILES):
        ch_as_list = list(ch)
        row = np.array([one_hot_SMILES_ch[SMILES_ch.index(x)] for x in ch_as_list]).flatten()
        F[i,:row.size] = row

    if w is not None:
        np.save(w, F)

    return F
    
def SMILES_repr_v2(monomer_SMILES, w=None):
    '''
    Deprecated. In this representation, each character of a SMILES 
    string is assigned a number (according to its rank that can be 
    looked up in get_SMILES_chars()). Representations are right-padded
    to the longest one. This representation is usually used together
    with embedding during model training.
    '''
    
    _, SMILES_d = _get_SMILES_chars(monomer_SMILES)
    x_dim = len(monomer_SMILES)
    y_dim = max([len(x) for x in monomer_SMILES])
    
    F = np.zeros((x_dim, y_dim)) - 1
    for i,ch in enumerate(monomer_SMILES):
        ch_as_list = list(ch)
        F[i,:len(ch_as_list)] = [SMILES_d[x] for x in ch_as_list]
    
    if w is not None:
        np.save(w, F)

    return F

    
def varimax_repr(alphabet, w=None):
    '''
    See J. Comput. Biol. (2009) 16, 703-723
        doi: 10.1089/cmb.2008.0173
        
    Deprecated. Represent each amino as a varimax vector.
    dimensions of F: number of aas x number of varimax features
    '''
    
    varimax_1 = {'A': 0.57, 
                 'C': 2.66,
                 'D': -2.46, 
                 'E': -3.08, 
                 'F': 3.12, 
                 'G': 0.15, 
                 'H': -0.39,
                 'I': 3.10, 
                 'K': -3.89,
                 'L': 2.72,
                 'M': 1.89, 
                 'N': -2.02, 
                 'P': -0.58, 
                 'Q': -2.54, 
                 'R': -2.80,
                 'S': -1.10, 
                 'T': -0.65, 
                 'V': 2.64, 
                 'W': 1.89, 
                 'Y': 0.79}

    varimax_2 = {'A': 3.37, 
                 'C': -1.52,
                 'D': -0.66, 
                 'E': 3.45, 
                 'F': 0.68, 
                 'G': -3.49, 
                 'H': 1.00,
                 'I': 0.37, 
                 'K': 1.47,
                 'L': 1.88,
                 'M': 3.88, 
                 'N': -1.92, 
                 'P': -4.33, 
                 'Q': 1.82, 
                 'R': 0.31,
                 'S': -2.05, 
                 'T': -1.60, 
                 'V': 0.03, 
                 'W': -0.09, 
                 'Y': -2.62}
    
    varimax_3 = {'A': -3.66, 
                 'C': -3.29,
                 'D': -0.57, 
                 'E': 0.05, 
                 'F': 2.40, 
                 'G': -2.97, 
                 'H': -0.63,
                 'I': 0.26, 
                 'K': 1.95,
                 'L': 1.92,
                 'M': -1.57, 
                 'N': 0.04, 
                 'P': -0.02, 
                 'Q': -0.82, 
                 'R': 2.84,
                 'S': -2.19, 
                 'T': -1.39, 
                 'V': -0.67, 
                 'W': 4.21, 
                 'Y': 4.11}
    
    varimax_4 = {'A': 2.34, 
                 'C': -3.77,
                 'D': 0.14, 
                 'E': 0.62, 
                 'F': -0.35, 
                 'G': 2.06, 
                 'H': -3.49,
                 'I': 1.04, 
                 'K': 1.17,
                 'L': 5.33,
                 'M': -3.58, 
                 'N': -0.65, 
                 'P': -0.21, 
                 'Q': -1.85, 
                 'R': 0.25,
                 'S': 1.36, 
                 'T': 0.63, 
                 'V': 2.34, 
                 'W': -2.77, 
                 'Y': -0.63}
    
    varimax_5 = {'A': -1.07, 
                 'C': 2.96,
                 'D': 0.75, 
                 'E': -0.49, 
                 'F': -0.88, 
                 'G': 0.70, 
                 'H': 0.05,
                 'I': -0.05, 
                 'K': 0.53,
                 'L': 0.08,
                 'M': -2.55, 
                 'N': 1.61, 
                 'P': -8.31, 
                 'Q': 0.09, 
                 'R': 0.20,
                 'S': 1.78, 
                 'T': 1.35, 
                 'V': 0.64, 
                 'W': 0.72, 
                 'Y': 1.89}
    
    varimax_6 = {'A': -0.40, 
                 'C': -2.23,
                 'D': 0.24, 
                 'E': 0.00, 
                 'F': 1.62, 
                 'G': 7.47, 
                 'H': 0.41,
                 'I': -1.18, 
                 'K': 0.10,
                 'L': 0.09,
                 'M': 2.07, 
                 'N': 2.08, 
                 'P': -1.82, 
                 'Q': 0.60, 
                 'R': -0.37,
                 'S': -3.36, 
                 'T': -2.45, 
                 'V': -2.01, 
                 'W': 0.86, 
                 'Y': -0.53}
    
    
    varimax_7 = {'A': 1.23, 
                 'C': 0.44,
                 'D': -5.15, 
                 'E': -5.66, 
                 'F': -0.15, 
                 'G': 0.41, 
                 'H': 1.61,
                 'I': -0.21, 
                 'K': 4.01,
                 'L': 0.27,
                 'M': 0.84, 
                 'N': 0.40, 
                 'P': -0.12, 
                 'Q': 0.25, 
                 'R': 3.81,
                 'S': 1.39, 
                 'T': -0.65, 
                 'V': -0.33, 
                 'W': -1.07, 
                 'Y': -1.30}
    
    varimax_8 = {'A': -2.32, 
                 'C': -3.49,
                 'D': -1.17, 
                 'E': -0.11, 
                 'F': -0.41, 
                 'G': 1.62, 
                 'H': -0.60,
                 'I': 3.45, 
                 'K': -0.01,
                 'L': -4.06,
                 'M': 1.85, 
                 'N': -2.47, 
                 'P': -1.18, 
                 'Q': 2.11, 
                 'R': 0.98,
                 'S': -1.21, 
                 'T': 3.43, 
                 'V': 3.93, 
                 'W': -1.66, 
                 'Y': 1.31}
    
 
    varimax = (varimax_1, 
               varimax_2,
               varimax_3, 
               varimax_4,
               varimax_5, 
               varimax_6,
               varimax_7, 
               varimax_8)    
    
    x_dim = len(alphabet)
    y_dim = len(varimax)
    F = np.zeros((x_dim, y_dim))
    for i,v in enumerate(varimax):
        for j,aa in enumerate(alphabet):
            F[j,i] = v[aa]

    if w is not None:
        np.save(w, F)

    return F



