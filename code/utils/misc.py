# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 12:39:30 2020
@author: Alex Vinogradov
"""
import numpy as np

def get_freqs(arr, tokens):
    #Create a matrix of zeros and iteratively fill it
    M = np.zeros((len(tokens), arr.shape[1]))
    
    for i, x in enumerate(tokens):
        M[i] = np.sum(arr == x, axis=0)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        freq = np.divide(M, arr.shape[0])
        
    return freq

def hamming_distance(P, pep, 
                     h=0,
                     cum=False,
                     return_count=False, 
                     return_index=False,
                     return_distance=False):
    '''
    h: desired hamming distance
    cum: cumulative h. i.e. cum_h=3 is h=0, 1, 2, 3
    '''

    D = P == pep
    
    if return_distance:
        return np.sum(~D, axis=1)
    
    match = pep.size - h
    if cum:
        ind = np.sum(D, axis=1) >= match
    else:
        ind = np.sum(D, axis=1) == match
        
    H = P[ind]
    
    if return_count:
        return H.shape[0]
        
    elif return_index:
        return np.where(ind)[0]

    return H


def shannon_entropy(P, norm=True, return_counts=True):
    '''
    Compute Shannon entropy for a peptide dataset.
    Note that unless norm is set to True, the 
    resulting value scales with the dataset size.
    log2 entropy computation is used.

    Parameters
    ----------
    P: P-matrix, any representation is fine
    norm : bool; if set to True, will calculate 
           "normalized entropy" (aka efficiency)

    Returns
    -------
    (Normalized) Shannon Entropy, float32

    '''
    
    #C - counts; n - dataset size
    C = np.unique(P, return_counts=True, axis=0)[1]
    n = C.sum()
    normC = np.divide(C, n)
    
    #E - a vector of individual entropy values
    E = -normC * np.log2(normC)
    if norm == True:
        E = np.divide(E.sum(), np.log2(n))
    
    else:
        E = E.sum()
        
    if return_counts:
        return E, C
    else:
        return E