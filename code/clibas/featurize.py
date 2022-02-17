# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 17:46:19 2022
@author: Alex Vinogradov
"""

import numpy as np

def from_matrix_v3(X, alphabet=None, F=None, reshape=True, repad=False):
    '''
    A general featurizer tool to represent datasets (list of peptides, DNA etc)
    as featurized matrices for downstream analyses applications. The implementation
    is fairly efficient, but can be further sped up (~10x) by using cupy to
    run matrix multiplication (not implemented here)
        
        X:     numerically represented categorical data (lists of peptides, DNA)
               dimensions: #entries x sequence length 
              
        F:     Matrix containing featurized representations.
               Each row is a new representation of a token in X
               Dimensions: alphabet size x #features to factor into
               Some F matrices may be padded if tokens featurioze into 
               features of different length. Paddings should have the values of -1.
              
               If left as None, data will be represented with one hot encoding
              
 alphabet:     Current representations of X
       
    repad:     Should be flagged True when F is internally padded.
               Many representations, for instance one hot, have equally
               long vectors corresponding to each amino acids, but some,
               for instance, SMILES_repr_v2 are all different.
               In that case, the SMILES_repr_v2 matrix is internally padded
               to the longest representation, which upon mapping to X will
               result in pads in the middle of the sequence. Flagging repad
               will push all pads to the right
               
   reshape:    Should be flagged if a peptide sequence is to be represented as
               a 2D array. REPADDED MATRIX SHOULD NOT BE RESHAPED 
               (it doesn't make sense but mathematically possible)
              
       out:    X-style matrix with factored representations.
               Dimensions: #peptides x sequence length * #features for each aa
    '''

    #make sure X is int-type    
    if not np.issubdtype(X.dtype, np.integer):
        raise ValueError('<featurize.from_matrix_v3> only takes int-type arrays as inputs.')
    
    #if no tokens are specified, attempt to infer
    if alphabet is None:
        alphabet = np.unique(X)
    else:
        alphabet = np.asarray(alphabet)
    
    #check whether the array is padded
    pad = X.dtype.type()
    _is_padded = pad in X
    
    #special treatment for padded arrays: try to rescue featurization
    if _is_padded:
        
        if alphabet.dtype.type() not in alphabet:
            alphabet = np.hstack((alphabet, alphabet.dtype.type()))
        
        #make sure F and X are matmul-compatible dimension-wise
        #F may not have encoding for pad tokens: in this case, 
        #just append an empty row to the matrix.
        if F is not None:
            if F.shape[0] != alphabet.size:
                if F.shape[0] == alphabet.size - 1:
                    F = np.vstack((F, np.zeros(F.shape[1])))
                else:
                    raise ValueError('<featurize.from_matrix_v3>: F matrix and \
X array have incompatible dimensions; cannot perform featurization. . .')
        
    x_shape = X.shape
    X = X.ravel()
    expansion = alphabet.size
    
    #convert the matrix to a one-hot encoding
    fX = np.zeros((X.size, expansion))
    fX[np.arange(X.size), X] = 1
    
    #matmul by a factorization matrix to get the featurized repr
    if F is not None:   
        fX = np.matmul(fX, F)
    
    fX = np.reshape(fX, (x_shape[0], -1))
    if repad:
        #find where values of X correspond to pads (-1)
        #and then do stable argsort to move them to the right
        ind = np.argsort(fX == -1, kind='stable')
        
        #reindex the matrix
        fX = fX[np.arange(fX.shape[0])[:, None], ind]
        
        #replace -1 with 0 for output
        fX[fX == -1] = 0
    
    if reshape:
        fX = np.reshape(fX, x_shape + (-1,))
    
    return fX


#TODO: needs to be rewritten; currently, dimension inference is all sorts of unreliable
def into_h5(X, y=None, path=None, F=None, alphabet=None, reshape=False, repad=False, chunks=20):
    '''
    Featurize the matrix and write it to an hdf5 file. Mainly used when
    the featurized X doesn't fit the memory. X will be split in chunks,
    factorized chunk by chunk and written to the file.
    
        X: data (np 2D array, numerically represented) to be featurized
        
        y: labels as a 1D p.ndarray
        
        path: full path to the .hdf5 file to be created
        
        F, reshape, repad, alphabet: featurization parameters passed to the 
                                     from_matrix_v3 routine.
                           
                           
        chunks: int; number of chunks to split the X in to. Each resulting 
                chunk should be small enough to fit the memory in the featurized
                form.
                    
        out:    None
                           
    '''

    import h5py
    
    if F is None:
        z = alphabet.size
    else:
        z = F.shape[-1]
    
    if reshape:
        dims = [X.shape[0], X.shape[-1], z]
    else:
        dims = [X.shape[0], X.shape[-1] * z]
        
    with h5py.File(path, 'w') as f:
        x_set = f.create_dataset("X", dims, dtype=np.float32)
        
        #n is the number of subarrays.
        indx = np.array_split(np.arange(X.shape[0]), chunks)

        for ind in indx:
            x_set[ind] = from_matrix_v3(X[ind],
                                        F=F, 
                                        alphabet=alphabet,
                                        reshape=reshape,
                                        repad=repad
                                       )
        
        if y.ndim > 0:
            y_set = f.create_dataset("y", (y.size,), dtype=np.float32)
            y_set[...] = y   
    
    return


