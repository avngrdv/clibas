# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 17:46:19 2022
@author: Alex Vinogradov
"""

import numpy as np
import os

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

class FeatureMatrix:
    
    def __init__(self, F=None, constants=None):
        self.F = F
        self.constants = constants
        self._validate()
        return
    
    def __repr__(self):
        return f'<FeatureMatrix object>\nF={self.F}'
    
    def _validate(self):
        if self.F is not None:
            if not isinstance(self.F, np.ndarray):
                raise ValueError(f'<FeatureMatrix>: feature matrix expected as a numpy array; received: {type(self.F)}')
            else:
                if not self.F.ndim == 2:
                    raise ValueError(f'<FeatureMatrix> expected to receive parameter F as a two-dimensional array; received: ndims={self.F.ndim}') 

        return

    @classmethod
    def make(cls, descr=None, constants=None):

        if descr is None:
            return cls(F=None, constants=constants)

        elif isinstance(descr, np.ndarray):
            return cls(F=descr, constants=constants)

        elif isinstance(descr, str):

            if os.path.isfile(descr):
                try:
                    F = np.load(descr)
                except:
                    raise ValueError('<FeatureMatrix.make> could not load the matrix from the specified location. . .')

            #if all easy setups failed, resort to string parsing
            else:
                
                #specific to proteinogenic amino acids
                if descr == 'varimax':
                    F = cls.varimax_repr(constants.aas)
    
                if descr[:descr.find('_')] == 'pep':
                    SMILES = constants.aa_SMILES
                    
                elif descr[:descr.find('_')] == 'dna':
                    SMILES = constants.base_SMILES

                else:
                    raise ValueError('<FeatureMatrix.make> did not understand feature matrix specification. . .')

                if 'ECFP' in descr:
                    r = int(descr[-1])
                    F, _, _ = cls.dense_morgan(r, SMILES)
                        
                elif 'SMILES' in descr:
                    F = cls.SMILES_repr_v2(SMILES)
                
                else:
                    raise ValueError('<FeatureMatrix.make> did not understand feature matrix specification. . .')

        return cls(F=F, constants=constants)


    def save(self, fname):
        if not self.F:
            raise ValueError('No feature matrix to be saved!')
        else:
            np.save(fname, self.F)

        return
    
    def single_linkage_tree(self, labels=None, fname=None):
        
        from scipy.cluster.hierarchy import linkage
        import clibas.plotters as P
               
        if self.F is not None:
            P.Miscellaneous.single_linkage_dendrogram(
                
                        linkage(self.F, method='ward'), 
                        labels=labels,
                        basename=fname
                                                     )
        return
    
    @staticmethod
    def dense_morgan(r, monomer_SMILES):
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
    
        return F, bit_features, info

    @staticmethod
    def SMILES_repr_v2(monomer_SMILES):
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
        
        return F
    
    @staticmethod
    def varimax_repr(alphabet):
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
    
        return F
    
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


