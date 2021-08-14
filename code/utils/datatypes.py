# -*- coding: utf-8 -*-

"""
Created on Tue Aug  3 17:45:24 2021
@author: Alex Vinogradov
"""

import numpy as np
class SequencingSample:
    '''
    A bundle of datasets (DNA, Q scores, peptides)
    for a sinle sequencing sample. Enforces everything
    to be stored in numpy arrays.
    '''
    def __init__(self, D=None, Q=None, P=None, name='unnamed'):
        self.D = D
        self.Q = Q
        self.P = P
        self.name = name
        self._validate()
        return
    
    def _validate(self):
        
        if not isinstance(self.name, str):
            raise ValueError('SequencingSample name was not understood. . .')
              
        self.D = np.array(self.D)
        self.Q = np.array(self.Q)
        self.P = np.array(self.P)    
        
        L = []
        for arr in (self.D, self.Q, self.P):
            try:
                L.append(arr.shape[0])
            except:
                pass
        
        if not len(set(L)) == 1:
            raise ValueError('Sequencing sample datasets contain data of different length!')

        self._internal_state = np.array(None)
        self._is_collapsed = False
        return
    
    #if a sample is iterated over, 
    #it will only expose D, Q and P attributes
    def __iter__(self):
        for tup in ((self.D, 'DNA'), (self.Q, 'Q scores'), (self.P, 'P')):
            yield tup

    def __len__(self):

        L = []
        
        for arr in (self.D, self.Q, self.P):
            try:
                L.append(arr.shape[0])
            except:
                pass
        
        if not len(set(L)) == 1:
            raise ValueError('Sequencing sample datasets contain data of different length!')

        return L[0]

    def __repr__(self):
        return f'<SequencingSample object containing {len(self)} entries>'

    def __call__(self, ind):
        '''
        Fancy indexing for datasets altogether.
        Re-indexing (filtration) is done in place.
        Not sure whether this is a best idea, but 
        it is easy to rewrite this.
        '''
        
        if self.D.ndim > 0:          
            self.D = self.D[ind]

        if self.Q.ndim > 0:          
            self.Q = self.Q[ind]

        if self.P.ndim > 0:          
            self.P = self.P[ind]
            
        if self._internal_state.ndim > 0:
            self._internal_state = self._internal_state[ind]
            
        return

    def _collapse_internal_state(self):
        '''
        the internal state array may need to be collapsed.
        what this means: some entries may perfectly match 
        at some stage multiple library designs. 
        if we are taking variable region number 1, and
        it happens to differ between two designs, and
        a particular entry matches both designs, which
        design do we go for? answer: we collapse the internal
        state to a single match per design and then go for it
        '''
        if self._internal_state.ndim > 0:
            collapsed = np.zeros_like(self._internal_state, dtype=np.bool)
            ind = np.arange(self._internal_state.shape[0]), self._internal_state.argmax(axis=1)
            collapsed[ind] = self._internal_state[ind]
            self._internal_state = collapsed
            self._is_collapsed = True
        return

    def transform_P(self, pad=None):
        '''
        An exposed method to transform 1D arrays into 2D.
        Upon transformation, ndarray's dtype is transformed
        to '<U1'. If 1D array contains sequences of unequal
        length, padding is necessary. If 'pad' is not specified, 
        the resulting 2D array will be padded with empty strings ('').
        '''
        
        if self.P.ndim == 1:
            shape = (self.P.size, max([len(x) for x in self.P]))
            P2d = np.zeros(shape, dtype='<U1')
            
            if pad is not None:
                P2d.fill(pad)
            
            for i, pep in enumerate(self.P):
                P2d[i,:len(pep)] = list(pep) 
        
            self.P = P2d
            
        return                

    def transform_D(self, pad=None):
        '''
        An exposed method to transform 1D arrays into 2D.
        Upon transformation, ndarray's dtype is transformed
        to '<U1'. If 1D array contains sequences of unequal
        length, padding is necessary. If 'pad' is not specified, 
        the resulting 2D array will be padded with empty strings ('').
        '''
        
        if self.D.ndim == 1:
            shape = (self.D.size, max([len(x) for x in self.D]))
            D2d = np.zeros(shape, dtype='<U1')
            
            if pad is not None:
                D2d.fill(pad)
            
            for i, read in enumerate(self.D):
                D2d[i,:len(read)] = list(read) 
        
            self.D = D2d
            
        return 

    def transform_Q(self, pad=None):
        '''
        An exposed method to transform 1D arrays into 2D.
        Upon transformation, ndarray's dtype is transformed
        to int16. If 1D array contains sequences of unequal
        length, padding is necessary. If 'pad' is not specified, 
        the resulting 2D array will be padded with empty strings ('').
        
        Converts Q datasets to numerical representations stored
        as 2D ndarrays.
        '''
        
        def ord_mapper(x):
            return [ord(a) - 33 for a in list(x)]
        
        if self.Q.ndim == 1:
            shape = (self.Q.size, max([len(x) for x in self.Q]))
            Q2d = np.zeros(shape, dtype=np.int16)
            
            if pad is not None:
                Q2d.fill(pad)
            
            for i, read in enumerate(self.Q):
                Q2d[i,:len(read)] = ord_mapper(read) 
        
            self.Q = Q2d
            
        return 

    def transform(self, pad=None):
        '''
        Perform all three (P, Q, D) transforms together.
        '''
        self.transform_D(pad=pad)
        self.transform_Q(pad=pad)
        self.transform_P(pad=pad)
        return

    def drop(self, where=None):
        '''
        Drop a dataset. Not sure if needed, but for a good measure.
        '''
        if where == 'pep':
            self.P = np.array(None)
            
        elif where == 'dna':
            self.D = np.array(None)
            
        elif where == 'q':
            self.Q = np.array(None)
            
        else:
            raise ValueError(f'SequencingSample does not have {where} attribute to drop. . .')
            
        return

    def unpad(self):
        
        def _unpad(arr):        
            
            if arr.ndim == 2:
                pad = np.zeros(1, arr.dtype)[0]
                return arr[:, ~np.all(arr == pad, axis=0)]
            else:
                return arr
        
        self.P = _unpad(self.P)
        self.D = _unpad(self.D)
        self.Q = _unpad(self.Q)
        return
    
class Data:
    '''
    Data is just a bunch of samples; all sorts of samples
    can all be held in a single container to facilitate
    pipeline creation.
    
    NOTE: what's really annoying is that when both TrainingSample
    and SequencingSample contain only optional fields, Union will
    recast everything into whatever type is specified first. For
    this reason, X in TrainingSample is not optional.
    '''
    def __init__(self, samples=None):
        self.samples = samples
        self._validate()
        return
    
    def _validate(self):
        if not isinstance(self.samples, list):
            msg = f'Data object expected to receive a list of samples; received {type(self.samples)} instead'
            raise ValueError(msg)

        for sample in self.samples:
            if not isinstance(sample, SequencingSample):
                msg = f'Unexpected sample types passed to the Data object; type: {type(sample)}'
                raise ValueError(msg)
        
        return

    def __iter__(self):
        for sample in self.samples:
            yield sample

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, item):
        return self.samples[item]

    def __repr__(self):
        return f'<Data container holding {len(self)} samples>'
