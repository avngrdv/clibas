# -*- coding: utf-8 -*-

"""
Created on Tue Aug  3 17:45:24 2021
@author: Alex Vinogradov
"""

import numpy as np
class SequencingSample:
    '''
    See below: Data object.
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
        dim = []
        for arr in (self.D, self.Q, self.P):
            try:
                L.append(arr.shape[0])
            except:
                pass
            
            if arr.ndim:
                dim.append(arr.ndim)
                
        if not len(set(L)) == 1:
            raise ValueError('Sequencing sample datasets contain data of different length!')

        if not len(set(dim)) == 1:
            raise ValueError('Sequencing sample datasets contain data of different dimensionality!')

        self._internal_state = np.array(None)
        self._is_collapsed = False
        return
    
    #if a sample is iterated over, 
    #it will only expose D, Q and P attributes
    def __iter__(self):
        for tup in ((self.D, 'DNA'), (self.Q, 'Q score'), (self.P, 'Peptide')):
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
        return f'<SequencingSample {self.name} containing {len(self)} entries>'

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

    def __getitem__(self, item):
        if item == 'pep':
            return self.P
        
        elif item == 'dna':
            return self.D
        
        elif item == 'Q':
            return self.Q

        else:
            raise ValueError(f'{item} indexing not understood by SequencingSample {self.name}')
        
        return
    
    def __setitem__(self, key, value):
        if key == 'pep':
            self.P = value
        
        elif key == 'dna':
            self.D = value
        
        elif key == 'Q':
            self.Q = value

        else:
            raise ValueError(f'{key} indexing not understood by SequencingSample {self.name}')
        
        return    
        
    def _collapse_internal_state(self):
        '''
        the internal state array may need to be collapsed sometimes.
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

    def get_ndims(self):
        '''
        Return the dimensionality of (P, Q, D) datasets
        '''
        dim = []
        for arr in (self.D, self.Q, self.P):
            if arr.ndim:
                dim.append(arr.ndim)
        
        if not len(dim):
            return 0
        
        if not len(set(dim)) == 1:
            raise ValueError('Sequencing sample datasets contain data of different dimensionality!')
        
        return dim[0]

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
    During analysis, data is stored as a Data object instance. Data is just a
    container for individual samples, which are stored as SequencingSample objects.
    Any number of DNA sequences can be a sample in principle, but in practice, 
    most of the time one sample = a single .fastq file. SequencingSample objects 
    have four public attributes: 
    
    	SequencingSample.name: sample name (as a str)
        SequencingSample.D: a list of DNA sequences (can be set as None)
    	SequencingSample.Q: a list of Q score sequences (can be set as None)
    	SequencingSample.P: a list of peptide sequences (can be set as None)
    	
    These lists are stored as numpy arrays: 1D prior to calling 
    FastqParser.transform() or FastqParser.translate(), and 2D arrays the 
    entire time after that; shape: (number of entries, sequence length).
    Because the sequences for different reads may have a different length,
    arrays are padded to the longest sequence.
    
    The number of entries in each array is kept equal throughout the process 
    unless one or more of the attributes are set to None. Although any given 
    filtration routine [for example, FastqParser.q_score_filt()] acts on a single
    array [SequencingSample.Q in this example], the entries for all three arrays
    are discarded/kept as a result.
    
    Depending on how many and what kinds of templates are specified in 
    LibraryDesign, any given entry in SequencingSample may in principle be 
    compatible with several templates simultaneously. Figuring out what entry
    should be assigned to what kind of template is one of the primary objectives
    of the parser. Initially, [i.e. right after calling FastqParser.translate()] 
    the parser deems every sequence to be compatible with every specified 
    template. As filtration goes on, op by op, this compatibility is refined.
    I call the state of the assignment for a particular SequencingSample
    _sample’s internal state_. Some ops [for instance, FastqParser.fetch_at()] 
    need to know exactly which template should be associated with which entry; 
    if they find an entry that is compatible with multiple possible templates, 
    they will “collapse” sample’s internal state, by choosing one compatible
    template and assigning everything else as incompatible. Refer to the list
    of ops below for details on which ops can collapse sample’s internal state. 
    In general, these should be called after filtration ops.
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
