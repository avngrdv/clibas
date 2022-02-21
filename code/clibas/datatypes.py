# -*- coding: utf-8 -*-

"""
Created on Tue Aug  3 17:45:24 2021
@author: Alex Vinogradov
"""
import numpy as np

class Sample:
    '''
    Base class for various data sample types. Should not be invoked directly
    '''
    def __init__(self, name='unnamed_sample', **kwargs):
        
        for key, value in kwargs.items():
            if key != 'name':
                self.__dict__.update({key: np.asarray(value)})           
    
        self.name = name
        return

    def __len__(self):

        L = []
        for arr in self._arrays():
            try:
                L.append(arr.shape[0])
            except:
                pass
            
        assert len(set(L)) < 2, 'Sample datasets cannot contain data of different size!'
        if L:
            return L[0]
        else:
            return 0

    #TODO: maybe worth it to make this a dynamic property:
    #if item is in self.__dict__.keys() than just return
    #self.__dict__[item] 
    
    #but if item is int or array, use np fancy indexing
    #to fetch a view of the sample; for now there really
    #isn't much use for the latter
    
    def __getitem__(self, item):   
        return getattr(self, item)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)
        return

    #make sure that when new attributes are set, they don't break 
    #the dataset validity (equal size for all arrays).
    def __setattr__(self, key, value):
        if key != 'name':
            super().__setattr__(key, np.asarray(value))
        else: 
            super().__setattr__(key, value)
        self._validate_size()

    @property
    def size(self):
        return self.__len__()

    #iterate over arrays; used by pipeline methods only
    def iterarrays(self):
        
        for attr in self.__dict__:
            if type(getattr(self, attr)) == np.ndarray:
                if not attr.startswith('_'):
                    arr = getattr(self, attr)        
                    yield (arr, attr)
        
    #easy way to fetch data arrays only (np.ndarrays)
    def _arrays(self):
        
        return [getattr(self, attr) for attr in self.__dict__
                if type(getattr(self, attr)) == np.ndarray]
                    
    def _validate_name(self):
        
        if not isinstance(self.name, str):
            raise ValueError(f'Sample name {self.name} was not understood. . .')        
        return

    def _validate_size(self):
        self.__len__()
        return
    
    def _validate_dims(self):

        dim = []
        for arr in self._arrays():
            if arr.ndim:
                dim.append(arr.ndim)       

        if len(set(dim)) > 1:
            raise ValueError('Sample datasets contain data of different dimensionality!')

        if len(set(dim)) == 1:
            return dim[0]
        
        else:
            return 0

    def ind_filter(self, ind):
        '''
        Fancy indexing for datasets altogether. Re-indexing (filtration) is 
        done in place.
        '''
        for attr in self.__dict__:
            if type(getattr(self, attr)) == np.ndarray:
                arr = getattr(self, attr)
                if arr.ndim > 0:
                    #note that setattr can't be used as it would fail 
                    #the _validate_size check 
                    self.__dict__[attr] = arr[ind]            
        return

    def transform(self, item, pad=None):

        if not type(getattr(self, item)) == np.ndarray:
            raise TypeError('Unsupported array type passed to Sample.transform function.')
        
        arr = getattr(self, item)
        if arr.ndim == 1:
            
            if np.issubdtype(arr.dtype, np.character) or np.issubdtype(arr.dtype, np.object):
                shape = (arr.size, max([len(x) for x in arr]))
                arr_2d = np.zeros(shape, dtype='<U1')
                
                if pad is not None:
                    shape.fill(pad)
                
                for i, x in enumerate(arr):
                    arr_2d[i,:len(x)] = list(x) 
            
                setattr(self, item, arr_2d)               
        return

    def drop(self, item):
        if not type(getattr(self, item)) == np.ndarray:
            return
        
        setattr(self, item, np.asarray(None))
        return

    def unpad(self, item):
        
        arr = getattr(self, item)
        if type(arr) == np.ndarray:            
            if arr.ndim == 2:
                pad = arr.dtype.type()
                setattr(self, item, arr[:, ~np.all(arr == pad, axis=0)])
        return
        
class SequencingSample(Sample):
    '''
    See below: Data object.
    '''
    def __init__(self, dna=None, Q=None, pep=None, name='unnamed'):
        super(SequencingSample, self).__init__(dna=dna, Q=Q, pep=pep, name=name)
        self._validate_name()
        self._validate_size()
        self._validate_dims()
        
        self._internal_state = None
        self._is_collapsed = False
        return
    
    def __repr__(self):
        return f'<SequencingSample {self.name} containing {len(self)} entries>'
          
    def _collapse_internal_state(self):
        '''
        The internal state array may need to be collapsed sometimes.
        What this means: some entries may perfectly match at some stage 
        multiple library designs. If we are taking variable region number 1,
        and it happens to differ between two designs, and a particular entry
        matches both designs, which design do we go for? answer: we collapse 
        the internal state to a single match per design and then go for it
        '''
        if self._internal_state.ndim > 0:
            collapsed = np.zeros_like(self._internal_state, dtype=np.bool)
            ind = np.arange(self._internal_state.shape[0]), self._internal_state.argmax(axis=1)
            collapsed[ind] = self._internal_state[ind]
            self._internal_state = collapsed
            self._is_collapsed = True
        return
 
    def transform_all(self):
        self.transform('pep')
        self.transform('dna')
        self.transform('Q')
        
        #fix up the Q score arrays
        encoding = [''] + [chr(i) for i in range(128)]
        sort_idx = np.argsort(encoding)
        Q = np.searchsorted(encoding, self.Q, sorter = sort_idx) - 34
        Q[Q == -34] = 0
        self.Q = Q
        return
        
    def get_ndims(self):
        return self._validate_dims()

    def unpad_arrs(self):
        self.unpad('pep')
        self.unpad('dna')
        self.unpad('Q')
        return
     
class AnalysisSample(Sample):
    '''
    A generic machine learning/data analysis sample type: a sample is 
    associated with a name, data array (X) and (optionally) labels array (y).
    '''
    
    def __init__(self, X=None, y=None, name='unnamed', **kwargs):
        super(AnalysisSample, self).__init__(X=X, y=y, name=name, **kwargs)
        self._validate_name()
        self._validate_size()
        
        self.transform('X')
        return

    def get_ndims(self):
        return self.X.ndim

    def __repr__(self):
        return f'<AnalysisSample {self.name} containing {len(self)} entries>'  

class Data:
    '''
    #TODO: rewrite documentation for Data to include other sample types
    
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
        if not isinstance(self.samples, (type(None), list, tuple, set)):
            msg = f'Data object expected to receive a list of samples; received {type(self.samples)} instead'
            raise ValueError(msg)

        #TODO: make sure they are all the same class
        if self.samples:
            for sample in self.samples:
                if not issubclass(sample.__class__, Sample):
                    msg = f'Unexpected sample types passed to the Data object; type: {type(sample)}'
                    raise ValueError(msg)
            
        return

    def __iter__(self):
        for sample in self.samples:
            yield sample

    def __len__(self):
        if self.samples:
            return len(self.samples)
        else:
            return 0

    @property
    def size(self):
        return self.__len__()

    @property
    def dtype(self):
        if not self.samples:
            return type(None)
        else:
            assert len(set([type(x) for x in self.samples])), 'All data samples \
should have the same sample type!'
            return type(self.samples[0])
                
    def __getitem__(self, item):
        return self.samples[item]

    def __repr__(self):
        s = f'<Data container holding {len(self)} samples>:\n'
        if self.size < 10:
            for sample in self.samples:
                s += f'{sample}\n'
        else:
            for sample in self.samples[:5]:
                s += f'{sample}\n'     
            s += '...\n'            
            for sample in self.samples[-5:]:
                s += f'{sample}\n'           

        return s

    #TODO: this class will really benefit from implementing methods
    #to wrap around ops over individual samples







