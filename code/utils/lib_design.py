# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 20:59:03 2020
@author: Vinogradov Alex
"""

import numpy as np    
class Template:

    def __init__(self, lib_seq='', monomers={}, wt=None, lib_type=None):
        '''
        The constructor that holds the information about the library for a single
        template, i.e. the generic peptide sequence, amino acid positions subject to 
        randomization, constant regions, etc. Can store information about
        the library design at the DNA OR the peptide levels. See the comments to 
        LibraryDesign class for more information. 
        '''
        self.lib_seq = lib_seq
        self.monomers = monomers
        self.lib_type = lib_type
        
        self._typecheck()
        self._build()
        return
    
    def _typecheck(self):
    
        if not isinstance(self.lib_seq, str):
            raise ValueError("Library design must be specified as a string (dtype=str). . .")

        #make sure that monomer classes match what's specified by the library sequence
        lib_seq_monomer_types = set(int(x) for x in self.lib_seq if x.isdigit())
        specified_monomer_types = set(self.monomers.keys())
        
        if not lib_seq_monomer_types.issubset(specified_monomer_types):
            raise ValueError("Variable region monomer names must match what's specified in monomers. . .")

        if self.lib_type != 'dna' and self.lib_type != 'pep':
            raise ValueError('Library type can only accept either "dna" or "pep" as valid values')
    
        import config
        if self.lib_type == 'dna':
            lookup_monomers = config.constants.bases
        else:
            lookup_monomers = config.constants.aas
        
        for m in self.lib_seq:
            if not m.isdigit():
                if not m in lookup_monomers:
                    raise ValueError('All library design monomers must be specified in the lookup tables. . .') 
                
        for key in self.monomers:
            for m in self.monomers[key]:
                 if not m in lookup_monomers:
                    raise ValueError('All library design monomers must be specified in the lookup tables. . .') 
        
        return
    
    def _build(self):

        #L is expected length of the library sequence
        #may be meaningless for libraries of variable length
        self.L = len(self.lib_seq)


        #build the key params: is_vr, region, mask
        is_vr = []
        region = []
        mask = []
        
        current_region = list()
        current_mask = list()
        
        for i,m in enumerate(self.lib_seq):
            if i == 0:
                is_vr.append(m.isdigit())
                
            if m.isdigit():     
                if is_vr[-1] is True:        
                    current_region.append(int(m))
                    current_mask.append(i)
                    
                else:
                    region.append(current_region)
                    mask.append(current_mask)
                    
                    current_region = [int(m)]
                    current_mask = [i]
                    is_vr.append(True)
                
            else:
                if not is_vr[-1] is True:        
                    current_region.append(m)
                    current_mask.append(i)
                    
                else:
                    region.append(current_region)
                    mask.append(current_mask)
                    
                    current_region = [m]
                    current_mask = [i]
                    is_vr.append(False)        
                
        region.append(current_region)
        mask.append(current_mask)
        
        self.is_vr = np.array(is_vr, dtype=np.bool)
        self.loc = np.arange(self.is_vr.size)
        self.mask = mask
        self.region = region
        return

    def __repr__(self):
        seq = ''.join(self.lib_seq)    
        return f'<Template container for {seq}  lib_type={self.lib_type}>'    

    def __call__(self, loc, return_mask=False):
        '''
        Run fancy indexing to get the sequence of regions at a specified location 'loc'
        '''
        def fancy_index(arr, loc):
            out = []
            for x in loc:
                out.extend(arr[x])
                
            return out
        
        if return_mask:
            arr = self.mask
        else:
            arr = self.region
        
        return fancy_index(arr, loc)

class LibraryDesign:
    '''
    LibraryDesign is a key object that specifies what kind of library the parser 
    should expect. The object can hold information about arbitrary DNA and 
    peptide libraries using a unified logic as follows. Randomized amino 
    acids/bases (hereafter _tokens_) are indicated as numerals (0-9), whereas
    tokens which are not subject to randomization (linker sequences, etc) are 
    indicated using the standard one letter encoding (A, C, T, G for DNA), 
    (A, C, D etc for peptides; the encoding must make sense according to the
     translation table in config.py). A continuous stretch of either random or
    fixed tokens makes up a “region” in the template sequence. For example:
    
                    seq:      ACDEF11133211AWVFRTQ12345YTPPK
                 region:      [-0-][---1--][--2--][-3-][-4-]
            is_variable:      False  True   False True False
    
    Region assignments are made automatically; in the example above, the 
    library contains 5 regions; 3 are “constant regions” and 2 are “variable 
    regions”. Numerals used for variable region tokens should be defined, with 
    one number corresponding to a particular token set. For example, an NNK
    codon encodes all 20 amino acids, whereas an NNC codon only 15. Thus, all
    amino acids derived from NNK codons should be encoded by one number, and
    another number for NNC-encoded positions. LibraryDesign can take several
    templates of different length to encode libraries with variable regions
    of variable size. Below is an example of a LibraryDesign initialization:
    
        lib = LibraryDesign(
            
                        templates=[
                                    '211113GSGSGS',
                                    '2111113GSGSGS',
                                    '21111113GSGSGS',
                                  ],
                
                        monomers={
                                  1: ('A', 'C', 'D', 'E', 'F', 'G', 'H'),
                                  2: ('M'),
                                  3: ('C')
                                 },
                        
                        lib_type='pep'
                            
                           )
    
    Note that variable positions can encode a single amino acid 
    (amino acids 2 and 3). In this way, there is a considerable flexibility in
    how a particular library can be represented. When initializing LibraryDesign 
    objects, several rules must be followed:
    
    1.	The topology of every passed template must be identical. Topology is 
    the total number of regions, and the total number of variable regions.
    Essentially, the templates should only differ in the internal composition 
    of variable regions.
    2.	All variable region monomers should be encoded in the translation table
    (or be one of the four standard DNA bases for DNA libraries; bases N, K etc
     should be converted to numerals).
    3.	Two LibraryDesign objects should be created for the parser 
    (lib_type=’dna’ and lib_type=’pep’). They should have the same number 
    of templates. The first DNA template should give rise to the first peptide
    template and so on. LibraryDesign objects are defined in the config file. 
    '''
    
    def __init__(self, templates=[], monomers={}, lib_type=None):
        '''
        templates should be passed as a list of library design strings
        '''
        self.monomers = monomers
        self.lib_type = lib_type
        self.templates = tuple(
                               Template(lib_seq=x, 
                                        monomers=self.monomers, 
                                        lib_type=self.lib_type) 
                               
                               for x in templates
                              )

        self.L = list(set([x.L for x in self.templates]))
        self._topology_check()
        
        self.loc = self.templates[0].loc
        self.is_vr = self.templates[0].is_vr
        return        
        
    def __repr__(self):
        return f'<Library design container for {len(self.templates)} templates. lib_type={self.lib_type}>'
        
    def __iter__(self):
        for template in self.templates:
            yield template
    
    def __len__(self):
        return len(self.templates)

    def __getitem__(self, item):
        return self.templates[item]

    def _topology_check(self):
        '''
        Make sure that all templates have one topology of
        variable and constant regions. (e.g., constant-variable-constants)
        In fact, it is enough to check that the number of crs and vrs 
        is identical for every template in the library.
        '''
        topologies = []
        for t in self.templates:
            topologies.append(tuple(t.is_vr))
        
        if len(set(topologies)) != 1:
            raise ValueError('All library templates should have the same topology. . .') 
        
        return
        
       