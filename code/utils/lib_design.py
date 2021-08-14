# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 20:59:03 2020
@author: Vinogradov Alex
"""

import numpy as np    
class Template:

    def __init__(self, lib_seq='', monomers={}, wt=None, lib_type=None):
        '''
        The constructor that holds the information about the library design,
        i.e. the parent peptide sequence, amino acid positions subject to 
        randomization, constant regions, etc. Can store information about
        the library design at the DNA OR the peptide levels.
    
        At this stage, the library design constructor only works for libraries
        of fixed length (i.e. every peptide is the same length). 
        
        Parameters
        ----------
        lib_seq :  generic library sequence as dtype=str.
                   Positions subject to randomization should be
                   indicated as numerals.
                   e.g. the string 'LPENGA1111111111111111YPYDVPDYAGELARP'
                   indicates a peptide sequence that where 'LPENGA' is fixed 
                   for all library members, whereas the following 16 amino acids
                   are subject to randomization (mutagenesis, etc).
                                          
        monomers:  dict holding all possible monomer encodings for numerals
                   representing variable regions.
                   e.g.: {1: ('A', 'C', 'D', 'E' etc),
                          2: ('S', 'T'),
                          3: ('W', 'Y')}
        
        lib_type:  either 'dna' or 'pep' to specify what kind of library this is.
        
        Returns
        -------
        A class instance.

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
    A wrapper around all possible template sequences.
    Basically, just holds a list of Template class instances
    and a few methods to iterate over them. Monomers should necessarily be
    shared between individual templates. Library type should be shared as well.
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
        
       