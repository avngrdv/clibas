# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 16:44:27 2022
@author: Alex Vinogradov
"""

import os, gzip, re, inspect
from clibas.datatypes import Data, SequencingSample
from clibas.baseclasses import Handler
import numpy as np
import pandas as pd

class FastqParser(Handler):
    '''
    A processor for fastq/fastq.gz data. Primary parser for the sequencing
    data. The class holds methods for applying sequential filters to DNA
    sequencing data to eliminate noise, etc, and to convert raw NGS output
    to a list of peptides for the downstream applications. 
    
    Most public ops act on Data objects (except IO data fetchers) to
    return a transformed instance of Data.

    The class also holds a number of ops for basic statistics gathering.
    These also take Data as input, describe it in some way, write an out
    file (.png or txt or both) and return Data as-is.
    '''
    
    def __init__(self, *args):
        super(FastqParser, self).__init__(*args)
    
        self._validate_designs()
        self._validate_constants()
        return
        
    def __repr__(self):
        return '<FastqParser object>'
    
    def _dna_to_pep(self, seq, force_at_frame=None, stop_readthrough=False):       
             
        def find_orf(seq):
            loc = re.search(self.utr5_seq, seq)
            if loc is not None:
                return seq[loc.end()-3:]
            else:
                return None
        
        def find_stop(peptide):
            if stop_readthrough:
                return peptide
            
            else:
                ind = peptide.find('*')
                if ind == -1:
                    return peptide + '+'
                else:
                    return peptide[:ind]

        #figure out what to use as orf
        if force_at_frame is None:                
            orf = find_orf(seq)
        else:
            orf = seq[force_at_frame:]
            
        #throughout, '+' is a reserved symbol to denote messed up sequences
        #no stop codon, weird codons, etc
        pep = ''
        if orf is not None: 
            for i in range(0, len(orf), 3):
                try:
                    pep += self.constants.translation_table[orf[i:i+3]]
                except:
                    if len(orf[i:i+3]) != 3:
                        pep += '_'
                    else:
                        pep += '+'

        return find_stop(pep)

    #--------------------------------------------
    #The methods below are public data transformers.
    #All of them modify data in some way.
    #--------------------------------------------
    def translate(self, force_at_frame=None, stop_readthrough=False):
        '''
    	For each sample in Data, perform in silico translation for DNA sequencing data. 
    	The op will return data containing translated peptide lists. The op is 
        intended for one-ORF-per-read NGS data, but not for long, multiple-ORFs-per-read
        samples.
             
        This op should be called after fetching the data and (optionally) running
        the FastqParser.revcom(), prior to any filtration ops.
        
        On top of running translation, this op will also transform the data 
        to a reprensentation suitable for downstream ops.
        
        Parameters:
                force_at_frame: if None, a regular ORF search will be performed. Regular ORF
                                search entails looking for a Shine-Dalgarno sequence upstream 
                                of an ATG codon (the exact 5’-UTR sequence signalling an 
                                ORF is specified in config.py).
                                								
                                if not None, can take values of 0, 1 or 2. This will force-start
                                the translation at the specified frame regardless of the 
                                presence or absence of the SD sequence.
                                
                                For example:
                                DNA: TACGACTCACTATAGGGTTAACTTTAAGAAGGA
                   force_at_frame=0  ----------> 
                    force_at_frame=1  ---------->
                     force_at_frame=2  ---------->
                                 
              stop_readthrough:	bool (True/False; default: False). if True, translation will
                                continue even after encountering a stop codon until the 3'-end
                                of the corresponding read. Note, that an "_" amino acid will
                                be appended to the peptide sequence at the C-terminus if the 
                                last encountred codon is missing 1 or 2 bases.
                                
                                if False, the op will return true ORF sequences. In this case,
                                peptide sequences coming from ORFs which miss a stop codon will
                                be labelled with a "+" amino acid at the C-terminus.
                                
                                Should be flagged True for ORFs with no stop codon inside the read.
				 
        Returns:
                Data object containing peptide sequence information
        '''
        
        if force_at_frame is not None:
            if not isinstance(force_at_frame, int):
                msg = f'<translate> op expected to receive param "force_at_frame" as dtype=int; received: {type(force_at_frame)}'
                self.logger.error(msg)
                raise ValueError(msg)
        else:     
            if not hasattr(self, 'utr5_seq'):
                msg = "5' UTR sequence is not set for the <translation> op. Can not perform ORF search. Aborting. . ."
                self.logger.error(msg)
                raise ValueError(msg)   
                
        if type(stop_readthrough) != bool:
            msg = f'<translate> op expected to receive param "stop_readthrough" as type=bool; received: {type(stop_readthrough)}'
            self.logger.error(msg)
            raise ValueError(msg)   
                                
        def translate_dna(data):
            for sample in data:
                
                sample.pep = np.array([self._dna_to_pep(
                                                      x, 
                                                      force_at_frame=force_at_frame,
                                                      stop_readthrough=stop_readthrough
                                                     ) 
                                     
                                       for x in sample.dna])
                
                #this transformation is not declared publicly; may be it should be
                sample.transform_all()
        
                #set the internal state for the first time 
                shape = (sample.size, len(self.P_design))
                sample._internal_state = np.ones(shape, dtype=np.bool)
                
            return data
        return translate_dna
    
    def transform(self):
        '''
        Deprecated in favor of using FastqParser.translate(). If used, should
        be called after fetching the data and (optionally) running the
        FastqParser.revcom() op. Transforms the data to a representation 
        suitable for downstream ops.
        
        Parameters:
                None
    
        Returns:
                Transformed Data object
        '''
        def transform_data(data):
            for i, sample in enumerate(data):
                
                sample.transform_all()
                
                #<transform> sets the internal state for the first time 
                shape = (sample.size, len(self.P_design))
                sample._internal_state = np.ones(shape, dtype=np.bool)
            return data
        return transform_data
    
    def revcom(self):
        '''
        For each sample in Data, get reverse complement of DNA sequences and 
        reverse sequences of the corresponding Q score. If used, should enqueued 
        right after the fetching op, and before any downstream ops.
        
        Parameters:
                None
    
        Returns:
                Transformed Data object holding reverse-complemented DNA
                and reversed Q score information
        '''        

        if not hasattr(self.constants, 'complement_table'):
            msg = "Complement table is not specified for <revcom> op. Aborting. . ."
            self.logger.error(msg)
            raise ValueError(msg)   

        @np.vectorize
        def _rc(seq):
            return seq.translate(self.constants.complement_table)[::-1]

        @np.vectorize
        def _r(seq):
            return seq[::-1]
        
        def revcom_data(data):
            for sample in data:
    
                if sample.get_ndims() != 1:
                    msg = f'<revcom> can only be called on samples holding 1D-represented DNA. Ignoring the op for {sample.name} sample. . .'
                    self.logger.error(msg)
                    raise ValueError(msg)
                
                if sample.pep:
                    msg = 'Attempting to to revcom a sample holding a pep dataset. Pep dataset will be ignored. . .'
                    self.logger.warning(msg)
                    
                sample.dna = _rc(sample.dna)
                sample.Q = _r(sample.Q)
    
            return data
        return revcom_data

    def len_filter(self, where=None, len_range=None):
        '''
        For each sample in Data, filter out sequences longer/shorter than the specified 
        library designs. Alternatively, a length range of sequences to take can be optionally 
        specified to filter out the entries (NGS reads) outside of this range.
        
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
               len_range: either None (filtration will be done according to
                          the library design rules), or a list of two ints 
                          that specifies the length range to fetch.						  
					 
        Returns:
                Transformed Data object containing length-filtered data
        '''        
        self._where_check(where)
        design = self._infer_design(where)
        
        if len_range is not None:
            if not isinstance(len_range, (list, tuple)):
                msg = f'<len_filter> op expected to receive len_range argument as a list; received: {type(len_range)}'
                self.logger.error(msg)
                raise ValueError(msg)
            
            if len(len_range) != 2:
                msg = f'<len_filter> op expected to receive len_range as a list with two values; received: len={len(len_range)}'
                self.logger.error(msg)
                raise ValueError(msg)                

        def length_filter(data):
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]  
   
                #L is a length summary array
                L = self._L_summary(arr)
                
                #change the sample internal state
                for i,template in enumerate(design):
                    row_mask = sample._internal_state[:,i]
                    
                    if len_range is None:
                        sample._internal_state[row_mask, i] = L[row_mask] == template.L
                    else:
                        sample._internal_state[row_mask, i] = (L[row_mask] > len_range[0]) & (L[row_mask] < len_range[1])
                    
                #keep every entry that has at least one positive
                #value in the internal state array
                ind = np.any(sample._internal_state, axis=-1)
                sample.ind_filter(ind)
            
            return data
        return length_filter

    def cr_filter(self, where=None, loc=None, tol=1):
        '''
        For each sample in Data, filter out sequences not containing intact constant
        regions. Entries (NGS reads) bearing constant regions with amino acids outside
    	of the library design specification will be discarded.    
	
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
                     loc: a list of ints to specify which constant regions 
                          the op should process. 

                     tol: int; specifies the maximum allowed number of mutations
                          constant region fetched with where/loc before the 
                          entry (NGS read) is discarded. For the library from above
                          
                seq:      ACDEF11133211AWVFRTQ12345YTPPK
             region:      [-0-][---1--][--2--][-3-][-4-]
        is_variable:      False  True   False True False
                          
                          calling cr_filter(where='pep', loc=[2], tol=1), will
                          discard all sequences containing more than 1 mutation
                          in the 'AWVFRTQ' region. Note that the insertions/deletions
                          in the constant region are not validated by the parser.					  
					 
        Returns:
                Transformed Data object containg entries with intact 
                constant regions
        '''        
        self._where_check(where)
        design = self._infer_design(where)           
        self._loc_check(loc, design)            
        
        if not isinstance(tol, int):
            msg = f'<constant_region_filter> expected to receive parameter tol as as int; received: {type(tol)}'
            self.logger.error(msg)  
            raise ValueError(msg)

        if np.any(design.is_vr[loc]):
            msg = '<constant_region_filter> expected a list of contant regions to operate on; some of the specified locations point to variable regions.'
            self.logger.error(msg)
            raise AssertionError(msg)                
            
        def constant_region_filter(data):        
            from clibas.misc import hamming_distance
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]
                
                #iterativelt fill in the indexing array
                for i, template in enumerate(design):
                    
                    cr = np.array(template(loc))
                    cr_mask = template(loc, return_mask=True)
                    
                    row_mask = sample._internal_state[:,i]
                    if np.sum(row_mask) > 0:
                        dist = hamming_distance(arr[row_mask][:, cr_mask], cr, return_distance=True)
                        sample._internal_state[row_mask, i] = dist <= tol
                    else:
                        continue                    

                #keep every entry that has at least one positive
                #value in the internal state array
                ind = np.any(sample._internal_state, axis=-1)
                sample.ind_filter(ind)
                    
            return data
        return constant_region_filter

    def vr_filter(self, where=None, loc=None, sets=None):
        '''
        For each sample in Data, filter out sequences not containing intact variable 
        regions. Entries (NGS reads) bearing variable regions with amino acids outside
    	of the library design specification will be discarded.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
                     loc: a list of ints to specify which variable regions 
                          the op should process. 

                    sets: a list of ints; a list of monomer subsets to
                          check. For the library from above
                          
                seq:      ACDEF11133211AWVFRTQ12345YTPPK
             region:      [-0-][---1--][--2--][-3-][-4-]
        is_variable:      False  True   False True False
                          
                          there are five distinct variable amino acids:
                          1, 2, 3, 4, 5. The config file specifies which specific
                          amino acids are allowed for each of these numbers.
                          <vr_filter> op will make sure that each variable position
                          contains only the "allowed" monomers.					

                          vr_filter(where='pep', loc=[1], sets=[1, 3]) will make
                          sure that in region loc=1, variable amino acids 1 and 3
                          match the specification; variable amino acid 2 will not
                          be checked against in this example. Passing loc=[2] to
                          <vr_filter> op will raise an error, because it isn't a
                          variable region.
					 
        Returns:
                Transformed Data object containg entries with intact 
                variable regions
        '''
        self._where_check(where)
        design = self._infer_design(where)  
        self._loc_check(loc, design)
        
        if not isinstance(sets, list):
            msg = f'variable_region_filter op expected to receive a list of monomer subsets to parse; received: {type(sets)}'
            self.logger.error(msg)
            raise ValueError(msg)            

        allowed = set(design.monomers.keys())
        passed = set(sets)
        if not passed.issubset(allowed):
            msg = 'Specified variable region sets for <variable_region_filter> op must constitute a subset of library design monomers.'
            self.logger.error(msg)
            raise AssertionError(msg)

        if not np.all(design.is_vr[loc]):
            msg = '<variable_region_filter> expected a list of variable regions to operate on; some of the specified locations point to constant regions.'
            self.logger.error(msg)
            raise AssertionError(msg)
            
        def variable_region_filter(data):                     
            for sample in data:
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]

                #first things first: temporarily expand the internal
                #state array by one dimension; will collapse back at the end
                sample._internal_state = np.repeat(sample._internal_state[:,:,np.newaxis], len(sets), axis=-1)
                
                for i, template in enumerate(design): 
                    
                    #use internal state to figure out which
                    #entries are worth focusing on
                    row_mask = sample._internal_state[:,i,0]
                  
                    for j,subset in enumerate(sets):
                        
                        #work out column-wise mask
                        col_mask = np.array(template(loc, return_mask=True))
                        col_mask = col_mask[np.array(template(loc)) == subset]
                                      
                        #get the matching array: check whether entries are all in the corresponding monomer subset
                        match = np.in1d(arr[row_mask][:,col_mask], design.monomers[subset])
                        
                        #np.in1d flattens the array, so it needs to be reshaped back
                        match = match.reshape(arr[row_mask][:,col_mask].shape)
                                     
                        #the entry is taken only if everything matches
                        sample._internal_state[row_mask, i, j] = np.all(match, axis=1)
               
                #reduce along the subset axis to return
                #internal state array in its original form
                sample._internal_state = np.all(sample._internal_state, axis=-1)
                
                #keep every entry that has at least one positive
                #value in the internal state array
                ind = np.any(sample._internal_state, axis=-1)
                sample.ind_filter(ind)
                
            return data
        return variable_region_filter

    def filt_ambiguous(self, where=None):
        '''
        For each sample in Data, filter out sequences not containing intact ambiguous 
        tokens. For DNA, these are "N" nucleotides, which Illumina NGS ops occasionally
        assign during base calling. For peptides, these are any sequences containing
        amino acids outside of the translation table specification.	
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
        Returns:
                Transformed Data object containg entries without ambiguous
                tokens
        '''
        self._where_check(where)  
        allowed_monomers = self._infer_alphabet(where, alphabet=None)
            
        def filter_ambiguous(data):      
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]
                
                #perform the check; a little annoying because pads are also technically not allowed
                ind = np.in1d(arr, allowed_monomers).reshape(arr.shape)
                ind = np.sum(ind, axis=1) == self._L_summary(arr)
            
                #filter the sample
                sample.ind_filter(ind)
            
            return data
        return filter_ambiguous

    def drop_data(self, where=None):
        '''
        For each sample in Data, delete datasets specified in 'where'. See documentation 
        on Data objects above for more information.
    
        Parameters:
                   where: 'dna', 'pep' or 'q' to specify which datasets 
                          should be dropped. 				
						  
        Returns:
                Transformed Data object without dropped datasets
        '''
        if where not in ('pep', 'dna', 'Q'):
            msg = f"Invalid argument passed to <drop_dataset> op. Expected where = any of ('pep', 'dna', 'Q'); received: {where}"
            self.logger.error(msg)
            raise ValueError(msg)
        
        def drop_dataset(data):            
            for sample in data:
                sample.drop(where)
            return data
        return drop_dataset

    def q_score_filt(self, minQ=None, loc=None):
        '''
        For each sample in Data, filter out sequences associated with Q scores below 
        the specified threshold minQ.
    
        Parameters:
                     loc: a list of ints to specify which regions 
                          the op should process. 

                    minQ: every Q score in the regions specified 
                          by loc should be greater or equal than 
						  this value; everything else will be discarded
                        						  
						  
        Returns:
                Transformed Data object
        '''
        
        if not isinstance(minQ, int):
            msg = f'<Q_score_filter> op expected to receive parameter minQ as as int; received: {type(minQ)}'
            self.logger.error(msg)  
            raise ValueError(msg)

        self._where_check('dna')
        self._loc_check(loc, self.D_design)
        def q_score_filter(data):
            
            for sample in data:
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample.Q                          
                for i, template in enumerate(self.D_design):
                    
                    row_mask = sample._internal_state[:,i]
                    col_mask = template(loc, return_mask=True)

                    sample._internal_state[row_mask, i] = np.all(arr[row_mask][:,col_mask] >= minQ, axis=1)
        
                #keep every entry that has at least one positive
                #value in the internal state array
                ind = np.any(sample._internal_state, axis=-1)
                sample.ind_filter(ind)
            
            return data
        return q_score_filter

    def fetch_at(self, where=None, loc=None):
        '''
        For each sample in Data, for a dataset specified by 'where', fetch the regions
        specified by 'loc' and discard other sequence regions.
        
        Collapses sample's internal state.
        See documentation on Data objects for more information.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
                     loc: a list of ints to specify regions to be fetched 
						  
        Returns:
                Transformed Data object		
        '''
        self._where_check(where)               
        design = self._infer_design(where)
        self._loc_check(loc, design)        

        def _fetch_reg(arr, sample, design, loc):
            #initialize the array to hold the results
            max_len = self._find_max_len(design, loc)
            result = np.zeros((arr.shape[0], max_len), dtype=arr.dtype)
            
            for i, template in enumerate(design):
                
                col_mask = template(loc, return_mask=True)                        
                row_mask = sample._internal_state[:,i]
                
                result[row_mask, :len(col_mask)] = arr[row_mask][:,col_mask]
                #sample[where] = result
            return result
            
            
        def fetch_region(data):
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                
                if not sample._is_collapsed:
                    msg = f"<fetch_region> op will collapse sample {sample.name}'s internal state"
                    self.logger.info(msg)
                    sample._collapse_internal_state()
                    
                #has to be done manually, because in this case dna and pep args 
                #are not equivalent: where='dna' also means Q scores need to
                #truncated
                if where == 'pep':
                    arr = sample.pep
                    result = _fetch_reg(arr, sample, design, loc)
                    sample.pep = result
                    
                if where == 'dna':
                    arr = sample.dna
                    result = _fetch_reg(arr, sample, design, loc)
                    sample.dna = result
                    
                    arr = sample.Q
                    result = _fetch_reg(arr, sample, design, loc)
                    sample.Q = result                    
                
            #reindex the library design accordingly so that the downstream ops
            #can still be called with originally defined loc pointers
            design.truncate_and_reindex(loc)
                    
            return data
        return fetch_region
    
    def unpad(self):
        '''
        For each sample in Data, unpads the dna, pep and Q arrays. For each 
        array, removes the columns where every value is a padding token. 
        See documentation on Data objects for more information.

        Parameters:
                None	
						  
        Returns:
                Transformed Data object		
        '''        
        def unpad_data(data):
            for sample in data:
                sample.unpad_arrs()            
                                           
            return data            
        return unpad_data

    #--------------------------------------------
    #Below are IO readers/writers
    #--------------------------------------------    
    def _fetch_fastq_file(self, reader):
        '''
        Fetch DNA and Q score sequence lists from a .fastq file.
        .fastq files are base call .fastqs from single pair reads
        on Illumina's MiSeq instrument.
        
        in:            
            reader: a buffered reader with a preloaded file
        
        out:            
            DNA: a list of strings each containing a single read DNA sequence
            Q:   Q-scores corresponding to individual base calls, in the same format            
        '''        
        basename = os.path.basename(reader.name)
        sample_name = os.path.splitext(basename)[0]        
 
        with reader as f:
            msg = f'Fetching {basename}. . .'
            self.logger.info(msg)
            content = f.readlines()
            
            DNA = content[1::4]
            DNA = np.array([x.rstrip('\n') for x in DNA])
            
            Q = content[3::4]
            Q = np.array([x.rstrip('\n') for x in Q])
            
            f.close()
        
        sample = SequencingSample(
                                  name=sample_name,
                                  dna=DNA,
                                  Q=Q,
                                  pep=None
                                 )
        return sample
 
    def stream_from_fastq_dir(self, *args):
        '''
        A generator that yields data from self.fastq_dir sample by sample.
        Good when the entirety of the folder does not fit the memory.
        '''
        fnames = [os.path.join(self.dirs.seq_data, x) for x in os.listdir(self.dirs.seq_data) if x.endswith(".fastq")]
        if not fnames:
            msg = f'No .fastq files were found in {self.dirs.seq_data}! Aborting.'
            self.logger.error(msg)
            raise IOError(msg)
                
        for f in fnames:
            reader = open(f, 'r')
            sample = self._fetch_fastq_file(reader)
            yield sample
    
    def stream_from_gz_dir(self, *args):
        '''
        Fetch all .fastq.gz files from the sequencing_data directory 
        (as specified in config.py). Should be called as the first op in the workflow.
        
            Parameters:
                    None
        
            Returns:
                    Fetched Fastq data as an instance of Data
        '''
        fnames = [os.path.join(self.dirs.seq_data, x) for x in os.listdir(self.dirs.seq_data) if x.endswith(".gz")]
        if not fnames:
            msg = f'No .fastq.gz files were found in {self.dirs.seq_data}! Aborting.'
            self.logger.error(msg)
            raise IOError(msg)
                
        for f in fnames:
            reader = gzip.open(f, "rt")
            sample = self._fetch_fastq_file(reader)
            yield sample
                          
    def fetch_fastq_from_dir(self):
        '''
        Fetch all .fastq files from the sequencing_data directory 
        (as specified in config.py). Should be called as the first op in the workflow.
        
            Parameters:
                    None
        
            Returns:
                    Fetched Fastq data as an instance of Data
        '''
        def fetch_dir_fastq(*args):
            samples = list()
            for sample in self.stream_from_fastq_dir():
                samples.append(sample)
                
            return Data(samples=samples)
        return fetch_dir_fastq
       
    def fetch_gz_from_dir(self):
        '''
        Analogous to self.fetch_fastq_dir
        '''
        def fetch_dir_gz(*args):
            samples = list()
            for sample in self.stream_from_gz_dir():
                samples.append(sample)
                
            return Data(samples=samples)
        return fetch_dir_gz
    
    def save(self, where=None, fmt=None):
        '''
        For each sample in Data, save the dataset specified by 'where'. The results are written 
        to a file in the parser output folder as specified by config.py.
        
        Parameters:					  
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.

                     fmt: the format of the output file. Supported values are
                          'npy', 'fasta' and 'csv'					 
                          						  							  
        Returns:
                Data object (no transformation)
        '''
        
        if fmt not in ('npy', 'csv', 'fasta'):
            msg = f"<save_data> op received invalid fmt argument. Acceted any of ('npy', 'csv', 'fasta'); received: {fmt}"
            self.logger.error(msg)
            raise ValueError(msg)
        
        self._where_check(where)
        def _writer(arr, fmt, path):
            
            if fmt == 'npy':
                np.save(path + '.npy', arr)    
                return
            
            arr_1d = [''.join(x) for x in arr]
            
            if fmt == 'csv':           
                with open(path + '.csv', 'w') as f:
                    for seq in arr_1d:
                        f.write(f'{seq},\n')    
                return
                
            if fmt == 'fasta':
                with open(path + '.fasta', 'w') as f:
                    for i,seq in enumerate(arr_1d):
                        f.write(f'>sequence_{i}\n')  
                        f.write(f'{seq}\n')
                return                
                
        def save_data(data):
                    
            self._prepare_destinations(data, self.dirs.parser_out)
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]  

                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f'{sample.name}_{where}'
                path = os.path.join(destination, fname)
                
                _writer(arr, fmt, path)
                
            return data
        return save_data    

    def count_summary(self, where=None, top_n=None, fmt=None):
        '''
        For each sample in Data, counts the number of times each unique 
        sequence is found in the dataset specified by 'where'. The results 
        are written to a file in the analysis folder as specified by config
        
        Parameters:					  
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
    
                   top_n: if None, full summary will be created. If
                          an int is passed, only top_n sequences (by count)
                          will be written to a file.
    
                     fmt: the format of the output file. Supported values are
                          'csv' and 'fasta'.					 
                          						  							  
        Returns:
                Data object (no transformation)
        '''
        
        self._where_check(where)
        if fmt not in ('csv', 'fasta'):
            msg = f"<fastq_count_summary> op received invalid fmt argument. Acceted any of ('csv', 'fasta'); received: {fmt}"
            self.logger.error(msg)
            raise ValueError(msg)
            
        if top_n is not None:
            if not isinstance(top_n, (int, np.int)):
                msg = f'<fastq_count_summary> op expected to receive parameter top_n as as int; received: {type(top_n)}'
                self.logger.error(msg)  
                raise ValueError(msg)
    
        def _writer(sample, og_ind, counts, fmt, path):
            if fmt == 'csv':           
                df = pd.DataFrame(columns=['Peptide', f'{where} count', 'DNA'])
                df['Peptide'] = [''.join(x) for x in sample.pep[og_ind]]
                df['DNA'] = [''.join(x) for x in sample.dna[og_ind]]
                df[f'{where} count'] = counts
                df.to_csv(path + '.csv', sep=',')
                
            if fmt == 'fasta':
                arr = sample[where][og_ind]
                arr_1d = [''.join(x) for x in arr]    
                
                with open(path + '.fasta', 'w') as f:
                    for i,seq in enumerate(arr_1d):
                        f.write(f'>seq_{i+1}_count_{counts[i]}\n')  
                        f.write(f'{seq}\n')
                return 
        
        def fastq_count_summary(data):
            
            self._prepare_destinations(data, self.dirs.parser_out)
            for sample in data:
            
                from clibas.misc import sorted_count
                _, og_ind, counts = sorted_count(sample[where],
                                                 top_n=top_n,
                                                 return_index=True
                                                )
                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f'{sample.name}_{where}_count_summary'
                path = os.path.join(destination, fname)
                
                _writer(sample, og_ind, counts, fmt, path)
                                
            return data
        return fastq_count_summary
    
    def library_design_match(self, where=None):
        '''
        For each sample in Data, compute the number of matches between the dataset 
        specified by 'where' and the corresponding library templates. The results 
        are written to a file in the analysis folder as specified by config
        
        In other words, summarize where dataset sequences come from (from which
        libraries). The op could also be called "_internal_state_summary"
        
        Parameters:					  
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on                        						  							  
        Returns:
                Data object (no transformation)
        '''        
    
        self._where_check(where)
        def library_design_summary(data):
            
            design = self._infer_design(where)
            if not isinstance(data, Data):
                msg = f'<library_design_match_analysis> op expected data as Data type; received: {type(data)}'
                self.logger.error(msg)
                raise TypeError(msg)
                
            #summarize straight into a pandas dataframe
            sample_names = [sample.name for sample in data]
            templates = [template.lib_seq for template in design]
        
            import pandas as pd
            df = pd.DataFrame(index=sample_names, columns=templates)
            
            #all this op is: axis=0-wide sum of the internal states
            for sample in data:
                df.loc[sample.name] = np.sum(sample._internal_state, axis=0)
                
            fname = f'{self.exp_name}_by_template_breakdown.csv'
            path = os.path.join(self.dirs.parser_out, fname)            
            df.to_csv(path + '.csv', sep=',')
    
            return data
        return library_design_summary    