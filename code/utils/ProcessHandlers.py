# -*- coding: utf-8 -*-
"""
Created on Fri May 21 20:28:10 2021
@author: Alex Vinogradov
"""

import time, os, logging, gzip, re, inspect
from utils import Plotter
from utils.datatypes import Data, SequencingSample

import numpy as np
import pandas as pd

class Logger:
    '''
    A decorated version of the standard python logger object.
    Can be setup from a config file. Two main customizations
    are implemented: verbosity (whether logger messages should be
    printed to the running stream) and log_to_file, which
    if set, will setup a dedicated handler to dump log info
    to a file.
    '''
    
    def __init__(self, config=None):
        self.conf = config
    
        self.__fallback()
        self.__configure_logger()
        return

    def __repr__(self):
        return f'<Logger {self.name}; verbose: {self.verbose}; log_to_file: {self.log_to_file}>'

    def __fallback(self):
        '''
        If no config is passed fallback to some innocuous defaults,
        which is basically a silent logger.
        '''
        
        if self.conf is None:
        
            self.name = 'unnamed ' + str(time.time())
            self.verbose = False
            self.log_to_file = False
            self.log_fname = None
            
        else:

            self.name = self.conf.name
            self.verbose = self.conf.verbose
            self.log_to_file = self.conf.log_to_file
            self.log_fname = self.conf.log_fname
        
        return

    def __configure_logger(self):

        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter("[%(levelname)s]: %(message)s")

        #clear any preexisting handlers to avoid stream duplication
        self.logger.handlers.clear()

        if self.verbose:
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            console_handler.setLevel(logging.INFO)

            self.logger.addHandler(console_handler)

        if self.log_to_file:
               
            filehandler = logging.FileHandler(self.log_fname)      
            filehandler.setFormatter(formatter)
            filehandler.setLevel(logging.INFO)
                    
            self.logger.addHandler(filehandler)       
        
        return
        
class DirectoryTracker:
    '''
    A simple object to keep track where the data,  
    logs, etc should belooked for. Config should be
    passed to customize the directories. Otherwise,
    everything will be looked up and dumped into cwd.
    '''
    
    def __init__(self, config=None):

        self._conf = config
        self.__fallback()
        self.__setup_dirs()

    def __repr__(self):
        return '<DirectoryTracker object>'

    def __fallback(self):
        '''
        If no config is passed, fallback to some preset defaults.
        At present, all directories will be set to cwd if config
        is not specified.
        '''
        cwd = os.getcwd()
        if self._conf is None:
            
            self.seq_data = cwd
            self.logs = cwd
            self.parser_out = cwd
            
        else:
            self.seq_data = self._conf.seq_data
            self.logs = self._conf.logs
            self.parser_out = self._conf.parser_out
            
        return
            
    def __setup_dirs(self):
                   
        for d in [x for x in dir(self) if not x.startswith('_')]:
            if not os.path.isdir(getattr(self, d)):
                os.makedirs(getattr(self, d))
        return        
    
class Handler:
    '''
    Base handler class. Should not be invoked directly.
    '''
    
    def __init__(self, *args):
        self.__dict__.update(*args)
        self.__logger_fallback()
        self.__tracker_fallback()
        return

    def __logger_fallback(self):
        '''
        If no logger is passed to a data handler, a default Logger 
        object will be invoked. The default logger is silent. 
        '''
        
        if not hasattr(self, 'logger'):
            self.logger = Logger().logger
    
        if self.logger is None:
            self.logger = Logger().logger
            
        return
    
    def __tracker_fallback(self):
        '''
        If no DirectorTracker object was passed to a handler, a default
        tracker will be invoked (everything in the cwd)
        '''        
        
        if not hasattr(self, 'dirs'):
            self.dirs = DirectoryTracker()
            pass
        
        return
        
    def _on_completion(self):
        msg = f'The following handler was succesfully initialized: {self}'
        self.logger.info(msg)
        return

class Pipeline(Handler):
    
    def __init__(self, *args):
        super(Pipeline, self).__init__(*args)
        self._on_startup()
        
        super(Pipeline, self)._on_completion()
        return

    def __repr__(self):
        return f'<Pipeline object; current queue size: {len(self.que)} routine(s)>'

    def _on_startup(self):
        self.que = []
        if not hasattr(self, 'exp_name'):
            self.exp_name = 'unnamed'
        return

    def _describe_data(self, data=None):
        '''
        Go over every dataset for every sample and
        log all array shapes. Used during dequeing
        to keep track of data flows.
        '''
        data_descr = []
        
        if data is None:
            return data_descr
        
        for sample in data:
            
            data_descr.append((sample.name, len(sample)))
            
            for tup in sample:
                
                if tup[0].shape:
                    shape = tup[0].shape
                else:
                    shape = None
                
                msg = f'{sample.name} {tup[1]} dataset shape: {shape}'
                self.logger.info(msg)
                
        msg = 65 * '-'
        self.logger.info(msg)
    
        return data_descr

    def _reassemble_summary(self, summary):
        
        ops = []
        times = []
        samples = []
        
        #code below is a mess, but the task is trivial, 
        #so whatever; fix if nothing better to do
        for x in summary:
            ops.append(x['op'])
            times.append(x['op_time'])
            for j in x['data_description']:
                samples.append(j[0])
        
        samples = list(set(samples))
        sizes = np.zeros((len(summary), len(samples)))
        for i,entry in enumerate(summary):
            for tup in entry['data_description']:
                for j, name in enumerate(samples):
                    if tup[0] == name:
                        sizes[i,j] = tup[1]
        
        df = pd.DataFrame(columns=['time'] + samples, index=ops)
        df['time'] = times
        for i,name in enumerate(samples):
            df[name] = sizes[:,i]
        
        return df
    
    def enque(self, routines):
        '''
        Takes a list of functions and adds them to the pipeline queue.
        self.deque will take some data as an argument and apply dump
        the queue on it, i.e. sequentially transform the data by applying
        the queued up routines. 

        Parameters
        ----------
        routines : a list of functions capable of acting on data.
                   every routine should take data as the only argument
                   and return transformed data in the same format (Data object)

        Returns
        -------
        None.

        '''
        
        for func in routines:
            self.que.append(func)
            
        msg = f'{len(routines)} routines appended to pipeline; current queue size: {len(self.que)}'
        self.logger.info(msg)
        
        return        

    def run(self, data=None, save_summary=True):
        '''
        Chainlinks the list of routines one by one to 
        sequentially transform the data. The method will
        basically execute the queued up experiment.

        Parameters
        ----------
        data : Data object or None
               if None, the first func in the que
               has to load the data

        save_summary: save a .csv summary file containing
                      the progress of the experiment and
                      the basic description of data at
                      every stage. location: logs

        Returns
        -------
        transformed data as a Data object

        '''
        summary = list()
        data_descr = self._describe_data(data)
        summary.append({'op': None, 'op_time': None, 'data_description': data_descr})
        
        for _ in range(len(self.que)):
        
            func = self.que.pop(0)
            msg = f'Queuing <{func.__name__}> routine. . .'
            self.logger.info(msg)
            
            t = time.time()
            data = func(data)
            op_time = np.round(time.time() - t, decimals=3)
            
            msg = f'The operation took {op_time} s'
            self.logger.info(msg)
            data_descr = self._describe_data(data)
            
            summary.append({'op': func.__name__, 'op_time': op_time, 'data_description': data_descr})
        
        if save_summary:
            summary = self._reassemble_summary(summary)
            fname = f'{self.exp_name}_pipeline_summary.csv'
            path = os.path.join(self.dirs.logs, fname)
            summary.to_csv(path)
        
        return data   

class FastqParser(Handler):
    '''
    A processor for fastq/fastq.gz data. Primary parser for the sequencing
    data. The class holds methods for applying sequential filters to DNA
    sequencing data to eliminate noise, etc, and to convert raw NGS output
    to a list of peptides for the downstream applications. 
    
    Most public routines act on Data objects (except IO data fetchers) to
    return a transformed instance of Data.

    The class also holds a number of ops for basic statistics gathering.
    These also take Data as input, describe it in some way, write an out
    file (.png or txt or both) and return Data as-is.
    '''
    
    def __init__(self, *args):
        super(FastqParser, self).__init__(*args)
    
        self._validate()
        super(FastqParser, self)._on_completion()
        return    
        
    def __repr__(self):
        return '<FastqParser object>'

    def _validate(self):
        
        if not (hasattr(self, 'P_design') and
                hasattr(self, 'D_design')
               ):
            msg = 'FastqParser requires peptide and DNA library design objects for setup. . .'
            self.logger.error(msg)
            raise ValueError(msg)
            
        if not len(self.P_design) == len(self.D_design):
            msg = 'Peptide and DNA library designs must contains the same number of templates; cannot inialize FastqProcessor. . .'
            self.logger.error(msg)
            raise ValueError(msg)
            
        if not hasattr(self, 'constants'):
            msg = 'FastqParser requires constants for setup. . .'
            self.logger.error(msg)
            raise ValueError(msg)
        return

    def _transform_check(self, sample, func):
        if not sample.get_ndims() == 2:
            raise ValueError(f'Sample {sample.name} holds arrays of unsupported dimensionality for {func} op. Expected: arrays of ndims=2, got: ndims={sample.get_ndims()}')
        
        return
    
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
                    pep += self.constants.codon_table[orf[i:i+3]]
                except:
                    if len(orf[i:i+3]) != 3:
                        pep += '_'
                    else:
                        pep += '+'

        return find_stop(pep)

    def _L_summary(self, arr):

        #infer what the pad token is
        pad = np.zeros(1, dtype=arr.dtype)[0]
        
        #fetch the indexes where dna/pep length == designed
        return np.sum(arr != pad, axis=1)

    def _find_max_len(self, design, loc):
        '''
        When trying to get a column-wise view of the array,
        The views for different designs can have a different
        shape (for example, different vr size). This will 
        find the largest possible column-wise view.
        Output m is used to as a shape parameter during
        array creation.
        '''
        m = 0 
        for template in design:
            if len(template(loc)) > m:     
                m = len(template(loc))
        return m

    def _where_check(self, where):
        
        if where == 'pep':
            if not hasattr(self, 'P_design'):
                msg = "Cannot run peptide filtration routines with unspecified library design."
                self.logger.error(msg)
                raise ValueError(msg)
            
        elif where == 'dna':
            if not hasattr(self, 'D_design'):
                msg = "Cannot run dna filtration routines with unspecified library design."
                self.logger.error(msg)
                raise ValueError(msg)                       
         
        else:
            msg = f'The parser did not understand which dataset it should operate on. Passed value: {where}; allowed values: pep/dna.'
            self.logger.error(msg)
            raise ValueError(msg)
        return

    def _loc_check(self, loc, design):
        
        if not isinstance(loc, list):
            msg = f'The Parser expected to receive a list of region indexes to parse; received: {type(loc)}'
            self.logger.error(msg)
            raise ValueError(msg)

        if max(loc) > design.loc.max():
            msg = f'{design.lib_type} library design does not contain enough regions. Library design contains {design.loc.max() + 1} regions; specified: up to {max(loc) + 1}'
            self.logger.error(msg)
            raise AssertionError(msg)
        return

    def _prepare_destinations(self, data):
        
        for sample in data:
                destination = os.path.join(self.dirs.parser_out, sample.name)
                if not os.path.isdir(destination):
                    os.makedirs(destination)        
        return    
    
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
        the FastqParser.revcom(), prior to any filtration routines.
        
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
            if force_at_frame not in (0, 1, 2):
                msg = f'<translate> routine expected to receive param "force_at_frame" as any of (0, 1, 2); received: {force_at_frame}'
                self.logger.error(msg)
                raise ValueError(msg)
        else:     
            if not hasattr(self, 'utr5_seq'):
                msg = "5' UTR sequence is not set for the <translation> routine. Can not perform ORF search. Aborting. . ."
                self.logger.error(msg)
                raise ValueError(msg)   
                
        if type(stop_readthrough) != bool:
            msg = f'<translate> routine expected to receive param "stop_readthrough" as type=bool; received: {type(stop_readthrough)}'
            self.logger.error(msg)
            raise ValueError(msg)   
                                
        def translate_dna(data):
            for sample in data:
                
                sample.P = np.array([self._dna_to_pep(
                                                      x, 
                                                      force_at_frame=force_at_frame,
                                                      stop_readthrough=stop_readthrough
                                                     ) 
                                     
                                     for x in sample.D])
                
                #this transformation is not declared publicly; may be it should
                sample.transform()
        
                #set the internal state for the first time 
                shape = (len(sample), len(self.P_design))
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
                
                sample.transform()
                
                #<transform> sets the internal state for the first time 
                shape = (len(sample), len(self.P_design))
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

        @np.vectorize
        def _rc(seq):
            return seq.translate(self.constants.complement_table)[::-1]

        @np.vectorize
        def _r(seq):
            return seq[::-1]
        
        def revcom_data(data):
            for sample in data:
    
                if sample.D.ndim != 1 or sample.Q.ndim != 1:
                    msg = f'<revcom> can only be called on samples holding 1D-represented DNA. Ignoring the routine for {sample.name} sample. . .'
                    self.logger.warning(msg)
                    continue
                
                if sample.P:
                    msg = 'Attempting to to revcom a sample holding a P dataset. P dataset will be ignored. . .'
                    self.logger.warning(msg)
                    
                sample.D = _rc(sample.D)
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

        if where == 'pep':
            design = self.P_design
            
        elif where == 'dna':
            design = self.D_design
        
        if len_range is not None:
            if not isinstance(len_range, list):
                msg = f'<len_filter> routine expected to receive len_range argument as a list; received: {type(len_range)}'
                self.logger.error(msg)
                raise ValueError(msg)
            
            if len(len_range) != 2:
                msg = f'<len_filter> routine expected to receive len_range as a list with two values; received: len={len(len_range)}'
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
                sample(ind)
            
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

        if where == 'pep':
            design = self.P_design
        
        elif where == 'dna':
            design = self.D_design
            
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
            from utils.misc import hamming_distance
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
                sample(ind)
                    
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
        if where == 'pep':
            design = self.P_design
        
        elif where == 'dna':
            design = self.D_design
        
        self._loc_check(loc, design)
        if not isinstance(sets, list):
            msg = f'variable_region_filter routine expected to receive a list of monomer subsets to parse; received: {type(sets)}'
            self.logger.error(msg)
            raise ValueError(msg)            

        allowed = set(design.monomers.keys())
        passed = set(sets)
        if not passed.issubset(allowed):
            msg = 'Specified variable region sets for <variable_region_filter> routine must constitute a subset of library design monomers.'
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
                sample(ind)
                
            return data
        return variable_region_filter

    def filt_ambiguous(self, where=None):
        '''
        For each sample in Data, filter out sequences not containing intact ambiguous 
        tokens. For DNA, these are "N" nucleotides, which Illumina NGS routines occasionally
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
        
        #fetch the relevant monomer sets
        if where == 'pep':
            allowed_monomers = self.constants.aas
            
        elif where == 'dna':
            allowed_monomers = self.constants.bases
            
        def filter_ambiguous(data):      
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]
                
                #perform the check; a little annoying because pads are also technically not allowed
                ind = np.in1d(arr, allowed_monomers).reshape(arr.shape)
                ind = np.sum(ind, axis=1) == self._L_summary(arr)
            
                #filter the sample
                sample(ind)
            
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
        if where not in ('pep', 'dna', 'q'):
            msg = f"Invalid argument passed to <drop_dataset> routine. Expected where = any of ('pep', 'dna', 'q'); got: {where}"
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
            msg = f'<Q_score_filter> routine expected to receive parameter minQ as as int; received: {type(minQ)}'
            self.logger.error(msg)  
            raise ValueError(msg)

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
                sample(ind)
            
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
        if where == 'pep':
            design = self.P_design
        
        elif where == 'dna':
            design = self.D_design
        
        self._loc_check(loc, design)        
        
        def fetch_region(data):
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]
                
                if not sample._is_collapsed:
                    msg = f"<fetch_region> routine will collapse sample {sample.name}'s internal state"
                    self.logger.info(msg)
                    sample._collapse_internal_state()
                
                #initialize the array to hold the results
                max_len = self._find_max_len(design, loc)
                result = np.zeros((arr.shape[0], max_len), dtype=arr.dtype)
                
                for i, template in enumerate(design):
                    
                    col_mask = template(loc, return_mask=True)                        
                    row_mask = sample._internal_state[:,i]
                    
                    result[row_mask, :len(col_mask)] = arr[row_mask][:,col_mask]
                    sample[where] = result
                    
            #reindex the library design accordingly so that the downstream ops
            #can still be called with originally defined loc pointers
            design.truncate_and_reindex(loc)
                    
            return data
        return fetch_region
    
    def unpad(self):
        '''
        For each sample in Data, unpads the D, Q, P arrays. For each array, removes 
    	the columns where every value is a padding token. See documentation on Data 
        objects for more information.

        Parameters:
                None	
						  
        Returns:
                Transformed Data object		
        '''        
        def unpad_data(data):
            for sample in data:
                sample.unpad()            
                                           
            return data            
        return unpad_data

    #--------------------------------------------
    #The methods below do not transform the data.
    #They are only used to assemble statistics, 
    #plots the results, etc.
    #--------------------------------------------
    def len_summary(self, where=None, save_txt=False):
        '''
        For each sample in Data, compute the distribution of peptide/DNA sequence lengths
        (specified by 'where') and plot the resulting histogram in the parser output folder
        as specified by config.py. Optionally, the data can also be written to a txt file.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
                          						  
                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots				
						  
        Returns:
                Data object (no transformation)
        '''
        self._where_check(where)
        def length_summary(data):
            
            self._prepare_destinations(data)        
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]
                
                L = self._L_summary(arr)            
                L, counts = np.unique(L, return_counts=True)
                
                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f'{sample.name}_{where}_L_distribution'
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.L_distribution(L, counts, where, basename)
                
                if save_txt:
                    np.savetxt(basename + '.csv',
                               np.array((L, counts)).T,
                               delimiter=',', 
                               header='Seq length,Count')
            
            return data
        return length_summary

    def convergence_summary(self, where=None):
        '''
        For each sample in Data, perform basic library convergence analysis on a sequence 
        level. Computes normalized Shannon entropy, and postition-wise sequence conservation. 
        Plots the results in the parser output folder as specified by config.py.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
                          						  							  
        Returns:
                Data object (no transformation)
        '''
        self._where_check(where)
        
        if where == 'pep':
            tokens = self.constants.aas
                    
        elif where == 'dna':               
            tokens = self.constants.bases
        
        from utils.misc import shannon_entropy, get_freqs
        def _seq_conservation(freq):
            '''
            NOTE: this computation doesn't really make sense for
            arrays containing sequences of uneven length. Rather,
            the meaning becomes somewhat counterintuitive, but
            what to do about it?
            '''
            with np.errstate(divide='ignore', invalid='ignore'):
                em = np.nan_to_num(np.multiply(freq, np.log2(freq)))

            return np.sum(em, axis=0) + np.log2(freq.shape[0])    
            
        def library_convergence_summary(data):

            self._prepare_destinations(data)
            for sample in data:

                self._transform_check(sample, inspect.stack()[0][3])                
                arr = sample[where]
                
                shannon, counts = shannon_entropy(arr, norm=True)
                freq = get_freqs(arr, tokens)
                seq_conservation = _seq_conservation(freq)
                
                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f'{sample.name}_{where}_library_convergence'
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.dataset_convergence(counts, shannon, where, basename)                
                
                fname = f'{sample.name}_{where}_sequence_conservation'
                basename = os.path.join(destination, fname)                
                Plotter.SequencingData.conservation(seq_conservation, where, basename)
                
            return data
        return library_convergence_summary

    def freq_summary(self, where=None, loc=None, save_txt=False):
        '''
        Perform basic library convergence analysis at a token level. For each sample in Data, 
        computes the frequency of each token in the dataset. Plots the results in the parser 
        output folder as specified by config.py. Optionally, the data can also be written to
        a txt file.

        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
						  
                     loc: a list of ints to specify regions to be analyzed;
                          in this case, the op will collapse sample's internal
                          state (see explanation for Data objects)
                          
                          OR
                          
                          'all': to get the same statistics over the entire sequence;
                                 in this case, the op will NOT collapse sample's 
                                 internal state

                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots						 
                          						  							  
        Returns:
                Data object (no transformation)
        '''
        self._where_check(where)
        if where == 'pep':
            design = self.P_design
            tokens = self.constants.aas
            
        elif where == 'dna':
            design = self.D_design        
            tokens = self.constants.bases
            
        if loc != 'all':
            self._loc_check(loc, design)
            
        from utils.misc import get_freqs
        
        def frequency_summary(data):
            self._prepare_destinations(data)
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]
                
                if loc == 'all':
                    freq = get_freqs(arr, tokens)
                    
                else:
                    #array internal state has to be collapsed for this calculation
                    if not sample._is_collapsed:
                        msg = f"<frequency_summary> routine will collapse sample {sample.name}'s internal state"
                        self.logger.info(msg)
                        sample._collapse_internal_state()
                    
                    #initialize the frequency array: 3D array to be reduced along axis 0 at the end
                    maxlen = self._find_max_len(design, loc)
                    freq = np.zeros((len(design), len(tokens), maxlen), dtype=np.float32)
                    
                    for i,template in enumerate(design):                    
                        
                        row_mask = sample._internal_state[:,i]
                        col_mask = template(loc, return_mask=True)
    
                        #calculated weighed contributions of each design
                        #to the overall frequency array
                        norm = np.divide(np.sum(row_mask), arr.shape[0])
                        freq[i,:,:len(col_mask)] = norm * np.nan_to_num(get_freqs(arr[row_mask][:,col_mask], tokens))
    
                    #reduce back to a 2D array and plot/save
                    freq = np.sum(freq, axis=0)
                
                if loc == 'all':
                    nloc = 'overall'
                    fname =f'{sample.name}_{where}_overall_tokenwise_frequency'
                else:
                    nloc =  ', '.join(str(x + 1) for x in loc)
                    fname = f'{sample.name}_{where}_reg{nloc}_tokenwise_frequency'
                
                destination = os.path.join(self.dirs.parser_out, sample.name)
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.tokenwise_frequency(freq, tokens, where, nloc, basename)  

                if save_txt:
                    
                    np.savetxt(basename + '.csv',
                               freq,
                               delimiter=',')                
                    
            return data
        return frequency_summary

    def q_summary(self, loc=None, save_txt=False):
        '''
        For each sample in Data, compute some basic Q score statistics.
    	For each position in regions specified by 'loc', computes the mean and standard deviation
        of Q scores. Plots the results in the parser output folder as specified by config.py.
        Optionally, the data can also be written to a txt file.
        	    	
        Parameters:					  
                     loc: a list of ints to specify regions to be analyzed;
                          in this case, the op will collapse sample's internal
                          state (see explanation for Data objects)
                          
                          OR
                          
                          'all': to get the same statistics over the entire Q 
                                 score arrays. in this case, the op will NOT 
                                 collapse sample's internal state
                          
                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots						 
                          						  							  
        Returns:
                Data object (no transformation)	
        '''
        
        if loc != 'all':
            self._loc_check(loc, self.D_design)
            
        def q_score_summary(data):
            
            self._prepare_destinations(data)
            for sample in data:
                self._transform_check(sample, inspect.stack()[0][3])
                
                if loc == 'all':
                    relevant_arr = sample.Q.astype(np.float32)
                    
                else:
                    arr = sample.Q
                    if not sample._is_collapsed:
                        msg = f"<q_score_summary> routine will collapse sample {sample.name}'s internal state"
                        self.logger.info(msg)
                        sample._collapse_internal_state()
                                    
                    maxlen = self._find_max_len(self.D_design, loc)
                    #iterate over templates and append all of the relevant arr views to this array
                    #relevant view: masked (row/columnwise) arr
                    relevant_arr = []
                    
                    for i,template in enumerate(self.D_design):
                        
                        row_mask = sample._internal_state[:,i]
                        col_mask = template(loc, return_mask=True)                
                        
                        arr_view = np.zeros((np.sum(row_mask), maxlen), dtype=np.float32)
                        arr_view[:,:len(col_mask)] = arr[row_mask][:,col_mask]
                        relevant_arr.append(arr_view)
    
                    #assemble into a single array
                    relevant_arr = np.vstack(relevant_arr)
                    
                #mask out pads (0) as nans for nanmean/nanstd statistics
                relevant_arr[relevant_arr == 0] = np.nan
                
                #get the stats; plot
                q_mean = np.nanmean(relevant_arr, axis=0)
                q_std = np.nanstd(relevant_arr, axis=0)

                destination = os.path.join(self.dirs.parser_out, sample.name)
                if loc == 'all':
                    nloc = 'overall'
                    fname = f'{sample.name}_overall_q_score_summary'
                else:
                    nloc =  ', '.join(str(x + 1) for x in loc)
                    fname = f'{sample.name}_reg{nloc}_q_score_summary'
                    
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.Q_score_summary(q_mean, q_std, nloc, basename)  

                if save_txt:
                    q = np.vstack((q_mean, q_std))
                    np.savetxt(basename + '.csv',
                               q.T,
                               delimiter=',',
                               header='Q mean, Q std')  
                
            return data
        return q_score_summary

    def count_summary(self, where=None, top_n=None, fmt=None):
        '''
        For each sample in Data, counts the number of times each unique sequence is found in the
        dataset specified by 'where'. The results are written to a file in the parser output  
        folder as specified by config.py.
        
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
            msg = f"<count_summary> routine received invalid fmt argument. Acceted any of ('csv', 'fasta'); received: {fmt}"
            self.logger.error(msg)
            raise ValueError(msg)
            
        if top_n is not None:
            if not isinstance(top_n, int):
                msg = f'<count_summary> routine expected to receive parameter top_n as as int; received: {type(top_n)}'
                self.logger.error(msg)  
                raise ValueError(msg)

        def _writer(sample, og_ind, counts, fmt, path):
            if fmt == 'csv':           
                df = pd.DataFrame(columns=['Peptide', f'{where} count', 'DNA'])
                df['Peptide'] = [''.join(x) for x in sample.P[og_ind]]
                df['DNA'] = [''.join(x) for x in sample.D[og_ind]]
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
        
        def full_count_summary(data):
            self._prepare_destinations(data)
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]            
                
                #count entries in the array
                unique, og_ind, counts = np.unique(arr, axis=0, 
                                                   return_counts=True,
                                                   return_index=True)                
                
                #if top_n is unset, ind array will index every entry in the sample
                ind = np.argsort(counts)[::-1][:top_n]

                og_ind = og_ind[ind]
                counts = counts[ind]
                
                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f'{sample.name}_{where}_count_summary'
                path = os.path.join(destination, fname)
                
                _writer(sample, og_ind, counts, fmt, path)
                                
            return data
        return full_count_summary
    
    def template_summary(self, where=None):
        '''
        For each sample in Data, compute the number of matches between the dataset 
        specified by 'where' and the corresponding library templates. The results 
        are written to a file in the parser output folder as specified by config.py.
        
        In other words, summarize where dataset sequences come from (from which
        libraries). The op could also be called "_internal_state_summary"
        
        Parameters:					  
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on                        						  							  
        Returns:
                Data object (no transformation)
        '''        
        
        self._where_check(where)
        if where == 'pep':
            design = self.P_design
        
        elif where == 'dna':
            design = self.D_design
            
        def template_breakdown(data):
            self._prepare_destinations(data)
       
            #summarize straight into a pandas dataframe
            sample_names = [sample.name for sample in data]
            templates = [template.lib_seq for template in design]
        
            #all this op is: axis=0-wide sum of the internal states
            import pandas as pd
            df = pd.DataFrame(index=sample_names, columns=templates)
            
            for sample in data:
                self._transform_check(sample, inspect.stack()[0][3])
                df.loc[sample.name] = np.sum(sample._internal_state, axis=0)
                
            fname = f'{self.exp_name}_by_template_breakdown.csv'
            path = os.path.join(self.dirs.logs, fname)            
            df.to_csv(path + '.csv', sep=',')
    
            return data
        return template_breakdown

    def tSNE_analysis(self, where=None, top_n=1000, cluster_fasta=False):
        '''
        For each sample in Data, compute tSNE embeddings for the dataset 
        specified by 'where'. Cluster the results (HDBSCAN) and summarize. 
        
        Parameters:					  
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on

                   top_n: the number of top entries (by count) to analyze   

           cluster_fasta: True/False. if True, create a .fasta file for
                          sequences comprising each cluster. 
          						  						   
        Returns:
                Data object (no transformation)
        '''        

        self._where_check(where)
        if not isinstance(top_n, int):
            msg = f'<tSNE_analysis> routine expected to receive parameter top_n as as int; received: {type(top_n)}'
            self.logger.error(msg)  
            raise ValueError(msg)        

        if not isinstance(cluster_fasta, bool):
            msg = f'<tSNE_analysis> routine expected to receive parameter cluster_fasta as as bool; received: {type(cluster_fasta)}'
            self.logger.error(msg)  
            raise ValueError(msg)        

        try:
            from sklearn.manifold import TSNE
            from sklearn.preprocessing import StandardScaler
            import hdbscan
        except:
            msg = 'Failed to import the libraries necessary for <tSNE_analysis> op. . .'
            self.logger.error(msg)
            raise ImportError(msg)

        def _cluster_writer(df, destination, sname):
            clusters = np.unique(df['cluster'])
            for cluster in clusters:
                fname = f'{sname}_{where}_HDBSCAN_cluster_{cluster}.fasta'
                path = os.path.join(destination, fname)
                
                with open(path, 'w') as f:
                    for entry in df[df['cluster'] == cluster].iterrows():
                        f.write(f'>cluster_{cluster}_count_{entry[1]["counts"]}\n')  
                        f.write(f'{entry[1]["sequence"]}\n')
            return        

        def tSNE_embedding(data):
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]
                
                #get an array of top_n entries, X; C - counts array
                X, C = np.unique(arr, return_counts=True, axis=0)
                ind = np.argsort(C)[::-1][:top_n]
                X = X[ind]
                C = C[ind]
                
                #get one-hot encodings of the datasets
                tokens = np.unique(arr)
                nX = np.vectorize(lambda x: np.where(tokens == x)[0][0])(X)
                nX = nX.astype(np.int8).ravel()
                
                fX = np.zeros((X.size, tokens.size))
                fX[np.arange(nX.size), nX] = 1
                fX = np.reshape(fX, (X.shape[0], -1))
                
                #embed the representations and scale
                embedder = TSNE(n_components=2,
                                perplexity=30, 
                                learning_rate=200,
                                init='random')
                
                Y = embedder.fit_transform(fX.astype(np.float32))
                Y = StandardScaler().fit_transform(Y)
                
                #run HDBSCAN clustering
                hdb = hdbscan.HDBSCAN(cluster_selection_epsilon=0.02,
                                      cluster_selection_method = 'eom')
                hdb.fit(Y)
                
                #plot the results
                destination = os.path.join(self.dirs.parser_out, sample.name, 'tSNE_analysis')
                if not os.path.isdir(destination):
                    os.makedirs(destination)       
                    
                fname = f'{sample.name}_{where}_tSNE'
                basename = os.path.join(destination, fname)         
                sizes = 8000 * np.power(np.divide(C, arr.shape[0]), 0.6)
                Plotter.Analysis.tSNE(Y, sizes, hdb.labels_, basename)
                
                #dump analysis into a .csv
                entries = np.array([''.join(x) for x in X])
                clusters = hdb.labels_ + 1
                
                d = {'sequence': entries,
                     'counts': C,
                     'cluster': clusters}
                
                df = pd.DataFrame.from_dict(d)
                df = df.sort_values(['cluster', 'counts'], ascending=False)
                fname = f'{sample.name}_{where}_HDBSCAN_clustering.csv'
                full_name = os.path.join(destination, fname)
                df.to_csv(full_name, sep=',', index=False)
                
                if cluster_fasta:
                    _cluster_writer(df, destination, sample.name)
                
            return data
        return tSNE_embedding

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
                                  D=DNA,
                                  Q=Q,
                                  P=None
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
            msg = f"<save_data> routine received invalid fmt argument. Acceted any of ('npy', 'csv', 'fasta'); received: {fmt}"
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
                    
            self._prepare_destinations(data)
            for sample in data:
                
                self._transform_check(sample, inspect.stack()[0][3])
                arr = sample[where]  

                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f'{sample.name}_{where}'
                path = os.path.join(destination, fname)
                
                _writer(arr, fmt, path)
                
            return data
        return save_data    