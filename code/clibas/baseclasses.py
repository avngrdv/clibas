# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 23:16:17 2022
@author: Alex Vinogradov
"""

import logging, os
import numpy as np
from clibas.datatypes import Data, AnalysisSample

class Logger:
    '''
    A decorated version of the standard python logger object. Can be setup
    from a config file. Two main customizations are implemented: verbosity
    (whether logger messages should be printed to the running stream) and 
    log_to_file, which if set, will setup a dedicated handler to dump log
    info to a file.
    '''
    
    def __init__(self, config=None):
        self.conf = config
    
        self.__fallback()
        self.__configure_logger()
        return

    def __repr__(self):
        return f'<Logger {self.name}; verbose: {self.verbose}; log_to_file: {self.log_to_file}; level: {self.level}>'

    def __fallback(self):
        '''
        If no config is passed fallback to some innocuous defaults, which is 
        basically a silent logger.
        '''
    
        attribs = {
                    'name':'unnamed',
                    'verbose':  False,
                    'log_to_file': False,
                    'log_fname': None,
                    'level': 'ERROR'     
                  }
   
        for attr in attribs:
            
            if not hasattr(self.conf, attr):
                setattr(self, attr, attribs[attr])
            else:
                setattr(self, attr, getattr(self.conf, attr))
        
        return

    def __configure_logger(self):

        levels = {
                  'DEBUG': 10,
                  'INFO': 20,
                  'WARNING': 30,
                  'ERROR': 40,
                 }
            
        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(levels[self.level])
        formatter = logging.Formatter("[%(levelname)s]: %(message)s")

        #clear any preexisting handlers to avoid stream duplication
        self.logger.handlers.clear()

        if self.verbose:
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            console_handler.setLevel(levels[self.level])

            self.logger.addHandler(console_handler)

        if self.log_to_file:
               
            if self.log_fname is None:
                cwd = os.getcwd()
                self.log_fname = os.path.join(cwd, 'unnamed_log_file.txt')
                
            filehandler = logging.FileHandler(self.log_fname)      
            filehandler.setFormatter(formatter)
            filehandler.setLevel(levels[self.level])
                    
            self.logger.addHandler(filehandler)
        
        return
        
class DirectoryTracker:
    '''
    A simple object to keep track where the data, logs, etc should be looked 
    for. Config should be passed to specify the locations. Otherwise,
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
        Any directories left unspecified by config will be set to cwd
        '''
        
        cwd = os.getcwd()
        
        #all concerned directories
        attribs = [
                    'seq_data', 
                    'logs', 
                    'parser_out', 
                    'prepro_out', 
                    'analysis_out' 
                  ]
        
        for attr in attribs:
            
            if not hasattr(self._conf, attr):
                setattr(self, attr, cwd)
            else:
                setattr(self, attr, getattr(self._conf, attr))
           
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

    def _validate_designs(self):
        
        if hasattr(self, 'P_design') and hasattr(self, 'D_design'):
            
            if not len(self.P_design) == len(self.D_design):
                msg = 'Peptide and DNA library designs must contains the same number of templates; cannot inialize FastqProcessor. . .'
                self.logger.error(msg)
                raise ValueError(msg)

        #actually, let's not require a design to instantiate parser:
        #maybe someone just wants to load the data and run translation
        # if not hasattr(self, 'P_design') and not hasattr(self, 'D_design'):                  
        #     msg = 'Either DNA or peptide library designs need to be specified for <FastqParser> to operate'
        #     self.logger.error(msg)
        #     raise ValueError(msg)
        return
    
    def _validate_constants(self):
        #but the constants do appear necessary: can't even run translation otherwise
        if not hasattr(self, 'constants'):
            msg = 'FastqParser requires constants for setup. . .'
            self.logger.error(msg)
            raise ValueError(msg)
            
        return
    
    def _transform_check(self, sample, func):
        if not sample.get_ndims() == 2:
            raise ValueError(f'Sample {sample.name} holds arrays of unsupported dimensionality for {func} op. Expected: arrays of ndims=2, got: ndims={sample.get_ndims()}')
        
        return
    
    def _where_check(self, where):

        if where == 'pep':
            if not hasattr(self, 'P_design'):
                msg = "Peptide library design not set: cannot analyze peptide datasets without unspecifying a library design."
                self.logger.error(msg)
                raise ValueError(msg)
            
        elif where == 'dna':
            if not hasattr(self, 'D_design'):
                msg = "DNA library design not set: cannot analyze DNA datasets without unspecifying a library design."
                self.logger.error(msg)
                raise ValueError(msg)                       
         
        else:
            msg = f'The parser did not understand which dataset it should operate on. Passed value: {where}; allowed values: pep/dna.'
            self.logger.error(msg)
            raise ValueError(msg)
        return

    def _infer_design(self, where):
        if where == 'pep':
            return self.P_design
        
        if where == 'dna':
            return self.D_design

    def _infer_alphabet(self, where, alphabet):
        
        #if alphabet is specified: override any implied specifications
        #by the where argument
        if isinstance(alphabet, (list, tuple)):
            return alphabet
        
        elif isinstance(alphabet, np.ndarray):
            if alphabet.ndim == 1: 
                return alphabet
            
            else: where = ''
            
        #if no alphabet is supplied, infer
        if where == 'pep':
            return self.constants.aas
            
        elif where == 'dna':      
            return self.constants.bases
        
        msg = 'Token alphabet was not supplied or was not understood. . .'
        self.logger.error(msg)
        raise ValueError(msg)
        
    def _loc_check(self, loc, design):
        
        if not isinstance(loc, (list, tuple)):
            msg = f'The Parser expected to receive a list of region indexes to parse; received: {type(loc)}'
            self.logger.error(msg)
            raise ValueError(msg)

        if max(loc) > design.loc.max():
            msg = f'{design.lib_type} library design does not contain enough regions. Library design contains {design.loc.max() + 1} regions; specified: up to {max(loc) + 1}'
            self.logger.error(msg)
            raise AssertionError(msg)
        return

    def _prepare_destinations(self, data, root):
        
        for sample in data:
            destination = os.path.join(root, sample.name)
            if not os.path.isdir(destination):
                os.makedirs(destination)        
        return    

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
    
    def _recast_data(self, maybe_arr, where=None):
        
        if not isinstance(maybe_arr, Data):
            
            if isinstance(maybe_arr, AnalysisSample):
                return Data([maybe_arr])
            try:
                return Data([AnalysisSample(X=maybe_arr)])
            except:
                msg = f'{self} did not understand input type. \
Expected data as type: Data object or an array; received: {type(maybe_arr)}'
                self.logger.error(msg)
                raise TypeError(msg)

        else:
            if maybe_arr.dtype != AnalysisSample:
                try:
                    return Data([AnalysisSample(X=sample[where], 
                                                name=sample.name) 
                                 for sample in maybe_arr]
                               )
                except:
                    msg = f'{self} could not recast {where} samples in data. . .'
                    self.logger.error(msg)
                    raise ValueError(msg)
                    
        return maybe_arr

    def _restack_and_repad(self, data):

        #TODO: add assertions and checks
        #TODO: think about where this method belongs and how to generalize it
        #location here is temporary
        new_x = sum([s.size for s in data.samples])
        new_y = max([s.X.shape[-1] for s in data.samples])
        dtype = data[0].X.dtype
        
        stacked = np.zeros((new_x, new_y), dtype=dtype)
        i = 0
        for s in data:
            
            x, y = s.X.shape
            stacked[i:i+x, :y] = s.X
            i += x
    
            padding_arr = np.zeros((x, new_y - y), dtype=dtype)
            s.X =  np.c_[s.X, padding_arr]
            
        return stacked, data





