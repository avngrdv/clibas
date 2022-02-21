# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 20:27:57 2021
@author: Alex Vinogradov
"""

import os, importlib
from clibas.baseclasses import (
                                  Logger,
                                  DirectoryTracker
)


_reserved_aa_names = (
                      '_',
                      '+', 
                      '*', 
                      '1',
                      '2', 
                      '3',
                      '4', 
                      '5', 
                      '6', 
                      '7', 
                      '8',
                      '9', 
                      '0'
)

class Dispatcher:
    '''
    An object used to coordinate and configure data handlers. Its primary 
    purpose is to parse config files and use the specified information to  
    dispatch various data handlers in a coordinated way. 
    '''
    
    def __init__(self, config):
        self.config = config
        self._parse_config()
        return

    def __repr__(self):
        return f'<Dispatcher object at {hex(id(self))}>'
    
    def _parse_dir_config(self):
        if not hasattr(self.config, 'TrackerConfig'):    
            self.dirs = DirectoryTracker()
        
        else:
            #don't have to check individual attribs becase DirectoryTracker
            #can handle it by itself
            self.dirs = DirectoryTracker(self.config.TrackerConfig)
        return
        
    def _parse_logger_config(self):

        #logger: just need to set up the log file location, if needed
        if hasattr(self.config, 'LoggerConfig'):                   
            if hasattr(self.config.LoggerConfig, 'log_to_file'):
                if self.config.LoggerConfig.log_to_file:
                    if not hasattr(self.config.LoggerConfig, 'log_fname'):
                        self.config.LoggerConfig.log_fname = os.path.join(
                        
                            self.dirs.logs, self.config.experiment + '_logs.txt'
                            )
        
            self.L = Logger(config=self.config.LoggerConfig).logger

        else:
            self.L = Logger()        

        return

    def _parse_constants_config(self):    

        #setup constants: some missing attribs might be inferred, but
        #at least something needs to be present
        if not hasattr(self.config, 'constants'):    
            msg = '<Dispatcher>: config is missing the minimal definitions of constants. Aborting. . . '
            self.L.error(msg)
            raise ValueError(msg)
                    
        self.constants = self.config.constants
        #deal with amino acids
        #amino acids can be inferred from the translation table
        if not hasattr(self.constants, 'aas'):
            if hasattr(self.constants, 'translation_table'):
            
                self.constants.aas = tuple(sorted(set(
                    
                                 x for x in self.constants.translation_table.values() 
                                 if x not in _reserved_aa_names))
                     )
        
        #deal with nucleic acids, if any  
        #bases can be inferred from the translation table
        if not hasattr(self.constants, 'bases'):
            if hasattr(self.constants, 'translation_table'):
            
                self.constants.bases = tuple(sorted(set(
                    
                                ''.join(self.constants.translation_table.keys())))
                      )        
        return
        
    def _parse_lib_design_config(self):
        
        from clibas.lib_design import LibraryDesign
        if not hasattr(self.config, 'LibaryDesigns'):
            msg = '<Dispatcher>: config has not specified any library designs. . .'
            self.L.warning(msg)
            return
        
        #deal with peptide libraries first
        if hasattr(self.config.LibaryDesigns, 'pep_templates') and \
           hasattr(self.config.LibaryDesigns, 'pep_monomers'):
            
            if hasattr(self.constants, 'aas'):             
                self.P_design = LibraryDesign(
                                              templates=self.config.LibaryDesigns.pep_templates,
                                              monomers=self.config.LibaryDesigns.pep_monomers,
                                              lib_type='pep',
                                              val_monomer=self.constants.aas
                                             ) 
            else:
                msg = 'Dispatcher: cannot setup a peptide library design without amino acid alphabet specification. . .'
                self.L.error(msg)
                raise ValueError(msg)
                
        #if only of the two necessary atrributes is present, raise
        #note that if both are absent, it's OK
        elif hasattr(self.config.LibaryDesigns, 'pep_templates') !=   \
             hasattr(self.config.LibaryDesigns, 'pep_monomers'):
                 
            msg = 'Dispatcher: config is missing necessary parameters to setup peptide library designs. . .'
            self.L.error(msg)
            raise ValueError(msg)

        #same, but now for DNA libs
        if hasattr(self.config.LibaryDesigns, 'dna_templates') and  \
           hasattr(self.config.LibaryDesigns, 'dna_monomers'):
            
            if hasattr(self.constants, 'bases'):
                self.D_design = LibraryDesign(
                                              templates=self.config.LibaryDesigns.dna_templates,
                                              monomers=self.config.LibaryDesigns.dna_monomers,
                                              lib_type='dna',
                                              val_monomer=self.constants.bases
                                             )
            else:
                msg = 'Dispatcher: cannot setup a DNA library design without nucleotide alphabet specification. . .'
                self.L.error(msg)
                raise ValueError(msg)
        
        #if only of the two necesasry atrributes is present, raise
        elif hasattr(self.config.LibaryDesigns, 'dna_templates') != \
             hasattr(self.config, 'dna_monomers'):
                 
            msg = 'Dispatcher: config is missing necessary parameters to setup DNA library designs. . .'
            self.L.error(msg)
            raise ValueError(msg)
            
        return
        
    def _parse_config(self):
        
        #check that all necessary config params are in place, one by one (tedious)
        #if not, fallback to some innocuous defaults where possible
        if hasattr(self.config, 'experiment'):
            self.experiment = self.config.experiment
        else: 
            self.experiment = 'untitled_experiment'
        
        self._parse_dir_config()
        self._parse_logger_config()
        self._parse_constants_config()
        self._parse_lib_design_config()
        
        self.common = {'logger': self.L,
                       'dirs': self.dirs, 
                       'exp_name': self.experiment,
                       'constants': self.constants
                      }
        return
    
    def _config_to_dict(self, conf):        
        return dict((name, getattr(conf, name)) for name 
                    in dir(conf) if not name.startswith('__'))        

    def _get_lib_params(self):
        params = dict()
        if hasattr(self, 'P_design'):
            params.update({'P_design': self.P_design})
    
        if hasattr(self, 'D_design'):
            params.update({'D_design': self.D_design})        

        return params
    
    def _dispatch_Pipeline(self):
        from clibas.pipelines import Pipeline
        return Pipeline(self.common)
    
    def _dispatch_FastqParser(self):

        from clibas.parsers import FastqParser
        params = dict()
        if hasattr(self.config, 'FastqParserConfig'):
            params = self._config_to_dict(self.config.FastqParserConfig)
            
        params.update(self._get_lib_params())
        params.update(self.common)
        return FastqParser(params)

    def _dispatch_DataPreprocessor(self):
        
        from clibas.datapreprocessors import DataPreprocessor
        if hasattr(self.config, 'PreproConfig'):
            params = self._config_to_dict(self.config.PreproConfig)  
            
        params.update(self.common)
        return DataPreprocessor(params)

    def _dispatch_DataAnalysisTools(self):
        
        from clibas.dataanalysis import DataAnalysisTools
        params = dict()
        if hasattr(self.config, 'DataAnalysisConfig'):
            params = self._config_to_dict(self.config.DataAnalysisConfig)
            
        params.update(self._get_lib_params())
        params.update(self.common)        
        return DataAnalysisTools(params)

    def _dispatch_generic(self, handler):
       
        handler_modules = ['parsers',
                           'dataanalysis',
                           'preprocessors',
                           'pipelines'
                          ]
        
        for module in handler_modules:
            try:
                m =  importlib.import_module(f'clibas.{module}')
                obj = getattr(m, handler.__name__)
                break
            except:
                obj = None
            
        if not obj: raise ImportError(f'<Dispatcher> could not set up {handler}')
         
        params = dict()
        params.update(self._get_lib_params())
        params.update(self.common)
        return obj(params)
    
    def dispatch(self, handlers):
        '''
        Take a tuple of handlers and initialize each by giving 
        them relevant configs params. Also makes sure that the
        information about file directories and logger is shared
        between all initialized objects
        
        Parameters
        ----------
        handlers : a tuple of handlers to be instantiated or an individual
                   handler
                   ex: (Pipeline, FastqParser, DataPreprocessor)

        Returns
        -------
        A tuple of initialized handler instances
        '''
        
        mapping = {'Pipeline': self._dispatch_Pipeline,
                   'FastqParser': self._dispatch_FastqParser,
                   'DataAnalysisTools': self._dispatch_DataAnalysisTools,
                   'DataPreprocessor': self._dispatch_DataPreprocessor
                  }
        
        #make sure either individual handlers or tuples of thereof work
        handlers = [handlers]
        try:
            handlers = [item for sublist in handlers for item in sublist]
        except:
            pass
     
        h = list()
        for handler in handlers:
            try:
                h.append(mapping[handler.__name__]())
            except:
                h.append(self._dispatch_generic(handler))
                
        for i in h:
            msg = f'{i} was succesfully initialized'
            self.L.info(msg)
            
        return tuple(h)
    







