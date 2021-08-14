# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 20:27:57 2021
@author: Alex Vinogradov
"""

from utils.ProcessHandlers import (
                                    Logger, 
                                    DirectoryTracker, 
                                    Pipeline, 
                                    FastqParser
)

class Dispatcher:
    '''
    An object used to coordinate and configure data handlers.
    Uses config.py to setup the handlers automatically making
    sure that different handlers share a common directory tracking
    status and a logger object. 
    '''
    
    def __init__(self, config):
        self.config = config
        self._setup()
        
    def _setup(self):
        T = DirectoryTracker(self.config.TrackerConfig)
        L = Logger(config=self.config.LoggerConfig)
        self.common = {'logger': L.logger,
                       'dirs': T, 
                       'exp_name': self.config.experiment,
                       'constants': self.config.constants}
        return
    
    def _config_to_dict(self, conf):        
        return dict((name, getattr(conf, name)) for name in dir(conf) if not name.startswith('__'))        

    #can be optimized, but whatever
    def _dispatch_pipeline(self):
        return Pipeline(self.common)
        
    def _dispatch_parser(self):
        meta = self._config_to_dict(self.config.ParserConfig)
        meta.update(self.common)
        return FastqParser(meta)

    def __repr__(self):
        return '<Dispatcher object>'
    
    def dispatch_handlers(self, handlers):
        '''
        The main class method. Takes a tuple of handlers and sets
        them up in a coordinated way.
        
        Parameters
        ----------
        handlers : a tuple of handlers to be instantiated.
                   ex: (Pipeline, FastqParser, DataPreprocessor)

        Returns
        -------
        a tuple of instantiated handlers in the same order.
        all meta parameters for individual handlers are fetched
        from the config file. logger and directory tracker are 
        shared between all instances for coordination.
        '''
        
        mapping = {Pipeline: self._dispatch_pipeline,
                   FastqParser: self._dispatch_parser,
                  }
        
        h = list()
        for handler in handlers:
            
            try:
                h.append(mapping[handler]())
            except:
                msg = f'Dispatcher failed to dispatch the following handler: {handler}. . .'
                self.common['logger'].error(msg)
                raise ValueError
    
        return tuple(h)
    







