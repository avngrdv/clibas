# -*- coding: utf-8 -*-
"""
Created on Sat Aug 8 20:28:10 2021
@author: Alex Vinogradov
"""

if __name__ == '__main__':
    
    #import prerequisities
    from utils.ProcessHandlers import Pipeline, FastqParser
    from utils.Dispatcher import Dispatcher
    
    #config file holds the information about library designs and
    #other parser instructions (where to look for data, where to save results etc)
    import config
    
    #initialize a dispatcher object; dispatcher is strictly speaking
    #not necessary, but it simplifies initialization of data handlers
    dispatcher = Dispatcher(config)    
    
    #a list of handlers to initialize; pipeline should always be included
    #if NGS data parsing is the goal, FastqParser will do most of the work
    handlers = (Pipeline, FastqParser)
    
    #initialize the handlers
    pip, par = dispatcher.dispatch_handlers(handlers)
        
    #enqueue the list of ops to run
    #note that at this stage no data processing will take place,
    #but the validity of specifications will be asserted.
    pip.enque([
                par.fetch_gz_from_dir(), 
                par.translate(),         
                par.filt_ambiguous(where='pep'),
                par.count_summary(where='pep', top_n=100, fmt='csv'),
                par.count_summary(where='pep', top_n=100, fmt='fasta')
             ])
    
    #this will execute the pipeline
    #if save_summary=True, summary will be saved in 
    #the logs folder as specified in the config file
    data = pip.run(save_summary=True)
        
        
    
    
    
    
    

















