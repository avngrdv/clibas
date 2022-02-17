# -*- coding: utf-8 -*-
"""
Created on Sat Aug 8 20:28:10 2021
@author: Alex Vinogradov
"""

if __name__ == '__main__':
    pass
    # #import prerequisities
    # from utils.ProcessHandlers import Pipeline, FastqParser
    # from utils.Dispatcher import Dispatcher
    
    # #config file holds the information about library designs and
    # #other parser instructions (where to look for data, where to save results etc)
    # import config_demo as config
    
    # #initialize a dispatcher object; dispatcher is strictly speaking
    # #not necessary, but it simplifies initialization of data handlers
    # dispatcher = Dispatcher(config)    
    
    # #a list of handlers to initialize; pipeline should always be included
    # #if NGS data parsing is the goal, FastqParser will do most of the work
    # handlers = (Pipeline, FastqParser)
    
    # #initialize the handlers
    # pip, par = dispatcher.dispatch_handlers(handlers)
        
    # pip.enque([
    #             par.fetch_gz_from_dir(), 
    #             par.translate(stop_readthrough=False),         
    #             par.len_filter(where='pep'),            
    #             par.cr_filter(where='pep', loc=[0, 2], tol=6),
    #             par.filt_ambiguous(where='pep'),
    #             par.q_score_filt(minQ=20, loc=[1]),
    #             par.fetch_at(where='pep', loc=[1]),
    #             par.unpad(),
    #             par.template_summary(where='pep'),
    #             par.tSNE_analysis(where='pep', top_n=1000, cluster_fasta=True)
    #           ])
    
    # #this will execute the pipeline
    # #if save_summary=True, summary will be saved in 
    # #the logs folder as specified in the config file
    # data = pip.run(save_summary=True)
        
        
    
    
    
    
    

















