# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 22:18:58 2022
@author: Alex Vinogradov
"""

from utils.ProcessHandlers import Handler
import numpy as np

try:
    import umap
    import hdbscan
except:
    msg = 'Failed to import the libraries necessary for HDBSCAN/UMAP analysis. . .'
    raise ImportError(msg)

class HDBUMAP(Handler):
    
    def __init__(self, *args):
        super(HDBUMAP, self).__init__(*args)
        super(HDBUMAP, self)._on_completion()
        return    
    
    def initialize(self, data):
        '''
        Determine what will take place by uptaking the data
        '''
        pass
        #determine data dtype
        #setup hyperparameters
        #here we import hopkins and calculate it
        #then we use these values to compute arrays of hypers
        #or single values
        
        '''
        n_neighbors
        min_dist
        min_cluster_size
        min_samples
        '''
        
        '''
        config param: F - should be in prepro
        config param: top_n

        '''
    
        '''
        package data -> recast data into new instance
        
        '''
    
    '''
    write the whole fucking class? 
    DataAnalysisTools(Handler):
        inherit from handler
        preprocess the set -> need a fucking data routine for this shit
        package the data: where arg or single array
        
    
    '''
    
    
    
    
    
    