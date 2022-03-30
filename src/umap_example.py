# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 19:36:14 2022
@author: av
"""

import config
import numpy as np

from clibas.dispatchers import Dispatcher
from clibas.dataanalysis import HDBUMAP
from clibas.datatypes import Data, AnalysisSample

''' 
Load pre-parsed .npy files as shown below. For this to work as shown below,
the files need to located in the same directory as this script. Optional, but
recommended is to append names to the individual samples.
'''

arr1 = np.load('Pf4_S46_L001_R1_001.fastq_pep.npy')
arr2 = np.load('Pf5_S47_L001_R1_001.fastq_pep.npy')
arr3 = np.load('Pf6_S48_L001_R1_001.fastq_pep.npy')
arr4 = np.load('Pf7_S49_L001_R1_001.fastq_pep.npy')

data = Data([
              AnalysisSample(X=arr1, name='Harry DUB round 4'),
              AnalysisSample(X=arr2, name='Harry DUB round 5'),
              AnalysisSample(X=arr3, name='Harry DUB round 6'),
              AnalysisSample(X=arr4, name='Harry DUB round 7'),
            ])

hdu, = Dispatcher(config).dispatch(HDBUMAP)

''' 
Below are hyperparameters for sequence embedding + clustering. 
These can be changed to tailor the results. I would say single_manifold and 
top_n params are the most important, everything else can be left as is 
(if the libraries don't have unusual amino acids, otherwise the config file
needs to be edited).
'''

analyzer = hdu.analysis(F='pep_ECFP4',
                        where='pep',
                        top_n=2000,
                        single_manifold=True,
                        return_modified=True,
                        )

data = analyzer(data)