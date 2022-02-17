# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 17:08:16 2022
@author: Alexander Vinogradov
"""

import os
from clibas.baseclasses import Handler
import numpy as np
from clibas.datatypes import Data, AnalysisSample
import clibas.featurize as featurize

class DataPreprocessor(Handler):
    '''
    A processor for assembling machine learning training/testing data
    from peptide sequence lists. DataHandler, i.e., every public routine
    acts on Data, transforms it, and returns the transformed version.
    
    The class is a container for ops commonly used for dataset preprocessing.    
    '''
    
    def __init__(self, *args):
        super(DataPreprocessor, self).__init__(*args)    
        self._on_completion()
        return

    def __repr__(self):
        return '<DataPreprocessor object>'

    def label_sequencing_samples(self, sele_name=None, anti_name=None):
        '''
        Label sequencing samples and convert them to TrainingSample instances
        All peptides from SequencingSample with name pos_name will be labelled 
        as 1; neg_name samples will be labelled as 0. 
    
        Parameters:
                sele_name: str. The name of the selection dataset
                anti_name: str. The name of the antiselection dataset
                                  
        Returns:
                Transformed Data object containing SequencingSample instances
        ''' 
        
        if not isinstance(sele_name, str):
            msg = '<label_sequencing_samples> op expected param sele_name as type=str; received: {type(sele_name)}'
            self.logger.error(msg)
            raise ValueError(msg)    
            
        if not isinstance(anti_name, str):
            msg = '<label_sequencing_samples> op expected param anti_name as type=str; received: {type(anti_name)}'
            self.logger.error(msg)
            raise ValueError(msg)    
            
        def label_datasets(data):
            
            #fetch datasets and reassemble Data
            
            X_pos = data[sele_name].P
            X_neg = data[anti_name].P
            
            y_pos = np.ones(X_pos.shape[0])
            y_neg = np.zeros(X_neg.shape[0])
            
            data = Data([
                         AnalysisSample(X=X_pos, y=y_pos, name='pos'),
                         AnalysisSample(X=X_neg, y=y_neg, name='neg'),
                       ])
            
            return data 
        return label_datasets


    def int_repr(self, alphabet=None):
        '''
        Transform whatever X dataset representation that is passed to 
        the "numerical representation", that is a X-array (shape: num_entries, max_len),
        where monomers are represented as integers. Pad tokens (usually '' for 
        string-type arrays) will be represented as -1, which is the de facto 
        standard int pad token mapping for the package.
        
        This is the preferred representation for most ops in the preprocessor.
        
        Parameters:
                alphabet:     current token alphabet; dtype = tuple, list or
                              np.ndarray are all OK.
    
        Returns:
                Transformed Data object
        '''
        
        if alphabet is None or not isinstance(alphabet, (tuple, list, np.ndarray)):
            msg = f'<int_repr> op expected param alphabet as type=(list, tuple, np.ndarray); received: {type(alphabet)}'
            self.logger.error(msg)
            raise ValueError(msg)
        
        #decided not to infer the alphabet: this may result in different mappings
        #for different samples in the dataset, which is definitely not ideal, because
        #the samples may get merged down the line.
        # else:
        #     msg = f'No token alphabet was specified for <X_as_int> op; will attempt to infer it. . .'
        #     self.logger.warning(msg)
                
        def arr_as_int(arr, alphabet):
            #New implementation - real fast!
            expanded = False
            pad = arr.dtype.type()    
            
            if pad not in alphabet:
                alphabet = np.hstack((alphabet, pad))
                expanded = True       
                
            sort_idx = np.argsort(alphabet)
            out = np.searchsorted(alphabet, arr, sorter = sort_idx)
            
            if expanded:
                out = out - 1
        
            return out.astype(np.int16)
        
        def X_as_int(data):
            for sample in data:
                
                if np.issubdtype(sample.X.dtype, np.character):
                    sample.X = arr_as_int(sample.X, alphabet)
                elif np.issubdtype(sample.X.dtype, np.int):
                    msg = f"<X_as_int>: sample {sample.name}'s X is already int type; op ignored. . ."
                    self.logger.warning(msg)
                else:
                    msg = f"<X_as_int> op could not repsent sample {sample.name}'s X as int; X datatype: {sample.X.dtype}"
                    self.logger.error(msg)
                    raise ValueError(msg)
                    
            return data
        return X_as_int

    def token_filter(self, alphabet=None, tokens_to_filter_by=None):
        '''
        Every X entry containing tokens specified in tokens_to_filter_by 
        will be removed from each dataset in Data. Entries should be 
        numerically represented (see self.int_repr).
        
        Parameters:
                tokens_to_filter_by: a list containing single-letter encoded
                                     tokens.
                                     
        Returns:
                Transformed Data object
        '''

        if tokens_to_filter_by is not None:
            if not isinstance(tokens_to_filter_by, (tuple, list, np.ndarray)):
                msg = '<filter_by_token> op expected param tokens_to_filter_by as type=(list, tuple, np.ndarray); received: {type(tokens_to_filter_by)}'
                self.logger.error(msg)
                raise ValueError(msg)                
        else:
                msg = '<filter_by_token> op expected aas_to_filter_by argument.'
                self.logger.error(msg)
                raise ValueError(msg)                  
            
        def filter_by_token(data):
            for sample in data:
                            
                if np.issubdtype(sample.X.dtype, np.int):
                    msg = 'X dataset dtype should be np.int type for <filter_by_token> op.'
                    self.logger.error(msg)
                    raise ValueError(msg)                    
                
                for token in tokens_to_filter_by:
                    
                    ind = np.all(sample.X != token, axis=1)
                    sample.ind_filter(ind)
                        
            return data
        return filter_by_token

    def pop_intraset(self):
        '''
        Remove duplicates within each sample in data. Duplicates are 
        looked up and popped in X datasets; corresponding y set labels
        are also removed. Entries are resorted during the process.
        
        Parameters:
                None
                                  
        Returns:
                Transformed Data object
        '''        

        def pop_intraset_duplicates(data):
            for sample in data:
                ind = np.unique(sample.X, axis=0, return_index=True)[1]
                sample.ind_filter(ind)               
            return data
        return pop_intraset_duplicates

    def pop_interset(self):
        '''
        If an entry is found in X1 and X2, pop it from X1 and X2.
        The function only makes sense if data contains two X sets.
        Hence, a warning will be issued if more/less than 2 datasets
        are passed, and the routine will be ignored.
        
        Note: the implementation is very fast, but it breaks if the 
              arrays are too big (particularly, their -1 axis)
              TODO: rewrite in a more reliable way    
              TODO: why should this be limited to 2 datasets? Generalize
              
        Parameters:
                None
                                  
        Returns:
                Transformed Data object
        '''        

        def pop_interset_duplicates(data):
            
            if len(data) != 2:
                msg = f'Exactly two X datasets should be passed to the <pop_interset_duplicates> routine. Passed {len(data)}. Operation will be ignored.'
                self.logger.warning(msg)
                return data
            
            if data.samples[0].X.dtype != np.int8 or data.samples[1].X.dtype != np.int8:
                msg = f'<pop_interset_duplicates> routine expected datasets of dtype=np.int8. Received: {data.samples[0].X.dtype, data.samples[1].X.dtype}. Operation will be ignored.'
                self.logger.warning(msg)
                return data                
    
            #Looks ridiculous, I know, but it's fast!
            cumdims = (np.maximum(data.samples[0].X.max(), data.samples[1].X.max()) + 1) ** np.arange(data.samples[0].X.shape[1])
            X1_ind = ~np.in1d(data.samples[0].X.dot(cumdims), data.samples[1].X.dot(cumdims))
            X2_ind = ~np.in1d(data.samples[1].X.dot(cumdims), data.samples[0].X.dot(cumdims))
               
            data.samples[0].ind_filter(X1_ind)
            data.samples[1].ind_filter(X2_ind)
                 
            return data
        return pop_interset_duplicates

    def pop_validation(self, min_hd=None):
        '''
        Remove validation set sequences from training/testing data, if any.
        All training/test entries must be at least Hamming distance=min_hd 
        away from any validation sequence.  
    
        Parameters:
                min_hd: int, Hamming distance threshold
                                  
        Returns:
                Transformed Data object
        '''     

        from utils.misc import hamming_distance
        
        if not hasattr(self, 'X_val'):
             msg = 'DataPreprocessor has no information about X_val to run <pop_validation_peptides> op.'
             self.logger.error(msg)
             raise ValueError(msg)
         
        if not isinstance(min_hd, int):
            msg = '<pop_validation_peptides> op expected min_hd as type=int; received: {type(min_hd)}'
            self.logger.error(msg)
            raise ValueError(msg)
             
        def pop_validation_peptides(data):
            for sample in data:            
                for pep in self.X_val:
                    
                    to_pop = hamming_distance(sample.X, pep, min_hd, cum=True, return_index=True)
        
                    ind = np.ones(sample.size, dtype=bool)
                    ind[to_pop] = False
                    sample.ind_filter(ind)
            
            return data
        return pop_validation_peptides
    
    #TODO: fix: will fail for arrays of different widths
    def merge(self):
        '''
        Merge all datasets in Data into a single X/y set.
    
        Parameters:
                None
                                  
        Returns:
                Transformed Data object
        '''    
        
        def merge_datasets(data):
            mX = []
            my = []
            for sample in data:
                mX.append(sample.X)
                if sample.y.ndim > 0:
                    my.append(sample.y)
            
            mX = np.vstack(mX)
            if my:
                my = np.hstack(my)
            else:
                my = None
            
            return Data(samples=[AnalysisSample(X=mX, y=my, name='merged_dataset')])
        return merge_datasets

    def sample(self, sample_size=None):
        '''
        Sample from datasets in data. The op will be performed over 
        every sample in the dataset.
    
        Parameters:
                sample_size: int or float. if sample_size <= 1, 
                             it is interpreted as a FRACTION of 
                             the dataset to keep.
                             
                             if self.sample is >1, it is interpreted 
                             as the number of peptides to sample.
                                  
        Returns:
                Transformed Data object
        '''           
        
        if not (isinstance(sample_size, int) or
                isinstance(sample_size, float)
               ):
            
            msg = '<sample_from_datasets> routine expected sample_size as type=int or float; received: {type(sample_size)}'
            self.logger.error(msg)
            raise ValueError(msg)
       
        def sample_from_datasets(data):
            for sample in data:
                
                if sample_size <= 1:
                    size = int(sample_size * sample.size)
                else:
                    size = int(sample.size)
                    
                if sample.size < size:
                    msg = f'Cannot take a sample that is bigger than the dataset. Sampling is ignored for {sample} sample.'
                    self.logger.warning(msg)
                    continue
                    
                ind = np.random.choice(sample.size, size=size, replace=False)
                sample.ind_filter(ind) 
                    
            return data
        return sample_from_datasets


    def shuffle(self):
        '''
        Reshuffle datapoints inside each dataset while keeping 
        theX/y mappings.
    
        Parameters:
                None
                                  
        Returns:
                Transformed Data object
        '''

        def shuffle_intraset(data):
            for sample in data:
                
                ind = np.arange(sample.X.shape[0])
                np.random.shuffle(ind)
                sample.ind_filter(ind)
                            
            return data
        return shuffle_intraset


    def tt_split(self, test_fraction=None):
        '''
        Perform test/train split.
        Data should contain a single sample at this stage.
    
        Parameters:
                test_fraction: int or float. The fraction of data
                               that will go to the test dataset. Everything
                               else to the train set.
                                  
        Returns:
                Transformed Data object
        '''
        
        if not (isinstance(test_fraction, int) or
                isinstance(test_fraction, float)
               ):
            
            msg = '<test_train_split> routine expected test_fraction as type=int or float; received: {type(test_fraction)}'
            self.logger.error(msg)
            raise ValueError(msg)
            
        def test_train_split(data):
            if data.size != 1:
                msg = 'A single dataset should be passed to the test_train_split routine. Operation is ignored.'
                self.logger.warning(msg)
                return data    
            
            full_set_size = data.samples[0].size
            test_set_size = int(test_fraction * full_set_size)
            test_set_ind = np.random.choice(full_set_size, size=test_set_size, replace=False)
            mask = np.ones(full_set_size, dtype=bool)
            mask[test_set_ind] = False
                  
            X_test = data.samples[0].X[~mask]
            X_train = data.samples[0].X[mask]
            
            #TODO: rewrite to make sure that whatever other arrays samples may contain
            #are not left behind
            if data.samples[0].y.ndim > 0:
                y_test = data.samples[0].y[~mask]
                y_train = data.samples[0].y[mask]
            else:
                y_test = None
                y_train = None
    
            data = Data(samples=[
                                 AnalysisSample(X=X_train, y=y_train, name='train_data'),
                                 AnalysisSample(X=X_test, y=y_test, name='test_data')
                                ])
            
            return data
        return test_train_split

    def to_h5(self, reshape=False, alphabet=None, repad=False, chunks=None, return_data=False):
        '''
        Featurize X datasets in Data to hdf5 files. Good when the featurized
        datasets don't fit the memory.
        
        Featurization lookup matrix must be specified in config.
        If not, one hot representations will be produced.
    
        Parameters:
            
                reshape: Should be flagged if a peptide sequence is to be 
                         represented as a mutltidimensional tensor. if reshape
                         is set False, the peptide representation will be unrolled
                         into a vector
                         
                         REPADDED MATRIX SHOULD NOT BE RESHAPED 
                         (it doesn't make sense but mathematically will work)
            
               alphabet: Current representations of X             
                          
                  repad: True/False. Should be flagged True when F is internally padded.
                         Many representations, for instance one hot, have equally
                         long vectors corresponding to each amino acids, but some,
                         for instance, SMILES_repr_v2 are all different.
                         In that case, the SMILES_repr_v2 matrix is internally padded
                         to the longest representation, which upon mapping to X will
                         result in pads in the middle of the sequence. repadding
                         will push all pads to the right
                                       
                 chunks: int. If featurized datasets don't fit the memory,
                         featurize data in chunks
                         
                         
            return_data: True/False. if True return untrasnformed Data object back
                         if False, returns None.
                         
                         
        Returns:
                Data or None
        '''        

        if not hasattr(self, 'F'):
            msg = 'DataPreprocessor did not receive an F matrix for featurization. Featurization will result in one-hot encoding. . .'
            self.logger.warning(msg)
            self.F = None
        
        if not isinstance(reshape, bool):
            msg = '<featurize_to_h5> op expected param reshape as type=bool; received: {type(reshape)}'
            self.logger.error(msg)
            raise ValueError(msg)
        
        if not isinstance(repad, bool):
            msg = '<featurize_to_h5> op expected param repad as type=bool; received: {type(repad)}'
            self.logger.error(msg)
            raise ValueError(msg)        
        
        if not isinstance(chunks, int):
            msg = '<featurize_to_h5> op expected param chunks as type=int; received: {type(chunks)}'
            self.logger.error(msg)
            raise ValueError(msg)         

        if not isinstance(return_data, bool):
            msg = '<featurize_to_h5> op expected param return_data as type=bool; received: {type(return_data)}'
            self.logger.error(msg)
            raise ValueError(msg)    

        if alphabet is not None:
            if not isinstance(alphabet, (tuple, list, np.ndarray)):
                msg = '<featurize_to_h5> routine expected param alphabet as type=(list, tuple, np.ndarray); received: {type(alphabet)}'
                self.logger.error(msg)
                raise ValueError(msg)        

        def featurize_to_h5(data):
            for sample in data:
    
                if not sample.y.ndim > 0:
                    msg = f'No y dataset for {sample.name} sample. Featurizing X only. . .'
                    self.logger.info(msg)
    
                path = os.path.join(self.dirs.ml_data, f'{sample.name}.hdf5')
                featurize.into_h5(sample.X, 
                                  y=sample.y,
                                  alphabet=alphabet,
                                  path=path,
                                  F=self.F, 
                                  reshape=reshape,
                                  repad=repad,
                                  chunks=chunks
                                 )
                
            if return_data:
                return data
            
            return
        return featurize_to_h5

    def featurize_X(self, reshape=False, alphabet=None, repad=False):
        '''
        Featurize X datasets in Data .
        Featurization lookup matrix must be specified in config.
        If not, one hot representations will be produced.
        
        The op can be used if the resulting datasets fit the memory. If not,
        self.to_h5 should be used instead.
    
        Parameters:
            
                reshape: Should be flagged if a peptide sequence is to be 
                         represented as a mutltidimensional tensor. if reshape
                         is set False, the peptide representation will be unrolled
                         into a vector
                         
                         REPADDED MATRIX SHOULD NOT BE RESHAPED 
                         (it doesn't make sense but mathematically will work)
            
               alphabet: Current representations of X              
            
                  repad: True/False. Should be flagged True when F is internally padded.
                         Many representations, for instance one hot, have equally
                         long vectors corresponding to each amino acids, but some,
                         for instance, SMILES_repr_v2 are all different.
                         In that case, the SMILES_repr_v2 matrix is internally padded
                         to the longest representation, which upon mapping to X will
                         result in pads in the middle of the sequence. repadding
                         will push all pads to the right

        Returns:
                Transformed Data object
        '''        

        if not hasattr(self, 'F'):
            msg = 'DataPreprocessor did not receive an F matrix for featurization. Featurization will result in one-hot encoding. . .'
            self.logger.warning(msg)
            self.F = None
        
        if not isinstance(reshape, bool):
            msg = '<featurize_X_datasets> op expected param reshape as type=bool; received: {type(reshape)}'
            self.logger.error(msg)
            raise ValueError(msg)
        
        if not isinstance(repad, bool):
            msg = '<featurize_X_datasets> op expected param repad as type=bool; received: {type(repad)}'
            self.logger.error(msg)
            raise ValueError(msg)
            
        if alphabet is not None:
            if not isinstance(alphabet, (tuple, list, np.ndarray)):
                msg = '<featurize_X_datasets> op expected param alphabet as type=(list, tuple, np.ndarray); received: {type(alphabet)}'
                self.logger.error(msg)
                raise ValueError(msg)        
    
        def featurize_X_datasets(data):
            for sample in data:
                sample.X = featurize.from_matrix_v3(
                                                    sample.X,
                                                    F=self.F,
                                                    alphabet=alphabet,
                                                    repad=repad,
                                                    reshape=reshape
                                                   )
            
            return data
        return featurize_X_datasets

    def drop(self, dataset_to_drop=None):
        '''
        Drop a dataset from Data.
    
        Parameters:
                dataset_to_drop: str. The name of the dataset to drop
                                  
        Returns:
                Transformed Data object
        ''' 

        if not isinstance(dataset_to_drop, 'str'):
            msg = '<drop_dataset> op expected param dataset_to_drop as type=str; received: {type(dataset_to_drop)}'
            self.logger.error(msg)
            raise ValueError(msg)   
 
        def drop_dataset(data):
            to_drop = []
            for i,sample in enumerate(data):
                if sample.name == dataset_to_drop:
                    to_drop.append(i)
                    
            if not to_drop:
                msg = f'<drop_dataset>: {dataset_to_drop} dataset specified for dropping could not be identified. . .'
                self.logger.warning(msg)
    
            for i in to_drop:                      
                del(data.samples[i])
            
            return data
        return drop_dataset
    
