# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 22:42:26 2022
@author: Alex Vinogradov
"""

import os, inspect
import numpy as np
import pandas as pd

import clibas.plotters as Plotter
from clibas.baseclasses import Handler
from clibas.datatypes import Data, AnalysisSample

class DataAnalysisTools(Handler):
    '''
    Crafted wrappers around various data analysis tools to simplify their
    calling during the pipeline creation process.
    '''
    
    def __init__(self, *args):
        super(DataAnalysisTools, self).__init__(*args)
    
        self._validate_designs()
        self._validate_constants()
        return

    def __repr__(self):
        return '<DataAnalysisTools object>'
    
    def length_analysis(self, where=None, save_txt=False):
        
        self._where_check(where)
        la = LengthAnalysis(self.__dict__)
        op = la.len_summary(where=where, save_txt=save_txt)
        return op
    
    def q_score_analysis(self, loc=None, save_txt=False):
        
        if loc is not None:
            self._where_check('dna')
            self._loc_check(loc, self.D_design)
            
        QA = QScoreAnalysis(self.__dict__)
        op = QA.q_summary(loc=loc, save_txt=save_txt)
        return op
    
    def sequence_convergence_analysis(self, where=None, alphabet=None):

        self._where_check(where)
        alphabet = self._infer_alphabet(where, alphabet)
        CA = ConvergenceAnalys(self.__dict__)
        op = CA.sequence_level_convergence(where=where, alphabet=alphabet)
        return op
        
    def token_convergence_analysis(self, where=None, 
                                   loc=None, 
                                   alphabet=None,
                                   save_txt=False):
    
        self._where_check(where)
        design = self._infer_design(where)
        alphabet = self._infer_alphabet(where, alphabet)
        if loc is not None:
            self._loc_check(loc, design)
        
        CA = ConvergenceAnalys(self.__dict__)
        op = CA.token_level_convergence(where=where, 
                                        loc=loc, 
                                        alphabet=alphabet,
                                        save_txt=save_txt
                                       )
        return op
    
    def umap_hdbscan_analysis(self, 
                              top_n=None,
                              where=None,
                              F=None, 
                              cluster_fasta=False,
                              alphabet=None,
                              return_modified=False,
                              single_manifold=False):

        hdu = HDBUMAP(self.__dict__)
        op = hdu.analysis(top_n=top_n,
                          where=where,
                          F=F, 
                          cluster_fasta=cluster_fasta,
                          alphabet=alphabet, 
                          return_modified=return_modified,
                          single_manifold=single_manifold
                         ) 
        return op
        
class LengthAnalysis(Handler):
    
    def __init__(self, *args):
        super(LengthAnalysis, self).__init__(*args)

    def __repr__(self):
        return '<LengthAnalysis object>'    
    
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
        
        def length_summary(data):
            
            recast_data = self._recast_data(data, where=where)
            self._prepare_destinations(recast_data, self.dirs.analysis_out)        
            for sample in recast_data:
                
                arr = sample.X
                L = self._L_summary(arr)            
                L, counts = np.unique(L, return_counts=True)
                
                destination = os.path.join(self.dirs.analysis_out, sample.name)
                fname = f'{sample.name}_{where}_L_distribution'
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.L_distribution(L, counts, 
                                                      where=where, 
                                                      basename=basename
                                                     )
                
                if save_txt:
                    np.savetxt(basename + '.csv',
                               np.array((L, counts)).T,
                               delimiter=',', 
                               header='Seq length,Count')
            
            return data
        return length_summary    
    
class ConvergenceAnalys(Handler):
    
    def __init__(self, *args):
        super(ConvergenceAnalys, self).__init__(*args)
        
        from clibas.misc import (
                                 shannon_entropy,
                                 get_freqs,
                                 positional_conservation
                                )

        self.shannon_entropy = shannon_entropy
        self.get_freqs = get_freqs
        self.positional_conservation = positional_conservation
        return

    def __repr__(self):
        return '<ConvergenceAnalys object>'

    def sequence_level_convergence(self, where=None, alphabet=None):
        '''
        For each sample in Data, perform basic library convergence analysis on 
        a sequence level. Computes normalized Shannon entropy and postition-wise 
        sequence conservation. Plots the results in the parser output folder as 
        specified by config.py.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
                          
                alphabet: token alphabet for the datasets to be analyzed.
                          will be automatically inferred if 'where' is 
                          specified. otherwise, dtype: any of (list, tuple, ndarray)
                          						  							  
        Returns:
                Data object (no transformation)
        '''
        
        def sequence_level_convergence_summary(data):
    
            recast_data = self._recast_data(data, where)
            self._prepare_destinations(recast_data, self.dirs.analysis_out)
            
            for sample in recast_data:
    
                self._transform_check(sample, inspect.stack()[0][3])                
                arr = sample.X
                
                shannon, counts = self.shannon_entropy(arr, norm=True)
                freq = self.get_freqs(arr, alphabet)
                seq_conservation = self.positional_conservation(freq)
                
                destination = os.path.join(self.dirs.analysis_out, sample.name)
                fname = f'{sample.name}_{where}_library_convergence'
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.dataset_convergence(counts, shannon, where, basename)                
                
                fname = f'{sample.name}_{where}_sequence_conservation'
                basename = os.path.join(destination, fname)                
                Plotter.SequencingData.conservation(seq_conservation, where, basename)
                
            return data
        return sequence_level_convergence_summary


    def token_level_convergence(self, where=None, loc=None, alphabet=None, save_txt=False):
        '''
        Perform basic library convergence analysis at a token level. For each sample in Data, 
        computes the frequency of each token in the dataset. Plots the results in the parser 
        output folder as specified by config.py. Optionally, the data can also be written to
        a txt file.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.

                alphabet: token alphabet for the datasets to be analyzed.
                          will be automatically inferred if 'where' is 
                          specified. otherwise, dtype: any of (list, tuple, ndarray)
    						  
                     loc: if not None: a list of ints to specify regions to be 
                          analyzed; in this case, the op will collapse sample's
                          internal state (see explanation for Data objects)
                          
                          if left None: get the same statistics over the entire
                          sequence; in this case, the op will NOT collapse  
                          sample's internal state
    
                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots						 
                          						  							  
        Returns:
                Data object (no transformation)
        '''
        def token_level_convergence_analysis(data):
            
            #TODO: see whether this code can be simplified and generalized
            if loc is None:
                recast_data = self._recast_data(data, where)
            else:
                recast_data = data

            self._prepare_destinations(recast_data, self.dirs.analysis_out)            
            
            for sample in recast_data:
                
                if loc is None:
                    arr = sample.X
                    freq = self.get_freqs(arr, alphabet)
                    nloc = 'overall'
                    fname = f'{sample.name}_{where}_tokenwise_frequency'

                else:     
                    arr = sample[where]
                    design = self._infer_design(where)
                    #array internal state has to be collapsed for this calculation
                    if not sample._is_collapsed:
                        msg = f"<frequency_summary> op will collapse sample {sample.name}'s internal state"
                        self.logger.warning(msg)
                        sample._collapse_internal_state()
                    
                    #initialize the frequency array: 3D array to be reduced 
                    #along axis 0 at the end
                    maxlen = self._find_max_len(design, loc)
                    freq = np.zeros((len(design), len(alphabet), maxlen),
                                     dtype=np.float32)
                    
                    for i,template in enumerate(design):                    
                        
                        row_mask = sample._internal_state[:,i]
                        col_mask = template(loc, return_mask=True)
    
                        #calculated weighed contributions of each design
                        #to the overall frequency array
                        norm = np.divide(np.sum(row_mask), arr.shape[0])
                        freq[i,:,:len(col_mask)] = norm * np.nan_to_num(
                            
                            self.get_freqs(arr[row_mask][:,col_mask], 
                                           alphabet) 
                        )
    
                    #reduce back to a 2D array and plot/save
                    freq = np.sum(freq, axis=0)
                    nloc =  ', '.join(str(x + 1) for x in loc)
                    fname = f'{sample.name}_{where}_reg_{nloc}_tokenwise_frequency'
                
                destination = os.path.join(self.dirs.analysis_out, sample.name)
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.tokenwise_frequency(freq, alphabet, where,
                                                           nloc, basename
                                                          )
                if save_txt:          
                    np.savetxt(basename + '.csv',
                               freq,
                               delimiter=',')
                
            return data
        return token_level_convergence_analysis
    
class QScoreAnalysis(Handler):
    
    def __init__(self, *args):
        super(QScoreAnalysis, self).__init__(*args)

    def __repr__(self):
        return '<QScoreAnalysis object>'

    def q_summary(self, loc=None, save_txt=False):
        '''
        For each sample in Data, compute some basic Q score statistics.
    	For each position in regions specified by 'loc', computes the mean and standard deviation
        of Q scores. Plots the results in the parser output folder as specified by config.py.
        Optionally, the data can also be written to a txt file.
        	    	
        Parameters:					  
                     loc: if not None: a list of ints to specify regions to be 
                          analyzed; in this case, the op will collapse sample's
                          internal state (see explanation for Data objects)
                          
                          if left None: get the same statistics over the entire
                          sequence; in this case, the op will NOT collapse  
                          sample's internal state
                          
                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots						 
                          						  							  
        Returns:
                Data object (no transformation)	
        '''            
        def q_score_summary(data):
            
            if loc is None:
                recast_data = self._recast_data(data, 'Q')
            else:
                recast_data = data        
                
            self._prepare_destinations(recast_data, self.dirs.analysis_out)
            
            for sample in recast_data:
            
                if loc is None:
                    relevant_arr = sample.X.astype(np.float32)
                    nloc = 'overall'
                    fname = f'{sample.name}_q_score_summary'
                    
                else:
                    arr = sample.Q
                    if not sample._is_collapsed:
                        msg = f"<q_score_summary> op will collapse sample {sample.name}'s internal state"
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
                    
                    nloc =  ', '.join(str(x + 1) for x in loc)
                    fname = f'{sample.name}_reg{nloc}_q_score_summary'
                    
                #mask out pads (0) as nans for nanmean/nanstd statistics
                relevant_arr[relevant_arr == 0] = np.nan
                
                #get the stats; plot
                q_mean = np.nanmean(relevant_arr, axis=0)
                q_std = np.nanstd(relevant_arr, axis=0)
    
                destination = os.path.join(self.dirs.analysis_out, sample.name)
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
           
class HDBUMAP(Handler):
    '''
    Performs UMAP embedding of the dataset and clusters the results using 
    HBSCAN clustering method. The object should be pimarily interacted with 
    using its "analysis" method. Full hdbumap analysis involves the following 
    steps:
    
    1. Sample preprocessing. A list of peptide/DNA sequences is converted
       to some numerical representation via featurization. Feature matrix 
       specifies exactly what kind of repr is used for embedding. For
       example, F=None will result in one hot encoding of the sequences, 
       and F='varimax' (works only for peptides) will convert each amino
       acid to an 8D vector encoding amino acid's biophysical parameters.
       
    2. UMAP and HDBSCAN hyperparameter inference. Hyperparameters for both
       are a function of the sample size and composition.
       
    3. UMAP embedding (done with n_components=2) to get 2D coordinates of
       sample's sequences in the embedded space.
       
    4. HDBSCAN clustering of sequences in the embedded space. Several
       clustering hyperparameter configurations will be tested and the
       one leading to the most robust clustering outcome will be chosen
       to compute the final "labels" array.
       
    5. Report the results: plot embeddings, clustering outcomes, etc.
    
    See more details for self.analysis
'''

    def __init__(self, *args):
        super(HDBUMAP, self).__init__(*args)
        self._validate_imports()
        return

    def _validate_imports(self):
        try:
            import umap
            import hdbscan
            from sklearn.preprocessing import StandardScaler  
            
            self.umap = umap
            self.hdbscan = hdbscan
            self.scaler = StandardScaler()
        except:
            msg = 'Failed to import the libraries necessary for HDBSCAN/UMAP analysis. . .'
            self.logger.error(msg)
            raise ImportError(msg)
        return

    def __repr__(self):
        return '<HDBSCAN UMAP analysis object>'    
    
    def _preprocess_arr(self, arr, top_n=None, alphabet=None, F=None):
        
        from clibas.datapreprocessors import DataPreprocessor
        from clibas.misc import sorted_count, arr_purity
        
        #this casting to AnalysisSample is done to ensure that a 2D array
        #in a proper format is passed downstream
        arr = AnalysisSample(X=arr)        
        X, C = sorted_count(arr.X, top_n=top_n)
        p_arr = self._recast_data(X)
        
        #featurize the array; everything is done in memory, so the array has to fit
        pre = DataPreprocessor(dict(F=F))
        p_arr = pre.int_repr(alphabet=alphabet)(p_arr)
        p_arr = pre.featurize_X(alphabet=alphabet,
                                reshape=False,
                                repad=False)(p_arr)

        #setup umap/clustering hyperparams right away to save computation
        purity = arr_purity(arr.X, alphabet)
        arr_size = X.shape[0]
        self._infer_hdbscan_hyperparameters(arr_size, purity)
        self._infer_umap_hyperparameters(arr_size, purity)
            
        return X, C, p_arr[0]['X']
    
    def _initialize_umap(self):
        
        if not hasattr(self, 'n_neighbors'):
            self.n_neighbors = 15
        
        if not hasattr(self, 'min_dist'):
            self.min_dist = 0.1
        
        self.embedder = self.umap.UMAP(n_neighbors=self.n_neighbors,
                                       min_dist=self.min_dist,
                                       n_components=2,
                                       spread=1.0
                                      )
        return
        
    def _infer_umap_hyperparameters(self, size, purity):
        '''
        Two UMAP hyperparameters are optimized: min_dist and n_neighbors
        
        min_dist is a function of array's purity; datasets with low purity need 
        to be repsented on a more "local" manifold is the gist.
        
        n_neighbors depends both on array's purity and on the total number
        of sequences in array
           
        Parameters:
                    size: number of sequences in the array to be embedded
                  purity: array's purity
        '''
        #4PL params are heuristic
        from clibas.misc import logistic_4_param
        self.min_dist = logistic_4_param(purity, 0.01, 0.65, 20, 0.25)
                
        if size < 200:
            self.n_neighbors = 5
        
        else:
            n_neighbors_min = max(0.001 * size, 5)
            n_neighbors_max = min(0.04 * size, 95)
            self.n_neighbors = np.ceil(logistic_4_param(purity,
                                                        n_neighbors_min,
                                                        n_neighbors_max,
                                                        15, 0.35
                                                       )
                                      ).astype(np.int)
        return
     
    def _infer_hdbscan_hyperparameters(self, size, purity):
        '''
        Two HDBSCAN hyperparameters are optimized: min_cluster_size and
                                                   min_samples
        
        Both are the function of array's size and purity. Rather than coming
        up with specific values, both hyperparams are held as ndarrays containing
        up to 8 values each; during clustering a hyperparameter grid will be
        constructed using these arrays and the best configuration will be
        utilized to perform the actual clustering.
           
        Parameters:
                    size: number of sequences in the array to be embedded
                  purity: array's purity
        '''
         
        def hbscan_param_curve(heuristic, size, purity):
            return np.ceil(np.divide(heuristic * np.log2(size), 
                                     -np.log10(purity))
                          )
          
        #for numerical stability
        purity += 0.001
           
        #HBDSCAN params also need to be approximated using array purity as a proxy
        min_cluster_heur = np.logspace(-4, 4, num=8, base=1.4)
        min_samples_heur= np.logspace(-10, 0, num=8, base=1.4)
        
        min_cluster_size = hbscan_param_curve(min_cluster_heur, size, purity)
        min_cluster_size[min_cluster_size < 2] = 2
        self.min_cluster_size = np.unique(min_cluster_size).astype(np.int)
        
        min_samples = hbscan_param_curve(min_samples_heur, size, purity)
        min_samples[min_samples < 1] = 1
        self.min_samples = np.unique(min_samples).astype(np.int)
        return

    def _cluster(self, X, Y, alphabet=None):
        '''
        Optimize and perform HDBSCAN clustering
                                                   
        Parameters:
                    X: array of sequences to be clustered (original str repr)
                    
                    Y: UMAP embeddings of sequences in X; shape=(X.shape[0], 2)
                    
             alphabet: token alphabet for the datasets to be analyzed
                       type: (tuple, list, ndarray)
                     
           Returns:         
               
               labels: ndarray (size=X.shape[0]) of cluster labels for sequences
                       in X
                       
           cluster_hp: summary of clustering hyperparameter optimization;
                       2D ndarray of shape=(min_cluster_size.size; 
                                            min_samples.size)
                       
                       holding overall clustering scores for every combination
                       of hyperparameters.
        '''
                 
        #Y: scaled umap embedding
        if not hasattr(self, 'min_cluster_size'):
            self.min_cluster_size = np.array([5])
            #as array to make it iterable
        
        if not hasattr(self, 'min_samples'):
            self.min_samples = np.array([3])
            
        from clibas.misc import naive_clustering_score
        cluster_hps = np.zeros((self.min_cluster_size.size,
                                self.min_samples.size
                              ))
    
        #scout for the best clustering hyperparameters
        for i,c in enumerate(self.min_cluster_size):
            for j,s in enumerate(self.min_samples):
                
                hdb = self.hdbscan.HDBSCAN(
                                           cluster_selection_epsilon=0.01,
                                           min_cluster_size=int(c),
                                           min_samples=int(s),
                                           cluster_selection_method = 'eom'
                                          )
                
                hdb.fit(Y)
                cluster_hps[i, j] = naive_clustering_score(X,
                                                              hdb.labels_,
                                                              alphabet=alphabet,
                                                              return_mean=True
                                                             ) 
        
        #use the best hyperparamaters to return the actual labels
        i, j = np.where(cluster_hps == cluster_hps.max())
        hdb = self.hdbscan.HDBSCAN(
                                   cluster_selection_epsilon=0.01,
                                   min_cluster_size=int(self.min_cluster_size[i[0]]),
                                   min_samples=int(self.min_samples[j[0]]),
                                   cluster_selection_method = 'eom'
                                  )

        hdb.fit(Y)
        return hdb.labels_, cluster_hps
    
    def embed_and_cluster(self, arr, top_n=None, alphabet=None, F=None):

        X, C, pX = self._preprocess_arr(arr, top_n=top_n, alphabet=alphabet, F=F)
        self._initialize_umap()
        
        Y = self.embedder.fit_transform(pX)
        Y = self.scaler.fit_transform(Y)
        labels, cluster_hps = self._cluster(X, Y, alphabet=alphabet)
        return (X, Y, C, labels, cluster_hps)

    def transform_and_cluster(self, arr, top_n=None, alphabet=None, F=None):

        X, C, pX = self._preprocess_arr(arr, top_n=top_n, alphabet=alphabet, F=F)
        Y = self.embedder.transform(pX)
        Y = self.scaler.transform(Y)
        labels, cluster_scores = self._cluster(X, Y, alphabet=alphabet)
        return (X, Y, C, labels, cluster_scores)

    def cluster_summary(self, X, labels, alphabet=None): 
        #labels should already be +1'd
        from clibas.misc import naive_clustering_score, arr_purity
        cluster_scores = naive_clustering_score(X, 
                                                labels-1,
                                                alphabet=alphabet, 
                                                return_mean=False
                                               )

        cluster_sizes = [X[labels == i].shape[0] 
                         for i in np.unique(labels)
                        ]
        
        #inefficient computation, but i can't imagine this ever being the bottleneck
        cluster_purities = [arr_purity(X[labels==i], alphabet=alphabet) 
                            for i in np.unique(labels)
                           ]
    
        d = {'Cluster purity': cluster_purities,
             'Cluster size': cluster_sizes,
             'Cluster number': np.unique(labels),
             'Cluster score': cluster_scores,
            }
        
        df = pd.DataFrame.from_dict(d)
        df = df.sort_values(['Cluster score'], ascending=False)    
        return df

    def entry_count_cluster(self, X, C, labels):
        '''
        Make a pandas dataframe with three columns: entry, count, cluster
        
        Parameters:
                    X: array of sequences (original str repr) (ndim=2)
                    
                    C: an array of counts for each sequence in X (ndim=1)
                    
               labels: an array of cluster assignments for each sequence in X 
                       (ndim=1)
                     
           Returns:         
               
                  df: pandas dataframe   
        '''
        
        d = {'Entry': np.array([''.join(x) for x in X]),
             'Count': C,
             'Cluster': labels}         
        
        df = pd.DataFrame.from_dict(d)
        df = df.sort_values(['Cluster', 'Count'], ascending=False)
        return df

    def analysis(self,
                 top_n=None,
                 where=None,
                 F=None, 
                 cluster_fasta=False,
                 alphabet=None,
                 return_modified=False,
                 single_manifold=False):

        '''
        Perform UMAP embedding of the dataset and clusters the results using 
        HBSCAN clustering method. Involves the following 
        steps:
        
        1. Sample preprocessing. A list of peptide/DNA sequences is converted
           to some numerical representation via featurization. Feature matrix 
           specifies exactly what kind of repr is used for embedding. For
           example, F=None will result in one hot encoding of the sequences, 
           and F='varimax' (works only for peptides) will convert each amino
           acid to an 8D vector encoding amino acid's biophysical parameters.
           
        2. UMAP and HDBSCAN hyperparameter inference. Hyperparameters for both
           are a function of the sample size and composition.
           
        3. UMAP embedding (done with n_components=2) to get 2D coordinates of
           sample's sequences in the embedded space.
           
        4. HDBSCAN clustering of sequences in the embedded space. Several
           clustering hyperparameter configurations will be tested and the
           one leading to the most robust clustering outcome will be chosen
           to compute the final "labels" array.
           
        5. Report the results: plot embeddings, clustering outcomes, etc.
        
        Parameters:
            
                top_n: int or None (default: None); top N sequences (by count) 
                       to be analyzed. if None, the entire dataset will be 
                       embedded and clustered. Note that if the dataset is large
                       this can take a long time. Recommended: 100 - 5000
                    
               where:  str or None (default: None) 'dna' or 'pep' to specify 
                       which dataset the op should work on. Only set if the op
                       is used as part of pipeline.
                    
             alphabet: tuple, list, ndarray or None; default: None
                       token alphabet for the datasets to be analyzed. if None,
                       the method will attempt to infer the alphabet based on
                       the 'where' keyword. If both are None, an error is raised.
                       
                   F:  str, ndarray or None; default: None
                       Feature matrix for data preprocessing. If None, one-hot
                       encodings will be utilized. A 2D ndarray can be passed
                       to specify the matrix explicitly. In this case, 
                       F.shape[0] = len(alphabet) is the necessary requirement.
                       
                       Several str values are supported: for example, 'pep_ECFP3',
                       'pep_ECFP4', 'pep_SMILES', 'varimax'. See feature matrix
                       documentation for more information 
                       (clibas.featurization.FeatureMatrix)
                       
       cluster_fasta:  bool; default: False. If set, fasta files for individual
                       clusters will be created.
                
     return_modified:  bool; default: False. If set, the op will return
                       modified data (includes top N entries for each sample,
                       UMAP embeddings, clustering labels, etc). If False, the 
                       op returns original data.
    
     single_manifold:  bool; default: False. By default, for data consisting 
                       of several samples, independent embeddings are performed.
                       Because UMAP is stochastic, two embeddings cannot be 
                       directly compared to each other. If set True, all samples
                       in data will be embedded to a single common manifold.

           Returns:         
               
                Data:  data as Data object
                
                       Additionally, several files will be created: for each
                       analyzed sample, umap hdbscan dashboard (an html file)
                       will be created to enable interactive expolation of the 
                       data. If single_manifold is True, an additional .html
                       file showing different sample embeddings side by side
                       will also be created. 
                            
                       Clustering outcomes will be summarized separately in 
                       two .csv files, and a .csv/.png figure will be created 
                       to show the results of clustering optimization.
    '''
         
        if top_n is not None:
            if not isinstance(top_n, int):
                msg = f'<HDBUMAP.analyze> op expected to receive parameter top_n as int; received: {type(top_n)}'
                self.logger.error(msg)  
                raise ValueError(msg)
        else:
            msg = '<HDBUMAP.analyze>: top_n param is unspecified. The op may take a long time to finish if the datasets are large. . .'
            self.logger.warning(msg)
   
        alphabet = self._infer_alphabet(where, alphabet)
        
        from clibas.featurization import FeatureMatrix
        F = FeatureMatrix.make(descr=F, constants=self.constants).F
        
        def _cluster_writer(df, destination, sname):
            
            clusters = np.unique(df['Cluster'])
            for cluster in clusters:
            
                fname = f'{sname}_{where}_cluster_{cluster}.fasta'
                path = os.path.join(destination, fname)    
                with open(path, 'w') as f:
                    for entry in df[df['Cluster'] == cluster].iterrows():
                        f.write(f'>cluster_{cluster}_count_{entry[1]["Count"]}\n')  
                        f.write(f'{entry[1]["Entry"]}\n')
            return           
        
        def umap_hdbscan_summary(data):
                            
            #make sure any input is ok -> recast everything to Data(dtype=AnalysisSample)
            recast_data = self._recast_data(data, where=where)
            func = self.embed_and_cluster
            mod_results = []

            if single_manifold:
                stacked, recast_data = self._restack_and_repad(recast_data)
                
                _, _, pX = self._preprocess_arr(stacked, 
                                                top_n=top_n,
                                                alphabet=alphabet,
                                                F=F
                                               )
                
                self._initialize_umap()
                Y = self.embedder.fit_transform(pX)
                self.scaler.fit(Y)
                func = self.transform_and_cluster
                
            #for nicer formatting downstream
            where_str = where
            if not where_str:
                where_str = ''
            
            self._prepare_destinations(recast_data, self.dirs.analysis_out)            
            for sample in recast_data:
            
                X, Y, C, labels, cluster_hps = func(sample.X,
                                                    top_n=top_n, 
                                                    alphabet=alphabet,
                                                    F=F
                                                   )                
                labels += 1
                destination = os.path.join(self.dirs.analysis_out, sample.name)
                
                #save cluster summary
                cluster_summary = self.cluster_summary(X, labels, alphabet=alphabet)
                fname = f'{sample.name}_{where_str}_clustering_summary.csv'
                full_name = os.path.join(destination, fname)
                cluster_summary.to_csv(full_name, sep=',', index=False)  
                
                sample_dict = {'X': X,
                               'Y': Y,
                               'C': C, 
                               'labels': labels, 
                               'cluster_summary': cluster_summary, 
                               'name': sample.name
                              }
                
                mod_results.append(sample_dict)
            
                #make + save an entry-count-label summary
                df = self.entry_count_cluster(X, C, labels)
                fname = f'{sample.name}_{where_str}_entry_count_cluster.csv'
                full_name = os.path.join(destination, fname)
                df.to_csv(full_name, sep=',', index=False)                     
                
                #plot clustering hyperparams
                basename = os.path.join(destination, 
                                        f'{sample.name}_{where_str}_clustering_optimization')                
                
                Plotter.Analysis.ClusteringHyperParams(self.min_cluster_size,
                                                       self.min_samples,
                                                       cluster_hps,
                                                       sample_name=sample.name,
                                                       basename=basename
                                                       )                
                
                #plot the results
                #maybe worth it to have a matplotlib fallback if holoviews cannot be imported
                from clibas.hv_plotters import hdbumap_analysis_dashboard
                fname = os.path.join(destination, 
                                        f'{sample.name}_{where_str}_umap_hdbscan_dashboard')
                                
                hdbumap_analysis_dashboard(sample_dict, 
                                           alphabet=alphabet,
                                           fname=fname
                                          )
                if cluster_fasta:
                    _cluster_writer(df, destination, sample.name)
                           
            if single_manifold:
                
               from clibas.hv_plotters import single_manifold_embedding_dashboard
               fname = os.path.join(self.dirs.analysis_out, 'single_manifold_embeddings')
               single_manifold_embedding_dashboard(mod_results, fname)
            
            if return_modified:

                for sample_dict in mod_results:
                    del sample_dict['cluster_summary']
                
                data = Data([
                             AnalysisSample(**sample_dict)       
                             for sample_dict in mod_results
                            ])
            
            return data
        return umap_hdbscan_summary
          
        
        
        
