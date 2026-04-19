"""
Data analysis tools for sequencing datasets.

Provides analysis methods for sequence length distributions, quality scores,
library convergence metrics, and dimensionality reduction with clustering.
"""

import inspect
import os

import numpy as np
import pandas as pd

import clibas.plotters as Plotter
from clibas.baseclasses import Handler
from clibas.datatypes import AnalysisSample, Data


class DataAnalysisTools(Handler):
    """
    Collection of analysis tools for sequencing data visualization and statistics.

    Most methods are pipeline operations: they return a callable for use in
    processing pipelines.

    Provides methods for analyzing sequence length distributions, quality scores,
    library convergence metrics, and dimensionality reduction with clustering.
    Typically accessed through the clibas facade after initialization.

    Example:
        >>> import clibas as C
        >>> C.initialize_from_config('config.yaml')
        >>> #analysis tools are now ready to use
        >>> #as C.analysis_tools

    Note:
        This class is not typically instantiated directly. Use the clibas
        initialization system to access analysis functionality.
    """

    def __init__(self, *args):
        super(DataAnalysisTools, self).__init__(*args)

        self._validate_designs()
        self._validate_constants()
        return

    def __repr__(self):
        return "<DataAnalysisTools object>"

    def _recast_for_analysis(self, maybe_arr, where=None):
        if not isinstance(maybe_arr, Data):
            if isinstance(maybe_arr, AnalysisSample):
                return Data([maybe_arr])
            else:
                return Data([AnalysisSample(X=maybe_arr)])

        else:
            if maybe_arr.dtype != AnalysisSample:
                try:
                    return Data(
                        [
                            AnalysisSample(X=sample[where], name=sample.name)
                            for sample in maybe_arr
                        ]
                    )
                except:
                    msg = f"{self} could not recast {where} samples in data. . ."
                    self.logger.error(msg)
                    raise ValueError(msg)
        return maybe_arr

    def length_analysis(self, where=None, save_txt=False):
        """
        Analyze sequence length distributions.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Computes and plots the distribution of sequence lengths for each sample.
        Results are saved as histogram plots in the analysis output directory.

        Args:
            where (str or None): Dataset selector when analyzing
                `SequencingSample` datasets. In such cases, must be either
                'dna' or 'pep'. For direct inputs such as lists, NumPy arrays, or
                pandas Series objects, leave as ``None``. In these cases, the provided
                input is analyzed directly.

            save_txt (bool): If True, save length distribution data to CSV file
                alongside the plots. Default is False.

        Returns:
            callable: Operation that accepts a Data object, generates length distribution analysis, and returns the unmodified Data object.

        Example:
            >>> #analyse the length distribution of peptide datasets
            >>> len_analysis = C.analysis_tools.length_analysis(where='pep', save_txt=True)
            >>> data = len_analysis(data)
        """
        if where is not None:
            self._where_check(where)
            
        la = LengthAnalysis(self.__dict__)
        op = la.len_summary(where=where, save_txt=save_txt)
        return op

    def q_score_analysis(self, loc=None, save_txt=False):
        """
        Analyze quality score distributions.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Computes mean and standard deviation of quality scores for each position
        in the specified regions or across entire sequences. Results are plotted
        in the analysis output directory.

        Args:
            loc (list, optional): List of integers specifying which regions of DNA
                to analyze. If ``None``, analyzes entire sequences. When specified,
                collapses sample internal state.

            save_txt (bool): If True, save quality score statistics to CSV file
                alongside the plots. Default is False.

        Returns:
            callable: Operation that accepts a Data object, generates quality score analysis, and returns the unmodified Data object.

        Example:
            >>> #summarize Q scores in regions 0 and 1
            >>> q_analysis = C.analysis_tools.q_score_analysis(loc=[0, 1], save_txt=True)
            >>> data = q_analysis(data)
        """
        if loc is not None:
            self._infer_design("dna")
            self._loc_check(loc, self.D_design)

        QA = QScoreAnalysis(self.__dict__)
        op = QA.q_summary(loc=loc, save_txt=save_txt)
        return op

    def sequence_convergence_analysis(self, where=None):
        """
        Analyze library convergence at the sequence level.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Computes normalized Shannon entropy and position-wise sequence conservation
        to assess library convergence. Results are plotted in the analysis output
        directory.

        Args:
            where (str or None): Dataset selector when analyzing
                `SequencingSample` datasets. In such cases, must be either
                'dna' or 'pep'. For direct inputs such as lists, NumPy arrays, or
                pandas Series objects, leave as ``None``. In these cases, the provided
                input is analyzed directly.

        Returns:
            callable: Operation that accepts a Data object, generates convergence analysis, and returns the unmodified Data object.

        Example:
            >>> conv_analysis = C.analysis_tools.sequence_convergence_analysis(where='pep')
            >>> data = conv_analysis(data)
        """
        if where is not None:
            self._where_check(where)
            
        CA = ConvergenceAnalys(self.__dict__)
        op = CA.sequence_level_convergence(where=where)
        return op

    def token_convergence_analysis(
        self, where=None, loc=None, alphabet=None, save_txt=False
    ):
        """
        Analyze library convergence at the token (amino acid/nucleobase) level.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Computes the frequency of each token at each position in the dataset.
        Generates frequency heatmaps and sequence logos. When analyzing specific
        regions, collapses sample internal state.

        Args:
            where (str or None): Dataset selector when analyzing
                `SequencingSample` datasets. In such cases, must be either
                'dna' or 'pep'. For direct inputs such as lists, NumPy arrays, or
                pandas Series objects, leave as ``None``. In these cases, the provided
                input is analyzed directly. If left as ``None``, the ``alphabet``
                argument must be specified.

            loc (list, optional): List of integers specifying which regions to
                analyze. If ``None``, analyzes entire sequences. When specified,
                collapses sample internal state.

            alphabet (list, tuple, ndarray, or str, optional): Token alphabet for
                analysis. Can be one of:
                
                - a sequence (list, tuple, or ndarray) specifying custom tokens, or
                - the string 'aa' for the amino acid alphabet as specified in the config, or
                - the string 'base' for the nucleotide base alphabet as specified in the config.
                
                If ``None``, the alphabet is automatically inferred from 'where'.

            save_txt (bool): If True, save frequency data to CSV file alongside
                the plots. Default is False.

        Returns:
            callable: Operation that accepts a Data object, generates token-level convergence analysis, and returns the unmodified Data object.

        Example:
            >>> token_analysis = C.analysis_tools.token_convergence_analysis(
            ...     where='pep', loc=[1, 3], save_txt=True
            ... )
            >>> data = token_analysis(data)
        """
        if where is not None:
            self._where_check(where)
            
        alphabet = self._infer_alphabet(where, alphabet)
        
        if loc is not None:
            design = self._infer_design(where)
            self._loc_check(loc, design)

        CA = ConvergenceAnalys(self.__dict__)
        op = CA.token_level_convergence(
            where=where, loc=loc, alphabet=alphabet, save_txt=save_txt
        )
        return op

    def umap_hdbscan_analysis(
        self,
        top_n=None,
        where=None,
        F=None,
        cluster_fasta=False,
        alphabet=None,
        return_modified=False,
        single_manifold=False,
    ):
        """
        Perform UMAP dimensionality reduction and HDBSCAN clustering.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Embeds sequences into 2D space using UMAP and clusters them using HDBSCAN.
        Automatically optimizes hyperparameters based on dataset size and composition.
        Generates interactive HTML dashboards and clustering summaries.

        The analysis pipeline:
        1. Featurize sequences (one-hot encoding or custom feature matrix)
        2. Infer optimal UMAP and HDBSCAN hyperparameters
        3. Embed sequences to 2D using UMAP
        4. Cluster embedded sequences with HDBSCAN
        5. Generate interactive visualizations and summaries

        Args:
            top_n (int, optional): Number of most abundant sequences to analyze.
                If None, analyzes entire dataset. Recommended: 100-5000 for
                computational efficiency.

            where (str or None): Dataset selector when analyzing
                `SequencingSample` datasets. In such cases, must be either
                'dna' or 'pep'. For direct inputs such as lists, NumPy arrays, or
                pandas Series objects, leave as ``None``. In these cases, the provided
                input is analyzed directly. If left as None, the ``alphabet``
                argument must be specified.

            F (str, ndarray, or None): Feature matrix specification. If ``None``,
                uses one-hot encoding. String options include 'varimax',
                'pep_ECFP3', 'pep_ECFP4', 'pep_SMILES'. Can also provide
                custom 2D array with ``F.shape[0] == len(alphabet)``.

            cluster_fasta (bool): If True, generates separate FASTA files for
                each cluster. Default is False.

            alphabet (list, tuple, ndarray, or str, optional): Token alphabet for
                analysis. Can be one of:
                
                - a sequence (list, tuple, or ndarray) specifying custom symbols, or
                - the string 'aa' for the amino acid alphabet as specified in the config, or
                - the string 'base' for the nucleotide base alphabet as specified in the config.
                
                If ``None``, the alphabet is automatically inferred from 'where'.

            return_modified (bool): If True, returns Data object with embeddings
                and cluster labels. If False, returns unmodified original data.
                Default is False.

            single_manifold (bool): If True, embeds all samples to a single
                shared manifold for direct comparison. If False, each sample
                gets independent embedding. Default is False.

        Returns:
            callable: Operation that accepts a Data object, performs UMAP/HDBSCAN analysis, and returns either the original or modified Data object.

        Note:
            Generates multiple output files per sample:
              - Interactive HTML dashboard for exploration
              - Clustering summary CSV
              - Entry-count-cluster mapping CSV
              - Static matplotlib plots
              - Optional: per-cluster FASTA files
              - Optional: single manifold comparison dashboard (if single_manifold=True)

        Example:
            >>> umap_analysis = C.analysis_tools.umap_hdbscan_analysis(
            ...     top_n=1000, where='pep', F='pep_ECFP4', cluster_fasta=True
            ... )
            >>> data = umap_analysis(data)
            >>>
            >>> #for comparing multiple samples on same manifold
            >>> umap_analysis = C.analysis_tools.umap_hdbscan_analysis(
            ...     top_n=500, where='pep', single_manifold=True
            ... )
            >>> data = umap_analysis(data)
        """
        # for historic reasons, here we defer to the child class to 
        # do the validation of arguments
        hdu = HDBUMAP(self.__dict__)
        op = hdu.analysis(
            top_n=top_n,
            where=where,
            F=F,
            cluster_fasta=cluster_fasta,
            alphabet=alphabet,
            return_modified=return_modified,
            single_manifold=single_manifold,
        )
        return op


class LengthAnalysis(DataAnalysisTools):
    """
    Sequence length distribution analysis.

    Internal class for length analysis operations. Not invoked directly.
    """

    def __init__(self, *args):
        super(LengthAnalysis, self).__init__(*args)

    def __repr__(self):
        return "<LengthAnalysis object>"

    def len_summary(self, where=None, save_txt=False):
        # Docstring already covered by DataAnalysisTools.length_analysis
        def length_summary(data):
            recast_data = self._recast_for_analysis(data, where=where)
            self._prepare_destinations(root=self.dirs.analysis_out, data=recast_data)
            for sample in recast_data:
                arr = sample.X
                self._empty_array_check(arr, inspect.stack()[0][3])

                L = self._L_summary(arr)
                L, counts = np.unique(L, return_counts=True)

                destination = os.path.join(self.dirs.analysis_out, sample.name)
                fname = f"{sample.name}_{where}_L_distribution" + self.timestamp
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.L_distribution(
                    L, counts, where=where, basename=basename
                )

                if save_txt:
                    np.savetxt(
                        basename + ".csv",
                        np.array((L, counts)).T,
                        delimiter=",",
                        header="Seq length,Count",
                    )

            return data

        return length_summary


class ConvergenceAnalys(DataAnalysisTools):
    """
    Library convergence analysis tools.

    Internal class for convergence analysis operations. Not invoked directly.
    """

    def __init__(self, *args):
        super(ConvergenceAnalys, self).__init__(*args)

        from clibas.misc import get_freqs, shannon_entropy

        self.shannon_entropy = shannon_entropy
        self.get_freqs = get_freqs
        return

    def __repr__(self):
        return "<ConvergenceAnalys object>"

    def sequence_level_convergence(self, where=None):
        # docstring already covered by DataAnalysisTools.sequence_convergence_analysis
        def sequence_level_convergence_summary(data):
            recast_data = self._recast_for_analysis(data, where)
            self._prepare_destinations(root=self.dirs.analysis_out, data=recast_data)

            for sample in recast_data:
                arr = sample.X
                self._empty_array_check(arr, inspect.stack()[0][3])

                shannon, counts = self.shannon_entropy(arr, norm=True)

                destination = os.path.join(self.dirs.analysis_out, sample.name)
                fname = f"{sample.name}_{where}_library_convergence" + self.timestamp
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.dataset_convergence(
                    counts, shannon, where, basename
                )

            return data

        return sequence_level_convergence_summary

    def token_level_convergence(
        self, where=None, loc=None, alphabet=None, save_txt=False
    ):
        # docstring already covered by DataAnalysisTools.token_convergence_analysis
        def token_level_convergence_analysis(data):
            if loc is None:
                recast_data = self._recast_for_analysis(data, where)
            else:
                recast_data = data

            self._prepare_destinations(root=self.dirs.analysis_out, data=recast_data)

            for sample in recast_data:
                if loc is None:
                    arr = sample.X
                    self._empty_array_check(arr, inspect.stack()[0][3])
                    # self._joint_alphabet_X_check(arr, alphabet) #TODO: think how to deal with this
                    # the problem is that a pipeline may feed this op arrays containing junk
                    # token not encoded in the alphabet. But on the other hand, it would make
                    # sense to make that the _joint_alphabet_X_check check passes
                    
                    freq = self.get_freqs(arr, alphabet)
                    nloc = "overall"
                    fname = (
                        f"{sample.name}_{where}_tokenwise_frequency" + self.timestamp
                    )

                else:
                    arr = self._cast(sample[where], "2d")
                    self._empty_array_check(arr, inspect.stack()[0][3])
                    # self._joint_alphabet_X_check(arr, alphabet) #TODO: think how to deal with this
                    design = self._infer_design(where)
                    # array internal state has to be collapsed for this calculation
                    if not sample._is_collapsed():
                        msg = f"<frequency_summary> op will collapse sample {sample.name}'s internal state"
                        self.logger.debug(msg)
                        sample._collapse_internal_state()

                    # initialize the frequency array: 
                    # the 3D array has to be reduced along axis 0 at the end
                    maxlen = self._find_max_len(design, loc)
                    freq = np.zeros(
                        (len(design), len(alphabet), maxlen), dtype=np.float32
                    )

                    for i, template in enumerate(design):
                        row_mask = sample._internal_state[:, i]
                        col_mask = template(loc, return_mask=True)

                        # calculated weighed contributions of each design
                        # to the overall frequency array
                        norm = np.divide(np.sum(row_mask), arr.shape[0])
                        freq[i, :, : len(col_mask)] = norm * np.nan_to_num(
                            self.get_freqs(arr[row_mask][:, col_mask], alphabet)
                        )

                    # reduce back to a 2D array and plot/save
                    freq = np.sum(freq, axis=0)
                    nloc = ", ".join(str(x + 1) for x in loc)

                fname = (
                    f"{sample.name}_{where}_reg_{nloc}_tokenwise_frequency"
                    + self.timestamp
                )
                flogo = f"{sample.name}_{where}_reg_{nloc}_seq_logo" + self.timestamp

                destination = os.path.join(self.dirs.analysis_out, sample.name)
                logobase = os.path.join(destination, flogo)
                basename = os.path.join(destination, fname)

                str_alphabet = [s.decode("ascii") for s in alphabet]
                Plotter.SequencingData.tokenwise_frequency(
                    freq, str_alphabet, where, nloc, basename
                )
                Plotter.Analysis.seq_logo(
                    X=None,
                    C=None,
                    freq=pd.DataFrame(
                        freq.T, index=range(freq.shape[1]), columns=str_alphabet
                    ),
                    alphabet=str_alphabet,
                    basename=logobase,
                )
                if save_txt:
                    np.savetxt(basename + ".csv", freq, delimiter=",")

            return data

        return token_level_convergence_analysis


class QScoreAnalysis(DataAnalysisTools):
    """
    Quality score analysis tools.

    Internal class for quality score analysis operations. Not invoked directly.
    """

    def __init__(self, *args):
        super(QScoreAnalysis, self).__init__(*args)

    def __repr__(self):
        return "<QScoreAnalysis object>"

    def q_summary(self, loc=None, save_txt=False):
        # docstring already covered by DataAnalysisTools.q_score_analysis
        def q_score_summary(data):
            import numpy.ma as ma

            self._prepare_destinations(root=self.dirs.analysis_out, data=data)

            for sample in data:
                if loc is None:
                    relevant_arr = sample.civilized_Q
                    self._empty_array_check(relevant_arr, inspect.stack()[0][3])
                    nloc = "overall"
                    fname = f"{sample.name}_q_score_summary" + self.timestamp

                else:
                    arr = sample.civilized_Q
                    self._empty_array_check(arr, inspect.stack()[0][3])

                    if not sample._is_collapsed():
                        msg = f"<q_score_summary> op will collapse sample {sample.name}'s internal state"
                        self.logger.debug(msg)
                        sample._collapse_internal_state()

                    maxlen = self._find_max_len(self.D_design, loc)
                    # iterate over templates and append all of the relevant arr views to this array
                    # relevant view: masked (row/columnwise) arr
                    relevant_arr = []

                    for i, template in enumerate(self.D_design):
                        row_mask = sample._internal_state[:, i]
                        col_mask = template(loc, return_mask=True)

                        arr_view = np.zeros(
                            (np.sum(row_mask), maxlen), dtype=np.float32
                        )
                        arr_view[:, : len(col_mask)] = arr[row_mask][:, col_mask]
                        relevant_arr.append(arr_view)

                    # assemble into a single array
                    relevant_arr = np.vstack(relevant_arr)

                    nloc = ", ".join(str(x + 1) for x in loc)
                    fname = f"{sample.name}_reg{nloc}_q_score_summary" + self.timestamp

                # mask out pads (0) as nans for nanmean/nanstd statistics
                relevant_arr = ma.masked_equal(relevant_arr, 0)

                # get the stats; plot
                q_mean = ma.mean(relevant_arr, axis=0)
                q_std = ma.std(relevant_arr, axis=0)

                destination = os.path.join(self.dirs.analysis_out, sample.name)
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.Q_score_summary(q_mean, q_std, nloc, basename)

                if save_txt:
                    q = np.vstack((q_mean, q_std))
                    np.savetxt(
                        basename + ".csv", q.T, delimiter=",", header="Q mean, Q std"
                    )
            return data

        return q_score_summary


class HDBUMAP(DataAnalysisTools):
    """
    UMAP dimensionality reduction with HDBSCAN clustering.

    Performs automated embedding and clustering of sequence data with optimized
    hyperparameters. Internal class accessed through DataAnalysisTools.umap_hdbscan_analysis.

    The workflow:
    1. Featurize sequences using specified or default encoding

    2. Infer optimal UMAP hyperparameters (n_neighbors, min_dist) based on
       dataset size and sequence diversity (purity)

    3. Infer optimal HDBSCAN hyperparameters (min_cluster_size, min_samples)

    4. Embed sequences to 2D using UMAP

    5. Test multiple HDBSCAN configurations and select best clustering

    6. Generate visualizations and summaries

    Note:
        Not invoked directly. Use through DataAnalysisTools.umap_hdbscan_analysis.
    """

    def __init__(self, *args):
        super(HDBUMAP, self).__init__(*args)
        self._validate_imports()
        return

    def _validate_imports(self):
        try:
            import hdbscan
            import umap
            from sklearn.preprocessing import StandardScaler

            self.umap = umap
            self.hdbscan = hdbscan
            self.scaler = StandardScaler()

            # as of May 2025, using hdbscan will raise sklearn's future warnings
            # let's silence these for now
            import warnings

            warnings.filterwarnings(
                "ignore", category=FutureWarning, module="sklearn.*"
            )
        except ImportError:
            msg = "Failed to import the packages necessary for HDBSCAN/UMAP analysis. Please install clibas with `pip install clibas[ml]`. . ."
            self.logger.error(msg)
            raise ImportError(msg)
        return

    def __repr__(self):
        return "<HDBSCAN UMAP analysis object>"

    def _preprocess_arr(self, arr, top_n=None, alphabet=None, F=None):
        from clibas.datapreprocessors import DataPreprocessor
        from clibas.misc import arr_purity, sorted_count

        # this casting to AnalysisSample is done to ensure that a 2D array
        # in a proper format is passed downstream
        arr = AnalysisSample(X=arr)
        X, C = sorted_count(arr.X, top_n=top_n)
        p_arr = self._recast_for_analysis(X)

        # featurize the array; everything is done in memory, so the array has to fit
        p_arr = DataPreprocessor(self.__dict__).featurize_X(
            F=F,
            alphabet=alphabet,
            reshape=False,
        )(p_arr)

        # setup umap/clustering hyperparams right away to save computation
        purity = arr_purity(arr.X, alphabet)
        arr_size = X.shape[0]
        self._infer_hdbscan_hyperparameters(arr_size, purity)
        self._infer_umap_hyperparameters(arr_size, purity)

        return X, C, p_arr[0]["X"]

    def _initialize_umap(self):
        if not hasattr(self, "n_neighbors"):
            self.n_neighbors = 15

        if not hasattr(self, "min_dist"):
            self.min_dist = 0.1

        self.embedder = self.umap.UMAP(
            n_neighbors=self.n_neighbors,
            min_dist=self.min_dist,
            n_components=2,
            spread=1.0,
        )
        return

    def _infer_umap_hyperparameters(self, size, purity):
        """
        Infer optimal UMAP hyperparameters based on dataset characteristics.

        Optimizes min_dist and n_neighbors parameters. min_dist is determined by
        array purity (lower purity requires more local manifold representation).
        n_neighbors depends on both purity and total sequence count.

        Args:
            size (int): Number of sequences in the array to be embedded.
            purity (float): Array purity metric (0-1 scale).

        Note:
            Uses heuristic 4-parameter logistic curves for optimization.
        """
        # 4PL params are heuristic
        from clibas.misc import logistic_4_param

        self.min_dist = logistic_4_param(purity, 0.01, 0.65, 20, 0.25)

        if size < 200:
            self.n_neighbors = 5

        else:
            n_neighbors_min = max(0.001 * size, 5)
            n_neighbors_max = min(0.04 * size, 95)
            self.n_neighbors = np.ceil(
                logistic_4_param(purity, n_neighbors_min, n_neighbors_max, 15, 0.35)
            ).astype(int)
        return

    def _infer_hdbscan_hyperparameters(self, size, purity):
        """
        Infer optimal HDBSCAN hyperparameters based on dataset characteristics.

        Optimizes min_cluster_size and min_samples parameters. Both are functions
        of array size and purity. Instead of single values, generates arrays of
        up to 8 candidate values for each parameter. During clustering, a grid
        search over these candidates identifies the optimal configuration.

        Args:
            size (int): Number of sequences in the array to be clustered.
            purity (float): Array purity metric (0-1 scale).

        Note:
            Sets self.min_cluster_size and self.min_samples as arrays for grid search.
        """

        def hbscan_param_curve(heuristic, size, purity):
            return np.ceil(np.divide(heuristic * np.log2(size), -np.log10(purity)))

        # for numerical stability
        purity += 0.001

        # HBDSCAN params also need to be approximated using array purity as a proxy
        min_cluster_heur = np.logspace(-4, 4, num=8, base=1.4)
        min_samples_heur = np.logspace(-10, 0, num=8, base=1.4)

        min_cluster_size = hbscan_param_curve(min_cluster_heur, size, purity)
        min_cluster_size[min_cluster_size < 2] = 2
        self.min_cluster_size = np.unique(min_cluster_size).astype(int)

        min_samples = hbscan_param_curve(min_samples_heur, size, purity)
        min_samples[min_samples < 1] = 1
        self.min_samples = np.unique(min_samples).astype(int)
        return

    def _cluster(self, X, Y, alphabet=None):
        """
        Optimize HDBSCAN hyperparameters and perform clustering.

        Performs grid search over min_cluster_size and min_samples parameters,
        evaluating clustering quality for each combination. Returns cluster labels
        using the best hyperparameter configuration.

        Args:
            X (ndarray): Array of sequences to cluster (original string representation).
            Y (ndarray): UMAP embeddings of sequences in X. Shape: (X.shape[0], 2).
            alphabet (tuple, list, or ndarray, optional): Token alphabet for the dataset.

        Returns:
            tuple: (labels, cluster_hps) where:
                - labels: ndarray of cluster assignments for sequences in X (size=X.shape[0])
                - cluster_hps: 2D array of clustering scores for each hyperparameter
                  combination. Shape: (min_cluster_size.size, min_samples.size)
        """
        # Y: scaled umap embedding
        if not hasattr(self, "min_cluster_size"):
            self.min_cluster_size = np.array([5])
            # as array to make it iterable

        if not hasattr(self, "min_samples"):
            self.min_samples = np.array([3])

        from clibas.misc import naive_clustering_score

        cluster_hps = np.zeros((self.min_cluster_size.size, self.min_samples.size))

        # scout for the best clustering hyperparameters
        for i, c in enumerate(self.min_cluster_size):
            for j, s in enumerate(self.min_samples):
                hdb = self.hdbscan.HDBSCAN(
                    cluster_selection_epsilon=0.01,
                    min_cluster_size=int(c),
                    min_samples=int(s),
                    cluster_selection_method="eom",
                )

                hdb.fit(Y)
                cluster_hps[i, j] = naive_clustering_score(
                    X, hdb.labels_, alphabet=alphabet, return_mean=True
                )

        # use the best hyperparamaters to return the actual labels
        i, j = np.where(cluster_hps == cluster_hps.max())
        hdb = self.hdbscan.HDBSCAN(
            cluster_selection_epsilon=0.01,
            min_cluster_size=int(self.min_cluster_size[i[0]]),
            min_samples=int(self.min_samples[j[0]]),
            cluster_selection_method="eom",
        )

        hdb.fit(Y)
        return hdb.labels_, cluster_hps

    def _embed_and_cluster(self, arr, top_n=None, alphabet=None, F=None):
        """
        Fit UMAP embedding and perform clustering.

        Preprocesses array, fits UMAP embedding, and performs HDBSCAN clustering
        with hyperparameter optimization.

        Args:
            arr: Array of sequences to embed and cluster.
            top_n (int, optional): Number of top sequences to analyze.
            alphabet: Token alphabet for the dataset.
            F: Feature matrix specification.

        Returns:
            tuple: (X, Y, C, labels, cluster_hps) where:
                - X: Original sequences (top N)
                - Y: 2D UMAP embeddings (scaled)
                - C: Counts for each sequence
                - labels: Cluster assignments
                - cluster_hps: Hyperparameter optimization scores
        """
        X, C, pX = self._preprocess_arr(arr, top_n=top_n, alphabet=alphabet, F=F)
        self._initialize_umap()

        Y = self.embedder.fit_transform(pX)
        Y = self.scaler.fit_transform(Y)
        labels, cluster_hps = self._cluster(X, Y, alphabet=alphabet)
        return (X, Y, C, labels, cluster_hps)

    def _transform_and_cluster(self, arr, top_n=None, alphabet=None, F=None):
        """
        Transform using existing UMAP embedding and perform clustering.

        Uses previously fitted UMAP embedder to transform new data, then performs
        HDBSCAN clustering. Used when analyzing multiple samples on single manifold.

        Args:
            arr: Array of sequences to transform and cluster.
            top_n (int, optional): Number of top sequences to analyze.
            alphabet: Token alphabet for the dataset.
            F: Feature matrix specification.

        Returns:
            tuple: (X, Y, C, labels, cluster_scores) where:
                - X: Original sequences (top N)
                - Y: 2D UMAP embeddings (scaled)
                - C: Counts for each sequence
                - labels: Cluster assignments
                - cluster_scores: Clustering quality scores
        """
        X, C, pX = self._preprocess_arr(arr, top_n=top_n, alphabet=alphabet, F=F)
        Y = self.embedder.transform(pX)
        Y = self.scaler.transform(Y)
        labels, cluster_scores = self._cluster(X, Y, alphabet=alphabet)
        return (X, Y, C, labels, cluster_scores)

    def _cluster_summary(self, X, labels, alphabet=None):
        """
        Generate summary statistics for each cluster.

        Computes cluster purity, size, and clustering quality score for each
        identified cluster.

        Args:
            X: Array of sequences.
            labels: Cluster assignments (already incremented by 1).
            alphabet: Token alphabet for the dataset.

        Returns:
            DataFrame: Summary with columns for cluster number, purity, size,
                and clustering score. Sorted by clustering score (descending).
        """
        from clibas.misc import arr_purity, naive_clustering_score

        cluster_scores = naive_clustering_score(
            X, labels - 1, alphabet=alphabet, return_mean=False
        )

        cluster_sizes = [X[labels == i].shape[0] for i in np.unique(labels)]

        # inefficient computation, but i can't imagine this ever being the bottleneck
        cluster_purities = [
            arr_purity(X[labels == i], alphabet=alphabet) for i in np.unique(labels)
        ]

        d = {
            "Cluster purity": cluster_purities,
            "Cluster size": cluster_sizes,
            "Cluster number": np.unique(labels),
            "Cluster score": cluster_scores,
        }

        df = pd.DataFrame.from_dict(d)
        df = df.sort_values(["Cluster score"], ascending=False)
        return df

    def _entry_count_cluster(self, X, C, labels):
        """
        Create DataFrame mapping sequences to counts and clusters.

        Args:
            X: Array of sequences (2D).
            C: Array of counts for each sequence (1D).
            labels: Array of cluster assignments (1D).

        Returns:
            DataFrame: Three columns (Entry, Count, Cluster), sorted by cluster
                and count (descending).
        """

        d = {
            "Entry": np.array([b"".join(x).decode("ascii") for x in X]),
            "Count": C,
            "Cluster": labels,
        }

        df = pd.DataFrame.from_dict(d)
        df = df.sort_values(["Cluster", "Count"], ascending=False)
        return df

    def analysis(
        self,
        top_n=None,
        where=None,
        F=None,
        cluster_fasta=False,
        alphabet=None,
        return_modified=False,
        single_manifold=False,
    ):
        # docstring covered by DataAnalysisTools.umap_hdbscan_analysis
        if top_n is not None:
            if not isinstance(top_n, int):
                msg = f"<HDBUMAP.analyze> op expected to receive parameter top_n as int; received: {type(top_n)}"
                self.logger.error(msg)
                raise ValueError(msg)
        else:
            msg = "<HDBUMAP.analyze>: top_n param is unspecified. The op may take a long time to finish if the datasets are large. . ."
            self.logger.warning(msg)

        if where is not None:
            self._where_check(where)

        alphabet = self._infer_alphabet(where, alphabet)
        from clibas.featurization import FeatureMatrix

        F = FeatureMatrix.make(descr=F, constants=self.constants).F

        def _cluster_writer(df, destination, sname):
            clusters = np.unique(df["Cluster"])
            for cluster in clusters:
                fname = f"{sname}_{where}_cluster_{cluster}.fasta"
                path = os.path.join(destination, fname)
                with open(path, "w") as f:
                    for entry in df[df["Cluster"] == cluster].iterrows():
                        f.write(f">cluster_{cluster}_count_{entry[1]['Count']}\n")
                        f.write(f"{entry[1]['Entry']}\n")
            return

        def umap_hdbscan_summary(data):
            # make sure any input is ok -> recast everything to Data(dtype=AnalysisSample)
            recast_data = self._recast_for_analysis(data, where=where)
            func = self._embed_and_cluster
            mod_results = []

            if single_manifold:
                stacked = recast_data.stack(in_place=False)[0].X

                _, _, pX = self._preprocess_arr(
                    stacked, top_n=top_n, alphabet=alphabet, F=F
                )

                self._initialize_umap()
                Y = self.embedder.fit_transform(pX)
                self.scaler.fit(Y)
                func = self._transform_and_cluster

            # for nicer formatting downstream
            where_str = where
            if not where_str:
                where_str = ""

            self._prepare_destinations(root=self.dirs.analysis_out, data=recast_data)
            for sample in recast_data:
                self._empty_array_check(sample.X, inspect.stack()[0][3])
                self._joint_alphabet_X_check(sample.X, alphabet)
                X, Y, C, labels, cluster_hps = func(
                    sample.X, top_n=top_n, alphabet=alphabet, F=F
                )
                labels += 1
                destination = os.path.join(self.dirs.analysis_out, sample.name)

                # save cluster summary
                cluster_summary = self._cluster_summary(X, labels, alphabet=alphabet)
                fname = (
                    f"{sample.name}_{where_str}_clustering_summary"
                    + self.timestamp
                    + ".csv"
                )
                full_name = os.path.join(destination, fname)
                cluster_summary.to_csv(full_name, sep=",", index=False)

                sample_dict = {
                    "X": X,
                    "Y": Y,
                    "C": C,
                    "labels": labels,
                    "cluster_summary": cluster_summary,
                    "name": sample.name,
                    "total_reads": sample.size,
                }

                mod_results.append(sample_dict)

                # make + save an entry-count-label summary
                df = self._entry_count_cluster(X, C, labels)
                fname = (
                    f"{sample.name}_{where_str}_entry_count_cluster"
                    + self.timestamp
                    + ".csv"
                )
                full_name = os.path.join(destination, fname)
                df.to_csv(full_name, sep=",", index=False)

                # plot clustering hyperparams IF the logger level is low enough
                if self.logger.level == 10:
                    basename = os.path.join(
                        destination,
                        f"{sample.name}_{where_str}_clustering_optimization"
                        + self.timestamp,
                    )

                    Plotter.Analysis.ClusteringHyperParams(
                        self.min_cluster_size,
                        self.min_samples,
                        cluster_hps,
                        sample_name=sample.name,
                        basename=basename,
                    )
                # plot the results
                from clibas.px_plotters import four_panel_dashboard

                fname = os.path.join(
                    destination,
                    f"{sample.name}_{where_str}_umap_hdbscan_dashboard"
                    + self.timestamp,
                )

                str_alphabet = [s.decode("ascii") for s in alphabet]
                four_panel_dashboard(sample_dict, alphabet=str_alphabet, fname=fname)

                # matplotlib version just in case...
                fname = os.path.join(
                    destination,
                    f"{sample.name}_{where_str}_umap_hdbscan_mpl" + self.timestamp,
                )

                Plotter.Analysis.UMAP_HDBSCAN(
                    Y,
                    labels,
                    C=C,
                    sample_name=sample.name,
                    total_reads=sample.size,
                    basename=fname,
                    show_annotations=False,
                )
                if cluster_fasta:
                    _cluster_writer(df, destination, sample.name)

            if single_manifold:
                from clibas.px_plotters import single_manifold_embedding_dashboard

                fname = self.exp_name + "_single_manifold_embeddings" + self.timestamp
                fname = os.path.join(self.dirs.analysis_out, fname)
                single_manifold_embedding_dashboard(mod_results, fname)

            if return_modified:
                # delete cluster summary because it won't fit as a sample dataset
                for sample_dict in mod_results:
                    del sample_dict["cluster_summary"]

                # return everything else
                data = Data(
                    [AnalysisSample(**sample_dict) for sample_dict in mod_results]
                )
            return data

        return umap_hdbscan_summary
