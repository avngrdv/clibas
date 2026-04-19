"""
Data preprocessing utilities for machine learning workflows.

Provides DataPreprocessor onject for preparing sequencing data for machine learning,
including filtering, featurization, train/test splitting, and data augmentation.
"""

import copy
import inspect
import os
from collections import defaultdict

import numpy as np

import clibas.featurization as featurize
from clibas.baseclasses import Handler
from clibas.datatypes import AnalysisSample, Data


class DataPreprocessor(Handler):
    """
    Preprocessing tools for machine learning dataset preparation.

    Provides methods for filtering, featurization, sampling, and data augmentation
    of sequence datasets. Operations transform Data objects for downstream machine
    learning workflows. Typically accessed through the clibas facade after initialization.

    Example:
        >>> import clibas as C
        >>> C.initialize_from_config('config.yaml')
        >>> # Preprocessor is now ready to use
        >>> # as C.preprocessor

    Note:
        This class is not typically instantiated directly. Use the clibas
        initialization system to access preprocessing functionality.
    """

    def __init__(self, *args):
        super(DataPreprocessor, self).__init__(*args)
        return

    def __repr__(self):
        return "<DataPreprocessor object>"


    def token_filter(self, tokens_to_filter_by=None):
        """
        Filter sequences containing specific tokens.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Removes all sequences (X array entries) that contain any of the specified
        tokens. Useful for filtering out sequences with ambiguous or unwanted
        characters.

        Args:
            tokens_to_filter_by (list, tuple, or ndarray): Single-letter encoded
                tokens to filter out (e.g., amino acids, bases).

        Returns:
            callable: Operation that accepts a Data object, filters sequences containing specified tokens, and returns the modified Data object.

        Example:
            >>> #remove sequences containing 'X' or 'Z' amino acids
            >>> token_filt = C.preprocessor.token_filter(tokens_to_filter_by=['X', 'Z'])
            >>> data = token_filt(data)
        """
        if tokens_to_filter_by is None:
            msg = "<filter_by_token> op expected aas_to_filter_by argument."
            self.logger.error(msg)
            raise ValueError(msg)

        if not isinstance(tokens_to_filter_by, (tuple, list, np.ndarray)):
            msg = "<filter_by_token> op expected param tokens_to_filter_by as type=(list, tuple, np.ndarray); received: {type(tokens_to_filter_by)}"
            self.logger.error(msg)
            raise ValueError(msg)

        def filter_by_token(data):
            for sample in data:
                arr = self._cast(sample.X, "2d")
                self._empty_array_check(arr, inspect.stack()[0][3])

                t = np.asarray(tokens_to_filter_by).astype(arr.dtype)

                ind = ~np.isin(arr, t).reshape(arr.shape)
                ind = np.sum(ind, axis=1) == arr.shape[1]

                sample.ind_filter(ind)

            return data

        return filter_by_token

    def intrasample_unique(self):
        """
        Remove duplicate sequences within each sample.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Removes duplicate entries within each sample's X dataset. Equivalent to
        calling np.unique(X, axis=0) on each sample. Entries are resorted during
        the process.

        Returns:
            callable: Operation that accepts a Data object, removes intra-sample duplicates, and returns the modified Data object.

        Example:
            >>> unique_op = C.preprocessor.intrasample_unique()
            >>> data = unique_op(data)
        """

        def get_intrasample_unique(data):
            for sample in data:
                self._empty_array_check(sample.X, inspect.stack()[0][3])

                ind = np.unique(sample.X, axis=0, return_index=True)[1]
                sample.ind_filter(ind)

            return data

        return get_intrasample_unique

    def intersample_unique(self):
        """
        Remove sequences found in multiple samples.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Compares X arrays across all samples and removes entries that appear in
        more than one sample. Only sequences unique to a single sample are retained.

        Returns:
            callable: Operation that accepts a Data object, removes inter-sample duplicates, and returns the modified Data object.

        Note:
            If X arrays have different widths, they will be padded to the maximum
            width before comparison.

        Example:
            >>> intersample_op = C.preprocessor.intersample_unique()
            >>> data = intersample_op(data)
        """

        def _pad_to_new_dim(arr, new_dim):
            to_pad = new_dim - arr.shape[1]
            if to_pad == 0:
                return arr

            new_arr = np.zeros((arr.shape[0], new_dim), dtype=arr.dtype)
            new_arr[:, : arr.shape[1]] = arr
            return new_arr

        def _repad_X(data):
            # unify the size of dimension -1 for all X arrays in the dataset
            W = [self._cast(sample.X, "2d").shape[-1] for sample in data]

            if len(set(W)) != 1:
                msg = "<intersample_unique>: the dataset contains X arrays of different width; they will be padded to unify. . ."
                self.logger.warning(msg)

                new_W = max(W)
                for sample in data:
                    sample.X = _pad_to_new_dim(self._cast(sample.X, "2d"), new_W)

            return data

        def to_void(arr):
            return arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[1])))

        def get_intersample_unique(data):
            # this ensures that X are cast to 2d
            data = _repad_X(data)

            V = []
            idx = []
            # turn to void type for quick comparison
            for i, sample in enumerate(data):
                self._empty_array_check(sample.X, inspect.stack()[0][3])

                V.append(to_void(sample.X))
                idx.append(np.full(sample.size, i, dtype=np.int8))

            V = np.concatenate(V)
            idx = np.concatenate(idx)

            _, inverse, counts = np.unique(V, return_inverse=True, return_counts=True)

            sample_tracker_dict = defaultdict(set)

            for i, row_idx in enumerate(inverse):
                sample_tracker_dict[row_idx.item()].add(idx[i])

            mask = np.array(
                [len(sample_tracker_dict[row.item()]) == 1 for row in inverse]
            )
            offset = 0

            for sample in data:
                n = sample.size
                sample.X = sample.X[mask[offset : offset + n]]
                offset += n

            return data

        return get_intersample_unique

    def filter_external(self, external=None, max_hd=None):
        """
        Filter sequences similar to an external dataset.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Removes sequences that are within a specified Hamming distance of any
        sequence in an external dataset. Useful for removing validation/test
        sequences or known contaminants from training data.

        Args:
            external: External dataset to compare against. Should be castable
                to np.ndarray or compatible with AnalysisSample.

            max_hd (int): Maximum Hamming distance threshold. Sequences with
                distance ≤ max_hd from any external sequence are removed.

        Returns:
            callable: Operation that accepts a Data object, filters sequences
                similar to external dataset, and returns the modified Data object.

        Example:
            >>> #remove sequences similar to validation set
            >>> external_filt = C.preprocessor.filter_external(
            ...     external=validation_sequences, max_hd=2
            ... )
            >>> data = external_filt(data)
        """
        from clibas.misc import hamming_distance

        if not isinstance(max_hd, int):
            msg = "<filter_external> op expected min_hd as type=int; received: {type(min_hd)}"
            self.logger.error(msg)
            raise ValueError(msg)

        # cast w/e is passed as an external dataset to AnalysisSample
        external = AnalysisSample(X=external).X

        def filter_external_X(data):
            for sample in data:
                arr = self._cast(sample.X, "2d")
                self._empty_array_check(sample.X, inspect.stack()[0][3])

                dtype = arr.dtype
                ext = self._cast(external, "2d").astype(dtype)

                if ext.shape[-1] != arr.shape[-1]:
                    msg = f"<filter_external>: the shape of external dataset along dimension 2 ({ext.shape[-1]}) does not match the data ({arr.shape[-1]})!"
                    self.logger.error(msg)
                    raise ValueError(msg)

                for pep in ext:
                    to_pop = hamming_distance(
                        arr, pep, max_hd, cum=True, return_index=True
                    )

                    ind = np.ones(sample.size, dtype=bool)
                    ind[to_pop] = False
                    sample.ind_filter(ind)

            return data

        return filter_external_X

    def merge(self):
        """
        Merge all samples into a single dataset.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Combines all samples in the Data object into a single merged sample.
        The merged sample is named 'merged_data'.

        Returns:
            callable: Operation that accepts a Data object, merges all samples, and returns the modified Data object containing a single sample.

        Example:
            >>> merge_op = C.preprocessor.merge()
            >>> data = merge_op(data)
            >>> # data now contains a single merged sample
        """

        def merge_datasets(data):
            data.stack()
            data.samples[0].name = "merged_data"

            return data

        return merge_datasets

    def sample(self, sample_size=None):
        """
        Randomly sample from each dataset.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Randomly samples a specified number or fraction of sequences from each
        sample in the dataset. Applied independently to each sample.

        Args:
            sample_size (int or float): Number or fraction of sequences to sample.
                If ≤ 1, interpreted as fraction of dataset to keep.
                If > 1, interpreted as absolute number of sequences to sample.

        Returns:
            callable: Operation that accepts a Data object, samples from each dataset, and returns the modified Data object.

        Example:
            >>> #sample 50% of each dataset
            >>> sample_op = C.preprocessor.sample(sample_size=0.5)
            >>> data = sample_op(data)
            >>>
            >>> #sample exactly 1000 sequences from each dataset
            >>> sample_op = C.preprocessor.sample(sample_size=1000)
            >>> data = sample_op(data)
        """
        if not (isinstance(sample_size, int) or isinstance(sample_size, float)):
            msg = "<sample_from_datasets> op expected sample_size as type=int or float; received: {type(sample_size)}"
            self.logger.error(msg)
            raise ValueError(msg)

        def take_a_sample(data):
            for sample in data:
                if sample_size <= 1:
                    size = int(sample_size * sample.size)
                else:
                    size = int(sample_size)

                if sample.size < size:
                    msg = f"Cannot take a sample that is bigger than the dataset. Sampling is ignored for {sample} sample."
                    self.logger.warning(msg)
                    continue

                ind = np.random.choice(sample.size, size=size, replace=False)
                sample.ind_filter(ind)

            return data

        return take_a_sample

    def shuffle(self):
        """
        Randomly shuffle sequences within each sample.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Randomly reorders sequences within each sample's X dataset. Applied
        independently to each sample.

        Returns:
            callable: Operation that accepts a Data object, shuffles sequences within each sample, and returns the modified Data object.

        Example:
            >>> shuffle_op = C.preprocessor.shuffle()
            >>> data = shuffle_op(data)
        """

        def shuffle_intraset(data):
            for sample in data:
                ind = np.arange(sample.X.shape[0])
                np.random.shuffle(ind)
                sample.ind_filter(ind)

            return data

        return shuffle_intraset

    # TODO: not very happy with having to deepcopy the sample in memory,
    # even if briefly; but it isn't a major issue at the moment.
    def tt_split(self, test_fraction=None):
        """
        Perform train/test split on dataset.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Splits a single sample into separate train and test datasets. The input
        Data object must contain exactly one sample.

        Args:
            test_fraction (float): Fraction of data to allocate to test set
                (between 0 and 1). Remaining data goes to training set.

        Returns:
            callable: Operation that accepts a Data object with one sample, splits it, and returns a Data object with two samples named 'train_data' and 'test_data'.

        Raises:
            ValueError: If input Data contains more than one sample.

        Example:
            >>> #split into 80% train, 20% test
            >>> split_op = C.preprocessor.tt_split(test_fraction=0.2)
            >>> data = split_op(data)
            >>> #data now contains train_data and test_data samples
        """
        if not (isinstance(test_fraction, int) or isinstance(test_fraction, float)):
            msg = "<test_train_split> op expected test_fraction as type=int or float; received: {type(test_fraction)}"
            self.logger.error(msg)
            raise ValueError(msg)

        def test_train_split(data):
            if data.size != 1:
                msg = "A single sample must be passed to the <test_train_split> op; received: {data.size} samples. . ."
                self.logger.error(msg)
                raise ValueError(msg)

            sample = data.samples[0]

            full_set_size = sample.size
            test_set_size = int(test_fraction * full_set_size)
            test_set_ind = np.random.choice(
                full_set_size, size=test_set_size, replace=False
            )

            mask = np.ones(full_set_size, dtype=bool)
            mask[test_set_ind] = False

            train_sample = sample
            test_sample = copy.deepcopy(sample)

            train_sample.name = "train_data"
            train_sample.ind_filter(mask)

            test_sample.name = "test_data"
            test_sample.ind_filter(~mask)

            data = Data(samples=[train_sample, test_sample])

            return data

        return test_train_split

    def to_h5(
        self, F=None, alphabet=None, reshape=False, chunks=None, return_data=False
    ):
        """
        Featurize sequences and save to HDF5 files.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Converts sequence datasets to numerical feature representations and saves
        them as HDF5 files. Useful when featurized datasets are too large to fit
        in memory. Processes data in chunks for memory efficiency.

        Args:
            F (str, ndarray, or None): Feature matrix specification. If ``None``,
                uses one-hot encoding. String options include 'varimax',
                'pep_ECFP3', 'pep_ECFP4', 'pep_SMILES'. Can provide custom
                2D array with ``F.shape[0] == len(alphabet)``.

            alphabet (tuple, list, or ndarray, optional): Token alphabet. If ``None``,
                uses peptide alphabet from configuration.

            reshape (bool): If True, represents sequences as multidimensional
                tensors. If False, unrolls to vectors. Default is False.

            chunks (int): Process data in chunks of this size for memory efficiency.

            return_data (bool): If True, returns unmodified Data object. If False,
                returns None. Default is False.

        Returns:
            callable: Operation that accepts a Data object, featurizes sequences to HDF5 files, and returns Data object or ``None`` based on return_data.

        Note:
            HDF5 files are saved in the ml_data directory specified in the config file.

        Example:
            >>> to_h5_op = C.preprocessor.to_h5(F='pep_ECFP4', chunks=20, reshape=True)
            >>> to_h5_op(data)
        """
        if not isinstance(reshape, bool):
            msg = "<featurize_to_h5> op expected param reshape as type=bool; received: {type(reshape)}"
            self.logger.error(msg)
            raise ValueError(msg)

        if not isinstance(chunks, int):
            msg = "<featurize_to_h5> op expected param chunks as type=int; received: {type(chunks)}"
            self.logger.error(msg)
            raise ValueError(msg)

        if not isinstance(return_data, bool):
            msg = "<featurize_to_h5> op expected param return_data as type=bool; received: {type(return_data)}"
            self.logger.error(msg)
            raise ValueError(msg)

        alphabet = self._infer_alphabet(alphabet=alphabet)

        # no need to pre-validate F -> FeatureMatrix will take care of it
        from clibas.featurization import FeatureMatrix

        F = FeatureMatrix.make(descr=F, constants=self.constants).F

        def featurize_to_h5(data):
            for sample in data:
                arr = self._joint_alphabet_X_check(sample.X, alphabet)
                self._empty_array_check(arr, inspect.stack()[0][3])
                self._prepare_destinations(root=self.dirs.ml_data)
                path = os.path.join(self.dirs.ml_data, f"{sample.name}.hdf5")
                
                featurize.into_h5(
                    arr,
                    y=sample.y,
                    alphabet=alphabet,
                    path=path,
                    F=F,
                    reshape=reshape,
                    chunks=chunks,
                )
            if return_data:
                return data

            return

        return featurize_to_h5

    def featurize_X(self, F=None, alphabet=None, reshape=False):
        """
        Featurize sequence datasets in memory.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Converts sequence datasets to numerical feature representations. Use when
        featurized datasets fit in memory. For large datasets, use to_h5() instead.

        Args:
            F (str, ndarray, or None): Feature matrix specification. If ``None``,
                uses one-hot encoding. String options include 'varimax',
                'pep_ECFP3', 'pep_ECFP4', 'pep_SMILES'. Can provide custom
                2D array with ``F.shape[0] == len(alphabet)``.

            alphabet (tuple, list, or ndarray, optional): Token alphabet. If ``None``,
                uses peptide alphabet from the config file.

            reshape (bool): If True, represents sequences as multidimensional
                tensors. If False, unrolls to vectors. Default is False.

        Returns:
            callable: Operation that accepts a Data object, featurizes X datasets, and returns the modified Data object.

        Example:
            >>> featurize_op = C.preprocessor.featurize_X(F='pep_ECFP4', reshape=True)
            >>> data = featurize_op(data)
        """
        if not isinstance(reshape, bool):
            msg = "<featurize_X_datasets> op expected param reshape as type=bool; received: {type(reshape)}"
            self.logger.error(msg)
            raise ValueError(msg)

        alphabet = self._infer_alphabet(alphabet=alphabet)

        # no need to pre-validate F -> FeatureMatrix will take care of it
        from clibas.featurization import FeatureMatrix

        F = FeatureMatrix.make(descr=F, constants=self.constants).F

        def featurize_X_datasets(data):
            for sample in data:
                arr = self._joint_alphabet_X_check(sample.X, alphabet)
                self._empty_array_check(arr, inspect.stack()[0][3])

                sample.X = featurize.from_matrix_v3(
                    arr, F=F, alphabet=alphabet, reshape=reshape
                )
            return data

        return featurize_X_datasets

    def featurize_for_RFA(self, alphabet=None, order=None):
        """
        Featurize sequences for Reference-Free Analysis (RFA) models as implemented
        in DOMEK workflows. See https://www.cell.com/chem/abstract/S2451-9294(25)00328-6

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Converts peptide sequences to feature representations compatible with
        RFA models. See https://www.nature.com/articles/s41467-024-51895-5
        for methodology details.

        Args:
            alphabet (list, tuple, or ndarray): Amino acid alphabet comprising
                the sequences. Alphabet size determines feature dimensions.

            order (str): RFA model order. Must be 'first' or 'second'. Higher
                orders not yet implemented.

        Returns:
            callable: Operation that accepts a Data object, featurizes sequences for RFA, and returns the modified Data object with flattened 2D feature arrays.

        Example:
            >>> rfa_op = C.preprocessor.featurize_for_RFA(alphabet='aa', order='second')
            >>> data = rfa_op(data)
        """
        if order != "first" and order != "second":
            msg = '<featurize_for_RFA>: the "order" keyword not understood: only "first" or "second" are allowed values!'
            self.logger.error(msg)
            raise ValueError(msg)

        alphabet = self._infer_alphabet(alphabet=alphabet)

        def RFA_featurization(data):
            for sample in data:
                arr = self._joint_alphabet_X_check(sample.X, alphabet)
                self._empty_array_check(arr, inspect.stack()[0][3])

                sample.X = featurize.RFA_featurization(
                    arr, alphabet=alphabet, order=order
                )
            return data

        return RFA_featurization

    def drop(self, sample_to_drop=None):
        """
        Remove a sample from the dataset by name.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Removes the specified sample from the Data object. If sample name is not
        found, logs a warning and returns data unchanged.

        Args:
            sample_to_drop (str): Name of the sample to remove.

        Returns:
            callable: Operation that accepts a Data object, removes the specified sample, and returns the modified Data object.

        Example:
            >>> drop_op = C.preprocessor.drop(sample_to_drop='validation_data')
            >>> data = drop_op(data)
        """
        if not isinstance(sample_to_drop, str):
            msg = "<drop_dataset> op expected param dataset_to_drop as type=str; received: {type(dataset_to_drop)}"
            self.logger.error(msg)
            raise ValueError(msg)

        def drop_dataset(data):
            to_drop = []
            for i, sample in enumerate(data):
                if sample.name == sample_to_drop:
                    to_drop.append(i)

            if not to_drop:
                msg = f"<drop_dataset>: no {sample_to_drop} in the dataset. . ."
                self.logger.warning(msg)

            for i in to_drop:
                del data.samples[i]

            return data

        return drop_dataset

    def pad_and_random_shift(self, new_x_dim=None):
        """
        Expand and randomly shift sequences for data augmentation.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Pads each sequence to a fixed width and applies a random circular shift.
        The shift is constrained so that non-pad values remain within bounds,
        effectively redistributing padding without truncating sequence data.
        Useful for data augmentation and positional invariance.

        Args:
            new_x_dim (int): Target width for padded sequences. Must be ≥ current
                sequence width.

        Returns:
            callable: Operation that accepts a Data object, pads and shifts sequences, and returns the modified Data object.
        Example:
            If new_x_dim=6, a sequence array::

                  [['A', 'B', 'C', 'D'],
                   ['E', 'F', 'G', 'H'],
                   ['I', 'J', 'K', 'L'],
                   ['M', 'N', 'O', 'P']]

            might become:

                  [['', '', 'A', 'B', 'C', 'D'],
                   ['', '', 'E', 'F', 'G', 'H'],
                   ['I', 'J', 'K', 'L', '', ''],
                   ['', 'M', 'N', 'O', 'P', '']]

            with padding randomly distributed on either side.

        >>> augment_op = C.preprocessor.pad_and_random_shift(new_x_dim=20)
        >>> data = augment_op(data)
        """
        if not isinstance(new_x_dim, int):
            msg = f'<x_expand_and_shift> op expended "new_x_dim" param as dtype=int; received: {type(new_x_dim)}'
            self.logger.error(msg)
            raise ValueError(msg)

        def expand_and_random_shift(data):
            for sample in data:
                arr = self._cast(sample.X, "2d")
                self._empty_array_check(arr, inspect.stack()[0][3])

                dtype = arr.dtype
                pad = dtype.type()

                expanded_arr = np.zeros((arr.shape[0], new_x_dim), dtype=dtype)
                expanded_arr[:, : arr.shape[1]] = arr

                L = np.sum(arr != pad, axis=1)
                max_shift = new_x_dim - L + 1

                # this isn't easy to understand, but it works. essentially,
                # we are doing a circular permutation by a random amount but
                # not larger than the number of pads on the right side of a given
                # row. what this does is reshuffle some of the pads to the left
                shift_ind = np.random.randint(0, high=max_shift, size=arr.shape[0])
                idx = np.mod(np.arange(new_x_dim) - shift_ind[:, None], new_x_dim)
                sample.X = expanded_arr[np.arange(expanded_arr.shape[0])[:, None], idx]

            return data

        return expand_and_random_shift
