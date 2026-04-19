"""
Core data structures for clibas pipelines.

Provides Sample and Data container classes for managing sequencing data,
machine learning datasets, and pipeline operations. These structures
maintain data integrity and enable consistent transformations throughout
processing workflows.
"""

import copy

import numpy as np


class Sample:
    """
    Base class for data sample types.

    Provides common functionality for different sample types (SequencingSample,
    AnalysisSample). Ensures all array attributes maintain consistent sizes and
    provides utilities for filtering, transforming, and manipulating data arrays.

    Sample object stores component datasets as numpy arrays, wiuth every row
    correponding to a dataset entry (NGS read, peptide sequence, etc).
    Datasets are kept track of together along axis=0. Removing a row from
    one array will remove a corresponding row from array2.

    TODO: this data structure can be replaced with xarrays:: consider porting

    Note:
        This is a base class and should not be invoked directly. Use
        subclasses like SequencingSample or AnalysisSample instead.

    Args:
        name (str): Sample identifier. Default is 'unnamed_sample'.
        **kwargs: Array data as keyword arguments. Each becomes a sample attribute.
    """

    def __init__(self, name="unnamed_sample", **kwargs):
        for key, value in kwargs.items():
            if key != "name":
                self.__dict__.update({key: np.asarray(value)})

        self.name = name
        return

    def __len__(self):
        L = []
        for arr in self._arrays():
            try:
                L.append(arr.shape[0])
            except:
                pass
        if L:
            return L[0]
        else:
            return 0

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, key, value):
        setattr(self, key, value)
        return

    # make sure that when new attributes are set, they don't break
    # the dataset validity (equal size for all arrays).
    def __setattr__(self, key, value):
        if key != "name":
            super().__setattr__(key, np.asarray(value))
        else:
            super().__setattr__(key, value)
        self._validate_size()

    # TODO: we still have __call__() method available

    @property
    def size(self):
        return self.__len__()

    def iterarrays(self):
        # iterate over publicly available arrays
        # will ignore the _internal_state array
        for attr in self.__dict__:
            if type(getattr(self, attr)) == np.ndarray:
                if not attr.startswith("_"):
                    arr = getattr(self, attr)
                    yield (arr, attr)

    # easy way to fetch data arrays only (np.ndarrays)
    def _arrays(self):
        return [
            getattr(self, attr)
            for attr in self.__dict__
            if type(getattr(self, attr)) == np.ndarray
        ]

    def _validate_name(self):
        if not isinstance(self.name, str):
            raise ValueError(f"Sample name {self.name} was not understood. . .")
        return

    def _validate_size(self):
        L = []
        for arr in self._arrays():
            try:
                L.append(arr.shape[0])
            except:
                pass

        assert len(set(L)) < 2, "Sample datasets cannot contain data of different size!"
        return

    def _transform(self, item, to_type=None):
        if not type(getattr(self, item)) == np.ndarray:
            raise TypeError(
                "Unsupported array type passed to Sample.transform function."
            )

        if to_type is not None:
            if to_type not in ("S", "U"):
                raise ValueError(
                    f"Sample.transform function expected 'to_type' as either 'S' or 'U'; received: {to_type}"
                )

        arr = getattr(self, item)

        if arr.ndim == 0:
            return

        arr = np.ascontiguousarray(arr)
        kind = arr.dtype.kind

        # currently, we can recast only SX, UX and object-type arrays
        if arr.ndim == 1 and arr.dtype.kind == "O":
            arr = arr.astype(f"U{max(len(s) for s in arr)}")

        if arr.ndim == 1 and kind in ("S", "U"):
            arr = arr.view(f"{kind}1").reshape((arr.shape[0], -1))
            if to_type is not None:
                arr = arr.astype(f"{to_type}1")

        elif arr.ndim > 1 and to_type is not None:
            if arr.dtype in ("|S1", "<U1"):
                arr = arr.astype(f"{to_type}1")

        setattr(self, item, arr)
        return

    def ind_filter(self, ind, in_place=True, new_name=None):
        """
        Filter sample arrays using boolean or integer indexing.

        Applies the same indexing operation to all array attributes
        simultaneously, maintaining consistency across datasets. Can operate
        in-place or return a new filtered copy. Useful for pipeline operations
        or subset selection.

        Args:
            ind (array-like): Boolean mask or integer indices for filtering.
            in_place (bool, optional): If True, modifies the current sample.
                If False, returns a new filtered sample. Default is True.
            new_name (str, optional): Name for the new sample. Only used if
                `in_place=False`. Defaults to original sample name.

        Returns:
            Sample or None: Returns a new Sample instance if `in_place=False`.
            Returns None if `in_place=True`.

        Example:
            >>> # In-place filtering (keep first two entries)
            >>> sample.ind_filter(np.arange(2))
            >>>
            >>> # Create a filtered copy with a new name
            >>> new_sample = sample.ind_filter(np.arange(2), in_place=False, new_name="subset")
            >>> new_sample.name
            'subset'
        """
        if in_place:
            # modify current object
            for attr, arr in self.__dict__.items():
                if isinstance(arr, np.ndarray) and arr.ndim > 0:
                    # note that setattr can't be used as it would fail
                    # the _validate_size check
                    self.__dict__[attr] = arr[ind]
            return

        # TODO: not super happy with in-memory copying but should do for now
        # may need to make SequencingSample object more mutable in the
        # future to avoid this
        new_obj = copy.copy(self)

        for attr, arr in self.__dict__.items():
            if isinstance(arr, np.ndarray):
                if arr.ndim > 0:
                    new_obj.__dict__[attr] = arr[ind]
                else:
                    new_obj.__dict__[attr] = arr.copy()

        if new_name is not None:
            new_obj.name = new_name

        new_obj._validate_size()
        return new_obj

    def drop(self, item):
        """
        Remove an array attribute from the sample.

        Args:
            item (str): Name of the array attribute to drop.
        """
        if not type(getattr(self, item)) == np.ndarray:
            return

        setattr(self, item, np.asarray(None))
        return

    def unpad(self, item):
        """
        Remove padding columns from a specific array.

        Removes columns where all values are padding tokens (dtype default value).

        Args:
            item (str): Name of the array attribute to unpad.
        """
        arr = getattr(self, item)
        if type(arr) == np.ndarray:
            pad = arr.dtype.type()

            if arr.ndim == 2:
                setattr(self, item, arr[:, ~np.all(arr == pad, axis=0)])

            elif arr.ndim == 1:
                if arr.dtype.kind == "S":
                    arr = np.ascontiguousarray(arr)
                    arr = arr.view("|S1").reshape((arr.shape[0], -1))
                    arr = arr[:, ~np.all(arr == pad, axis=0)]
                    arr = np.ascontiguousarray(arr)
                    setattr(
                        self, item, arr.view(f"|S{arr.shape[1]}").reshape(arr.shape[0])
                    )
        return

    def unpad_arrays(self):
        """Remove padding from all public array attributes."""
        for _, attr in self.iterarrays():
            self.unpad(attr)
        return

    @classmethod
    def from_dict(cls, data_dict):
        """
        Create sample from dictionary of arrays.

        Args:
            data_dict (dict): Dictionary with array data as key-value pairs.

        Returns:
            Sample: New sample instance with data from dictionary.
        """
        return cls(**data_dict)


class SequencingSample(Sample):
    """
    Container for DNA sequencing data.

    Stores DNA sequences, quality (Q) scores, and translated peptide sequences
    with associated metadata. Maintains internal state for tracking library
    design compatibility during filtering operations.

    Args:
        dna (array-like, optional): DNA sequences.
        Q (array-like, optional): Quality scores (Phred format).
        pep (array-like, optional): Translated peptide sequences.
        name (str): Sample identifier. Default is 'unnamed'.

    Attributes:
        dna (ndarray): DNA sequence data (1D or 2D).
        Q (ndarray): Quality score data (1D or 2D).
        pep (ndarray): Peptide sequence data (1D or 2D).
        _internal_state (ndarray): Tracks library design compatibility.

    Note:
        Arrays are usually stored as 1D byte-strings. Pipeline operations
        have the ability to quickly recast there between 1D and 2D representations
        for processing. 2D arrays are padded to accommodate variable-length
        sequences. All arrays maintain equal first dimension (number of entries).
    """

    def __init__(self, dna=None, Q=None, pep=None, name="unnamed"):
        super(SequencingSample, self).__init__(dna=dna, Q=Q, pep=pep, name=name)
        self._validate_name()
        self._validate_size()
        self._validate_dims()

        self._internal_state = None
        return

    def __repr__(self):
        return f"<SequencingSample {self.name} containing {len(self)} entries>"

    def _collapse_internal_state(self):
        """
        Resolve multi-template matches to single template per entry.

        Samples maintain internal state tracking library design compatibility
        during filtering. Some operations (like fetch_at) need to "collapse" this
        state, assigning each entry to a single template.

        When an entry matches multiple library design templates, selects one
        template and marks others as incompatible. Used during pipeline operations
        requiring unambiguous template assignments.
        """
        if self._internal_state.ndim > 0:
            collapsed = np.zeros_like(self._internal_state, dtype=bool)
            ind = (
                np.arange(self._internal_state.shape[0]),
                self._internal_state.argmax(axis=1),
            )
            collapsed[ind] = self._internal_state[ind]
            self._internal_state = collapsed
        return

    def _is_collapsed(self):
        """
        Check whether the internal state array is collapsed.
        The idea is that depending on the pipeline, a sample's
        internal state may collapse without explicitly calling
        self._collapse_internal_state (which is a good thing).
        """
        if self._internal_state.ndim > 1:
            if np.all(self._internal_state.sum(axis=1) == 1):
                return True
        return

    def _validate_dims(self):
        dim = []
        for arr in [self.dna, self.Q, self.pep]:
            if arr.ndim:
                dim.append(arr.ndim)

        if len(set(dim)) > 1:
            raise ValueError(
                "Sample datasets contain data of different dimensionality!"
            )

        if len(set(dim)) == 1:
            return dim[0]

        else:
            return 0

    @property
    def civilized_Q(self):
        """
        Convert quality scores to numerical Phred values.

        Transforms ASCII-encoded quality scores to numerical values by
        subtracting 33 (Illumina/Sanger encoding). Padding bytes (0) are
        preserved as 0.

        Returns:
            ndarray: Numerical quality scores (uint8).
        """
        if self.Q.dtype.kind == "S" and self.Q.ndim == 2:
            civ = self.Q

        elif self.Q.dtype.kind == "S" and self.Q.ndim == 1:
            civ = self.Q.view("|S1").reshape((self.Q.shape[0], -1))

        else:
            raise ValueError("Could not civilize the Q score array. . .")

        return np.where(civ.view(np.uint8) == 0, 0, civ.view(np.uint8) - 33)

    def get_ndims(self):
        """
        Get dimensionality of sample arrays.

        Returns:
            int: Number of dimensions (0, 1, or 2).
        """
        return self._validate_dims()


class AnalysisSample(Sample):
    """
    Container for machine learning and analysis data.

    Generic sample type for data analysis and machine learning workflows.
    Stores feature arrays (X) with optional labels (y) and additional metadata.

    Args:
        X (array-like, optional): Feature data (sequences or numerical features).
        y (array-like, optional): Labels or target values.
        name (str): Sample identifier. Default is 'unnamed'.
        **kwargs: Additional array data as keyword arguments.

    Attributes:
        X (ndarray): Feature array (automatically converted to byte strings).
        y (ndarray): Label array (if provided).

    Example:
        >>> sample = AnalysisSample(X=sequences, y=labels, name='test_data')
    """

    def __init__(self, X=None, y=None, name="unnamed", **kwargs):
        super(AnalysisSample, self).__init__(X=X, y=y, name=name, **kwargs)
        self._validate_name()
        self._validate_size()

        # at initialization, convert to an S-type
        self._transform("X", to_type="S")
        return

    def get_ndims(self):
        """
        Get dimensionality of X array.

        Returns:
            int: Number of dimensions in X.
        """
        return self.X.ndim

    def __repr__(self):
        return f"<AnalysisSample {self.name} containing {len(self)} entries>"


class Data:
    """
    Container for collections of Sample objects.

    Primary data structure for clibas pipelines. Holds multiple samples of the
    same type (SequencingSample or AnalysisSample) and ensures type consistency.
    All pipeline operations accept and return Data objects.

    Sample types typically represent individual FASTQ files or analysis datasets.
    Data objects maintain sample integrity while enabling batch operations across
    all contained samples.

    Args:
        samples (list, optional): List of Sample objects (all same type).

    Attributes:
        samples (list): List of Sample objects.
        size (int): Number of samples in container.
        dtype (type): Type of samples (e.g., SequencingSample).

    Example:
        >>> data = Data(samples=[sample1, sample2, sample3])
        >>> len(data)  #Number of samples
        3
    """

    def __init__(self, samples=None):
        self.samples = samples
        self._validate()
        return

    def _validate(self):
        if not isinstance(self.samples, (type(None), list, tuple, set)):
            msg = f"Data object expected to receive a list of samples; received {type(self.samples)} instead"
            raise ValueError(msg)

        if self.samples:
            for sample in self.samples:
                if not issubclass(sample.__class__, Sample):
                    msg = f"Unexpected sample types passed to the Data object; type: {type(sample)}"
                    raise ValueError(msg)

            if len(set([type(x) for x in self.samples])) > 1:
                msg = "All samples should have the same sample type!"
                raise ValueError(msg)

        return

    def __iter__(self):
        for sample in self.samples:
            yield sample

    def __len__(self):
        if self.samples:
            return len(self.samples)
        else:
            return 0

    @property
    def size(self):
        """Number of samples in the Data container."""
        return self.__len__()

    @property
    def dtype(self):
        """Type of samples contained (e.g., SequencingSample, AnalysisSample)."""
        if not self.samples:
            return type(None)
        else:
            return type(self.samples[0])

    def __getitem__(self, item):
        return self.samples[item]

    def __repr__(self):
        s = f"<Data container holding {len(self)} samples>:\n"
        if self.size < 10:
            for sample in self.samples:
                s += f"{sample}\n"
        else:
            for sample in self.samples[:5]:
                s += f"{sample}\n"
            s += "...\n"
            for sample in self.samples[-5:]:
                s += f"{sample}\n"

        return s

    # TODO: make sure that the situations where only some y arrays
    # are present are dealt with; it shouldn't stack in such cases
    def stack(self, in_place=True):
        """
        Combine all samples into a single merged sample.

        Stacks arrays from all samples vertically (row-wise). Arrays with
        different widths are padded to match the maximum width. All samples
        must have the same set of array attributes.

        Args:
            in_place (bool): If True, modifies this Data object to contain only
                the stacked sample. If False, returns a new Data object with the
                stacked sample. Default is True.

        Returns:
            Data or None: New Data object with stacked sample if in_place=False,
                otherwise None (modifies in place).

        Raises:
            ValueError: If samples have different array attributes or incompatible
                array dimensions.

        Example:
            >>> #merge all samples into one
            >>> data.stack(in_place=True)
            >>> len(data)  #now contains 1 sample
            1
            >>> #create new merged dataset without modifying original
            >>> merged_data = data.stack(in_place=False)
        """

        def _new_width_finder(arrays):
            if len(set([arr.ndim for arr in arrays])) != 1:
                raise ValueError(
                    "All arrays of one type should be same-dimensional to stack!"
                )

            # SX and UX arrays are out of the question, because of the transform below
            if arrays[0].ndim == 0:
                return

            if arrays[0].ndim == 1:
                return "concat"

            if arrays[0].ndim > 2:
                raise ValueError(
                    "Stacking higher dimensional samples is not supported at the moment."
                )

            # this is the two-dimensional case
            new_dim = [arr.shape[-1] for arr in arrays]
            return max(new_dim)

        def _pad_to_new_dim(arr, new_dim):
            to_pad = new_dim - arr.shape[1]
            if to_pad == 0:
                return arr

            new_arr = np.zeros((arr.shape[0], new_dim), dtype=arr.dtype)
            new_arr[:, : arr.shape[1]] = arr
            return new_arr

        if not len(self.samples):
            return None

        # need to make sure that every sample has the same set of arrays,
        # otherwise, how would it work?
        attrs = [tuple(sorted(sample.__dict__.keys())) for sample in self.samples]

        if len(set(attrs)) != 1:
            raise ValueError(
                "Samples can only be stacked if they all hold the same set of arrays!"
            )

        stacked = {}
        for attr in attrs[0]:
            if not type(getattr(self.samples[0], attr)) == np.ndarray:
                continue

            to_stack = []
            for sample in self.samples:
                sample._transform(attr)
                to_stack.append(getattr(sample, attr))

            new_width = _new_width_finder(to_stack)

            if not new_width:
                stacked[attr] = np.array(None)

            elif new_width == "concat":
                stacked[attr] = np.concatenate(to_stack, axis=0)

            else:
                to_stack = [_pad_to_new_dim(arr, new_width) for arr in to_stack]
                stacked[attr] = np.vstack(to_stack)

        stacked["name"] = "stacked_samples"
        sample_type = self.samples[0].__class__
        if in_place:
            self.samples = [sample_type.from_dict(stacked)]
            return

        return Data([sample_type.from_dict(stacked)])
