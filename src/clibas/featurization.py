"""
Sequence featurization utilities for machine learning.

Provides tools for converting sequence data (peptides, DNA) to numerical
feature representations using one-hot encoding, chemical descriptors (SMILES,
ECFPs), or biophysical properties (varimax).
"""

import os

import numpy as np


def _get_SMILES_chars(monomer_SMILES):
    """
    Extract and enumerate unique characters from SMILES strings.

    Generates sorted list and enumeration dictionary of all unique characters
    found in the SMILES tuple.

    Args:
        monomer_SMILES (tuple): SMILES strings for monomers.

    Returns:
        tuple: (char_list, char_dict) where:
            - char_list: Sorted list of unique characters
            - char_dict: Dictionary mapping characters to indices (1-based)
    """
    ch = "".join(x for x in monomer_SMILES)
    ch = set(list(ch))
    ch = list(ch)
    ch.sort()
    return ch, {c: i + 1 for i, c in enumerate(ch)}


class FeatureMatrix:
    """
    Feature matrix manager for sequence featurization.

    Handles creation, validation, and management of feature matrices that
    transform sequences into numerical representations for machine learning.
    Supports various featurization schemes including one-hot encoding,
    SMILES-based representations, Morgan fingerprints, and varimax features.

    Args:
        F (ndarray, optional): Feature matrix. Shape: (alphabet_size, n_features).
        constants: Config constants object with alphabet and (optionally) SMILES data.

    Attributes:
        F (ndarray): The feature matrix used for transformations.
    """

    def __init__(self, F=None, constants=None):
        self.F = F
        self.constants = constants
        self._validate()
        return

    def __repr__(self):
        return f"<FeatureMatrix object>\nF={self.F}"

    def _validate(self):
        if self.F is not None:
            if not isinstance(self.F, np.ndarray):
                raise ValueError(
                    f"<FeatureMatrix>: feature matrix expected as a numpy array; received: {type(self.F)}"
                )
            else:
                if not self.F.ndim == 2:
                    raise ValueError(
                        f"<FeatureMatrix> expected to receive parameter F as a two-dimensional array; received: ndims={self.F.ndim}"
                    )

        return

    @classmethod
    def make(cls, descr=None, constants=None):
        """
        Create a feature matrix from various specifications.

        Factory method that constructs feature matrices from descriptions,
        file paths, or returns one-hot encoding (None).

        Args:
            descr (str, ndarray, or None): Feature matrix specification.
                - None: One-hot encoding (identity matrix)
                - ndarray: Use provided matrix directly
                - File path: Load matrix from .npy file
                - 'varimax': Varimax biophysical features (peptides only)
                - 'pep_ECFP3', 'pep_ECFP4': Morgan fingerprints for peptides
                - 'pep_SMILES': SMILES character encoding for peptides
                - 'dna_SMILES': SMILES character encoding for DNA

            constants: config file constants with alphabet and (optionally) SMILES data.

        Returns:
            FeatureMatrix: Configured feature matrix instance.

        Example:
            >>> #one-hot encoding
            >>> F = FeatureMatrix.make(descr=None, constants=config.constants)
            >>>
            >>> #varimax features
            >>> F = FeatureMatrix.make(descr='varimax', constants=config.constants)
            >>>
            >>> #ECFP
            >>> F = FeatureMatrix.make(descr='pep_ECFP4', constants=config.constants)
        """
        if descr is None:
            return cls(F=None, constants=constants)

        elif isinstance(descr, np.ndarray):
            return cls(F=descr, constants=constants)

        elif isinstance(descr, str):
            if os.path.isfile(descr):
                try:
                    F = np.load(descr)
                except:
                    raise ValueError(
                        "<FeatureMatrix.make> could not load the matrix from the specified location. . ."
                    )

            # if all easy setups failed, resort to string parsing
            else:
                # specific to proteinogenic amino acids                
                if descr == "varimax":
                    F = cls.varimax_repr(constants.aas)
                    return cls(F=F, constants=constants)
                
                prefix = descr.split("_")[0]
                if prefix == "pep":
                    SMILES = constants.aa_SMILES
                elif prefix == "dna":
                    SMILES = constants.base_SMILES
                else:
                    raise ValueError(
                        f"<FeatureMatrix.make> did not understand feature matrix specification: {descr}"
                    )

                if "ECFP" in descr:
                    r = int(descr[-1])
                    F, _, _ = cls.dense_morgan(r, SMILES)

                elif "SMILES" in descr:
                    F = cls.SMILES_repr_v2(SMILES)

                else:
                    raise ValueError(
                        "<FeatureMatrix.make> did not understand feature matrix specification. . ."
                    )

        return cls(F=F, constants=constants)

    def save(self, fname):
        """
        Save feature matrix to file.

        Args:
            fname (str): Output filename (.npy format).
        """
        if not self.F:
            raise ValueError("No feature matrix to be saved!")
        else:
            np.save(fname, self.F)

        return

    def single_linkage_tree(self, labels=None, fname=None):
        """
        Generate hierarchical clustering dendrogram of features.

        Args:
            labels (list, optional): Labels for features.
            fname (str, optional): Output filename for plot.
        """

        from scipy.cluster.hierarchy import linkage
        import clibas.plotters as P

        if labels is not None:
            try:
                labels = [a.decode("ascii") for a in labels]
            except KeyError:
                pass

        if self.F is not None:
            P.Miscellaneous.single_linkage_dendrogram(
                linkage(self.F, method="ward"), labels=labels, basename=fname
            )
        return

    @staticmethod
    def dense_morgan(r, monomer_SMILES):
        """
        Generate Morgan fingerprint (ECFP) feature matrix.

        Creates dense bit-string representation of Morgan (circular) fingerprints
        for molecular structures defined by SMILES strings.

        Args:
            r (int): Maximum fingerprint radius.
            monomer_SMILES (tuple): SMILES strings for each monomer.

        Returns:
            tuple: (F, bit_features, info) where:
                - F: Feature matrix (n_monomers, n_features)
                - bit_features: Fingerprint indices (internal RDKit representation)
                - info: Fingerprint metadata (for structure mapping)
        """
        try:
            from rdkit import Chem          
            from rdkit.Chem import rdFingerprintGenerator
            
        except ImportError:
            msg = "Failed to import rdkit packages. Please install clibas with `pip install clibas[ml]`. . ."
            raise ImportError(msg)

        aas = [Chem.MolFromSmiles(x) for x in monomer_SMILES]
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=r)

        # construct a list of all bit features
        bit_features = []
        for aa in aas:
            fingerprints = generator.GetSparseCountFingerprint(aa)
            fingerprints_dict = fingerprints.GetNonzeroElements()
            keys = list(fingerprints_dict.keys())

            for k in keys:
                bit_features.append(k)

        bit_features = list(set(bit_features))

        # assemble the F matrix, encoding fingerprints as a dense bit string
        F = np.zeros((len(monomer_SMILES), len(bit_features)))
        info = []
        for i, aa in enumerate(aas):
            fp_info = rdFingerprintGenerator.AdditionalOutput()
            fp_info.CollectBitInfoMap()

            fingerprints_sparse = generator.GetSparseCountFingerprint(
                aa, additionalOutput=fp_info
            )
            fingerprints = fingerprints_sparse.GetNonzeroElements()

            for f in fingerprints:
                F[i, bit_features.index(f)] = 1

            info.append(fp_info.GetBitInfoMap())

        return F, bit_features, info

    @staticmethod
    def SMILES_repr_v2(monomer_SMILES):
        """
        Create SMILES character-based feature matrix.

        Deprecated. Each SMILES character assigned numeric index, right-padded
        to longest SMILES string.

        Args:
            monomer_SMILES (tuple): SMILES strings for each monomer.

        Returns:
            ndarray: Feature matrix (n_monomers, max_smiles_length).
        """
        _, SMILES_d = _get_SMILES_chars(monomer_SMILES)
        x_dim = len(monomer_SMILES)
        y_dim = max([len(x) for x in monomer_SMILES])

        F = np.zeros((x_dim, y_dim)) - 1
        for i, ch in enumerate(monomer_SMILES):
            ch_as_list = list(ch)
            F[i, : len(ch_as_list)] = [SMILES_d[x] for x in ch_as_list]

        return F

    @staticmethod
    def varimax_repr(alphabet):
        """
        Create Varimax biophysical feature matrix for amino acids.

        Based on J. Comput. Biol. (2009) 16, 703-723. Represents each amino
        acid as an 8D vector encoding biophysical properties. Only works for
        standard proteinogenic amino acids.

        Args:
            alphabet (tuple): Amino acid single-letter codes.

        Returns:
            ndarray: Feature matrix (n_amino_acids, 8).
        """
        varimax_1 = {
            "A": 0.57,
            "C": 2.66,
            "D": -2.46,
            "E": -3.08,
            "F": 3.12,
            "G": 0.15,
            "H": -0.39,
            "I": 3.10,
            "K": -3.89,
            "L": 2.72,
            "M": 1.89,
            "N": -2.02,
            "P": -0.58,
            "Q": -2.54,
            "R": -2.80,
            "S": -1.10,
            "T": -0.65,
            "V": 2.64,
            "W": 1.89,
            "Y": 0.79,
        }

        varimax_2 = {
            "A": 3.37,
            "C": -1.52,
            "D": -0.66,
            "E": 3.45,
            "F": 0.68,
            "G": -3.49,
            "H": 1.00,
            "I": 0.37,
            "K": 1.47,
            "L": 1.88,
            "M": 3.88,
            "N": -1.92,
            "P": -4.33,
            "Q": 1.82,
            "R": 0.31,
            "S": -2.05,
            "T": -1.60,
            "V": 0.03,
            "W": -0.09,
            "Y": -2.62,
        }

        varimax_3 = {
            "A": -3.66,
            "C": -3.29,
            "D": -0.57,
            "E": 0.05,
            "F": 2.40,
            "G": -2.97,
            "H": -0.63,
            "I": 0.26,
            "K": 1.95,
            "L": 1.92,
            "M": -1.57,
            "N": 0.04,
            "P": -0.02,
            "Q": -0.82,
            "R": 2.84,
            "S": -2.19,
            "T": -1.39,
            "V": -0.67,
            "W": 4.21,
            "Y": 4.11,
        }

        varimax_4 = {
            "A": 2.34,
            "C": -3.77,
            "D": 0.14,
            "E": 0.62,
            "F": -0.35,
            "G": 2.06,
            "H": -3.49,
            "I": 1.04,
            "K": 1.17,
            "L": 5.33,
            "M": -3.58,
            "N": -0.65,
            "P": -0.21,
            "Q": -1.85,
            "R": 0.25,
            "S": 1.36,
            "T": 0.63,
            "V": 2.34,
            "W": -2.77,
            "Y": -0.63,
        }

        varimax_5 = {
            "A": -1.07,
            "C": 2.96,
            "D": 0.75,
            "E": -0.49,
            "F": -0.88,
            "G": 0.70,
            "H": 0.05,
            "I": -0.05,
            "K": 0.53,
            "L": 0.08,
            "M": -2.55,
            "N": 1.61,
            "P": -8.31,
            "Q": 0.09,
            "R": 0.20,
            "S": 1.78,
            "T": 1.35,
            "V": 0.64,
            "W": 0.72,
            "Y": 1.89,
        }

        varimax_6 = {
            "A": -0.40,
            "C": -2.23,
            "D": 0.24,
            "E": 0.00,
            "F": 1.62,
            "G": 7.47,
            "H": 0.41,
            "I": -1.18,
            "K": 0.10,
            "L": 0.09,
            "M": 2.07,
            "N": 2.08,
            "P": -1.82,
            "Q": 0.60,
            "R": -0.37,
            "S": -3.36,
            "T": -2.45,
            "V": -2.01,
            "W": 0.86,
            "Y": -0.53,
        }

        varimax_7 = {
            "A": 1.23,
            "C": 0.44,
            "D": -5.15,
            "E": -5.66,
            "F": -0.15,
            "G": 0.41,
            "H": 1.61,
            "I": -0.21,
            "K": 4.01,
            "L": 0.27,
            "M": 0.84,
            "N": 0.40,
            "P": -0.12,
            "Q": 0.25,
            "R": 3.81,
            "S": 1.39,
            "T": -0.65,
            "V": -0.33,
            "W": -1.07,
            "Y": -1.30,
        }

        varimax_8 = {
            "A": -2.32,
            "C": -3.49,
            "D": -1.17,
            "E": -0.11,
            "F": -0.41,
            "G": 1.62,
            "H": -0.60,
            "I": 3.45,
            "K": -0.01,
            "L": -4.06,
            "M": 1.85,
            "N": -2.47,
            "P": -1.18,
            "Q": 2.11,
            "R": 0.98,
            "S": -1.21,
            "T": 3.43,
            "V": 3.93,
            "W": -1.66,
            "Y": 1.31,
        }

        varimax = (
            varimax_1,
            varimax_2,
            varimax_3,
            varimax_4,
            varimax_5,
            varimax_6,
            varimax_7,
            varimax_8,
        )
        
        try:
            alphabet = [a.decode("ascii") for a in alphabet]
        except KeyError:
            pass

        x_dim = len(alphabet)
        y_dim = len(varimax)
        F = np.zeros((x_dim, y_dim))
        for i, v in enumerate(varimax):
            for j, aa in enumerate(alphabet):
                F[j, i] = v[aa]

        return F


def from_matrix_v3(X, alphabet=None, F=None, reshape=True):
    """
    Featurize sequence array using a feature matrix.

    General-purpose featurizer for converting categorical sequence data to
    numerical representations for machine learning. Efficiently handles large
    datasets using matrix multiplication.

    Args:
        X (ndarray): Sequences to featurize. Shape: (n_sequences, seq_length).

        alphabet (array-like): Token alphabet matching X dtype.

        F (ndarray, optional): Feature matrix (alphabet_size, n_features).
            If None, uses one-hot encoding.

        reshape (bool): If True, returns sequences as 2D feature arrays
            (n_sequences, seq_length, n_features). If False, flattens to
            1D vectors. Default is True.

    Returns:
        ndarray: Featurized sequences. Shape depends on reshape parameter.

    Example:
        >>> #one-hot encoding
        >>> X_feat = from_matrix_v3(X, alphabet=aas, F=None)
        >>>
        >>> #with custom features
        >>> X_feat = from_matrix_v3(X, alphabet=aas, F=varimax_matrix, reshape=True)
    """
    try:
        alphabet = np.asarray(alphabet)
    except:
        raise ValueError(
            '<featurize.from_matrix_v3> received an invalid "alphabet" variable that could not be coerced into an array.'
        )

    # special treatment for padded arrays: try to rescue featurization
    # the problem is that the feature matrix will probably not
    # have a row to encode pads
    _is_padded = X.dtype.type() in X
    if _is_padded:
        # just in case, let's assert that this is true
        assert X.dtype.type() == alphabet.dtype.type(), "<from_matrix_v3>: X and alphabet dtypes must match!"

        if alphabet.dtype.type() not in alphabet:
            alphabet = np.hstack((alphabet.dtype.type(), alphabet))

        if F is not None:
            if F.shape[0] != alphabet.size:
                if F.shape[0] == alphabet.size - 1:
                    F = np.vstack((np.zeros(F.shape[1]), F))
                else:
                    raise ValueError(
                        "<featurize.from_matrix_v3>: F matrix and X array have incompatible dimensions; cannot perform featurization. . ."
                    )

    x_shape = X.shape
    X = X.ravel()
    expansion = alphabet.size

    # convert the matrix to a one-hot encoding
    sort_idx = np.argsort(alphabet)
    X = np.searchsorted(alphabet, X, sorter=sort_idx)
    fX = np.zeros((X.size, expansion))
    fX[np.arange(X.size), X] = 1

    # matmul by a factorization matrix to get the featurized repr
    if F is not None:
        fX = np.matmul(fX, F)

    fX = np.reshape(fX, (x_shape[0], -1))
    if reshape:
        fX = np.reshape(fX, x_shape + (-1,))

    return fX


def RFA_featurization(X, alphabet=None, order=None):
    """
    Featurize sequences for Reference-Free Analysis models.

    Generates first or second-order positional features for RFA models.
    See Nat Commun 2024, 15, 7953 for methodology details.

    Args:
        X (ndarray): Sequences to featurize (2D array).

        alphabet (array-like): Monomer alphabet.

        order (str): 'first' for first-order features or 'second' for
            first + second-order features.

    Returns:
        ndarray: Flattened feature array (n_sequences, n_features).
    """

    def _first_order_features(X, alphabet):
        x_shape = X.shape
        X = X.ravel()
        expansion = alphabet.size

        # convert the matrix to a one-hot encoding
        sort_idx = np.argsort(alphabet)
        X = np.searchsorted(alphabet, X, sorter=sort_idx)
        fX = np.zeros((X.size, expansion))
        fX[np.arange(X.size), X] = 1

        fX = fX.reshape(x_shape + (expansion,))
        return np.reshape(fX, (x_shape[0], -1)).astype(np.float32)

    def _second_order_features(fX, alphabet):
        # fX: first order features

        # we first need to reshape the one-hot array back to a suitable form
        fX = fX.reshape((fX.shape[0], -1, alphabet.size))

        # broadcasting is amazing
        sX = fX[:, :, None, :, None] * fX[:, None, :, None, :]

        # return only non-redundant enties,
        # i. e., the upper triangle of the hypercube
        i, j = np.triu_indices(sX.shape[1], k=1, m=sX.shape[1])
        return sX[:, i, j].reshape((sX.shape[0], -1))

    if order != "first" and order != "second":
        raise ValueError(
            '<RFA_featurization>: the "order" keyword not understood: only "first" or "second" are allowed values!'
        )

    try:
        alphabet = np.asarray(alphabet)
    except:
        raise ValueError(
            '<featurize.from_matrix_v3> received an invalid "alphabet" variable that could not be coerced into an array.'
        )

    # compute the first order features
    fX = _first_order_features(X, alphabet)
    if order == "second":
        sX = _second_order_features(fX, alphabet)
        return np.hstack([fX, sX]).astype(np.float32)

    return fX.astype(np.float32)


def into_h5(X, y=None, path=None, F=None, alphabet=None, reshape=False, chunks=20):
    """
    Featurize sequences and save to HDF5 file.

    Memory-efficient featurization for large datasets. Processes data in
    chunks and writes directly to disk.

    Args:
        X (ndarray): Sequences to featurize (2D array).

        y (ndarray, optional): Labels (1D array).

        path (str): Output HDF5 file path.

        F (ndarray, optional): Feature matrix. If None, uses one-hot encoding.

        alphabet (array-like): Token alphabet.

        reshape (bool): If True, keeps 2D sequence structure. If False,
            flattens to vectors. Default is False.

        chunks (int): Number of chunks to split X into for processing.
            Default is 20.

    Example:
        >>> into_h5(X, y=labels, path='features.h5', F=None,
        ...         alphabet=aas, reshape=False, chunks=10)
    """
    try:
        import h5py   

    except ImportError:
        msg = "Failed to import h5py. Please install clibas with `pip install clibas[ml]`. . ."
        raise ImportError(msg)

    with h5py.File(path, "w") as f:
        # chunks is the number of subarrays.
        indx = np.array_split(np.arange(X.shape[0]), chunks)

        # get the first chunk before iterating -> this allows
        # the creation of h5 datasets
        x = from_matrix_v3(
            X[indx[0]],
            F=F,
            alphabet=alphabet,
            reshape=reshape,
        )

        dims = (X.shape[0],) + x.shape[1:]
        x_set = f.create_dataset("X", dims, dtype=np.float32)

        # write this first chunk
        x_set[indx[0][0] : indx[0][-1] + 1] = x

        if len(indx) > 1:
            for ind in indx[1:]:
                x_set[ind[0] : ind[-1] + 1] = from_matrix_v3(
                    X[ind],
                    F=F,
                    alphabet=alphabet,
                    reshape=reshape,
                )
        if y.ndim > 0:
            y_set = f.create_dataset("y", (y.size,), dtype=np.float32)
            y_set[...] = y

    return
