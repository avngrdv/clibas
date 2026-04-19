"""
Utility functions for sequence analysis and statistics.

Provides functions for frequency calculations, entropy metrics, Hamming
distance, clustering scores, and other common sequence analysis operations.
"""

import numpy as np


def get_freqs(arr, alphabet):
    """
    Compute positional token frequencies.

    Args:
        arr (ndarray): Dataset (2D array).
        alphabet (array-like): Tokens to count. Must match arr.dtype.

    Returns:
        ndarray: Frequency matrix (n_tokens, sequence_length).
    """
    assert len(set(alphabet)) == len(alphabet), (
        "Token alphabet should not \
contain duplicated tokens!"
    )

    # C: count matrix for tokens over positions in arr
    C = np.zeros((len(alphabet), arr.shape[1]))

    # iteratively fill it
    for i, x in enumerate(alphabet):
        C[i] = np.sum(arr == x, axis=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        freq = np.divide(C, arr.shape[0])

    return freq


def get_Y_star(f1, f2, alphabet, f_out=None):
    """
    Compute enrichment Y scores from positive/negative datasets.

    Args:
        f1 (str): Path to positive dataset (.npy file).

        f2 (str): Path to negative dataset (.npy file).

        alphabet: Token alphabet.

        f_out (str, optional): Output path for Y scores.

    Returns:
        ndarray: Log-ratio enrichment scores (n_tokens, n_positions).
    """
    # load positive and negative P matrices
    pos = np.load(f1).astype(str)
    neg = np.load(f2).astype(str)

    freq_pos = get_freqs(pos, alphabet)
    freq_neg = get_freqs(neg, alphabet)

    # calculate Y matrix from it and save it
    Y = np.log(np.divide(freq_pos, freq_neg))
    if f_out is not None:
        np.save(f_out, Y)

    return Y


def positional_conservation(freq):
    """
    Compute position-wise sequence conservation from a frequency
    matrix array.

    Args:
        freq (ndarray): Frequency matrix from get_freqs().

    Returns:
        ndarray: Conservation values per position (1D array).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        E = np.nan_to_num(np.multiply(freq, np.log2(freq)))

    n = np.log2(freq.shape[0])
    return np.divide(np.sum(E, axis=0) + n, n)


def arr_purity(arr, alphabet):
    """
    Compute average positional conservation (purity metric) without
    doing multiple sequence alignment; the dataset is used as is for
    the computation.

    Args:
        arr (ndarray): Dataset (2D array).

        alphabet (array-like): Token alphabet.

    Returns:
        float: Average conservation across all positions (0-1 scale).
    """
    freq = get_freqs(arr, alphabet)
    conservation = positional_conservation(freq)
    return np.divide(np.sum(conservation), arr.shape[1])


def naive_clustering_score(X, labels, alphabet=None, return_mean=False):
    """
    Evaluate clustering quality based on within-cluster conservation.

    Args:
        X (ndarray): Sequences (2D array).

        labels (ndarray): Cluster assignments (1D array, -1 for noise).

        alphabet (array-like): Token alphabet.

        return_mean (bool): If True, returns mean score. If False, returns
            per-cluster scores. Default is False.

    Returns:
        float or ndarray: Clustering score(s).
    """
    starting_purity = arr_purity(X, alphabet)
    n_labels = np.unique(labels)

    cluster_purities = np.zeros(n_labels.size)
    cluster_sizes = np.zeros(n_labels.size)

    # compute purities for every cluster
    for i, label in enumerate(n_labels):
        # everything labelled as noise, receives the score of 0 automatically
        if label == -1:
            cluster_purities[i] = starting_purity
            cluster_sizes[i] = X[labels == label].shape[0]

        else:
            cluster_purities[i] = arr_purity(X[labels == label], alphabet)
            cluster_sizes[i] = X[labels == label].shape[0]

    # compute individual cluster scores
    weighed_added_purity = np.multiply(
        cluster_purities - starting_purity, np.log2(cluster_sizes)
    )

    if return_mean:
        # the final score is the average of the individual scores
        return weighed_added_purity.mean()
    else:
        return weighed_added_purity


def hamming_distance(
    P,
    pep,
    h=0,
    cum=False,
    return_count=False,
    return_index=False,
    return_distance=False,
):
    """
    Calculate Hamming distances between sequences.

    Flexible function for Hamming distance calculations with various output
    options.

    Args:
        P (ndarray): Sequences to compare (2D array).

        pep (ndarray): Reference sequence (1D array).

        h (int): Hamming distance threshold. Default is 0.

        cum (bool): If True, returns sequences at distance ≤ h. If False,
            returns sequences at exactly distance h. Default is False.

        return_count (bool): If True, returns count instead of sequences.

        return_index (bool): If True, returns indices instead of sequences.

        return_distance (bool): If True, returns array of distances.

    Returns:
        ndarray: Filtered sequences, indices, count, or distances depending
            on parameters.
    """
    D = pep == P

    if return_distance:
        return np.sum(~D, axis=1)

    match = pep.size - h
    if cum:
        ind = np.sum(D, axis=1) >= match
    else:
        ind = np.sum(D, axis=1) == match

    H = P[ind]

    if return_count:
        return H.shape[0]

    elif return_index:
        return np.where(ind)[0]

    return H


def shannon_entropy(arr, norm=True, return_counts=True):
    """
    Compute Shannon entropy of a dataset.

    Args:
        arr (ndarray): Data array (any representation).

        norm (bool): If True, returns normalized entropy (efficiency).
            If False, returns raw entropy. Default is True.

        return_counts (bool): If True, also returns unique entry counts.

    Returns:
        float or tuple: Entropy value, or (entropy, counts) if return_counts=True.
    """
    # C - counts; n - dataset size
    C = np.unique(arr, return_counts=True, axis=0)[1]
    n = C.sum()
    normC = np.divide(C, n)

    # E - a vector of individual entropy values
    E = -normC * np.log2(normC)
    if norm == True:
        E = np.divide(E.sum(), np.log2(n))

    else:
        E = E.sum()

    if return_counts:
        return E, C
    else:
        return E


def sample_random_sequences(n, y, monomers):
    """
    Generate random (peptide/DNA) sequences .

    Args:
        n (int): Number of sequences to generate.

        y (int): Sequence length.

        monomers (list): Tokens (monomers) to sample from.

    Returns:
        ndarray: Random sequences (n, y).
    """
    P = np.random.choice(monomers, size=(n, y), replace=True)
    return P


def sample_from_template(template, n, monomers):
    """
    Generate partially randomized sequences from template.

    Args:
        template (ndarray): Template sequence (1D, dtype='<U1').
            'X' denotes positions to randomize.

        n (int): Number of sequences to generate.

        monomers (list): Monomers to sample for 'X' positions.

    Returns:
        ndarray: Partially randomized sequences (n, template_length).
    """

    P = sample_random_sequences(n, len(template), monomers)
    for i, aa in enumerate(template):
        if aa != "X":
            P[:, i] = [aa] * P.shape[0]

    return P


def sample_from_template_v2(template, n, monomers):
    """Generate sequences from template with numeral-based specification.

    Args:
        template (ndarray): Template with numerals indicating monomer sets.

        n (int): Number of sequences to generate.

        monomers (dict): Mapping of numerals to monomer lists.

    Returns:
        ndarray: Generated sequences (n, template_length).
    """
    P = np.zeros((n, template.size), dtype="<U1")

    for pos in range(template.size):
        if not template[pos].isdigit():
            P[:, pos] = template[pos]

        else:
            P[:, pos] = np.random.choice(monomers[template[pos]], size=n, replace=True)
    return P


def sorted_count(arr, top_n=None, return_index=False):
    """
    Count and sort unique entries by abundance.

    Args:
        arr (ndarray): Data array (any dimensionality).

        top_n (int, optional): Return only top N most abundant entries.

        return_index (bool): If True, also return original indices.

    Returns:
        tuple: (unique_entries, counts) or (unique_entries, indices, counts).
    """
    if arr.ndim > 1:
        X, og_ind, C = np.unique(arr, return_counts=True, return_index=True, axis=0)
    else:
        X, og_ind, C = np.unique(arr, return_counts=True, return_index=True)

    ind = np.argsort(C)[::-1][:top_n]
    if return_index:
        return (X[ind], og_ind[ind], C[ind])

    return (X[ind], C[ind])


def get_S(P, Y):
    """
    Compute enrichment S scores for sequences.

    Args:
        P (ndarray): Sequences (2D array, dtype=int).

        Y (ndarray): Enrichment score matrix from get_Y_star().

    Returns:
        ndarray: S scores (sum of Y scores per sequence).
    """

    def F(pep):
        return sum([Y[x, i] for i, x in enumerate(pep)])

    return np.array(list(map(F, P)))


def logistic_4_param(x, A, B, C, D):
    """
    Evaluate 4-parameter logistic curve.

    y = A + B / (1 + exp(-C*(x-D)))

    Args:
        x: Input value(s).
        A: Minimum asymptote.
        B: Curve span.
        C: Hill slope.
        D: Inflection point.

    Returns:
        Logistic curve value(s).
    """
    exp = np.exp(np.multiply(-C, x - D))
    return A + np.divide(B, 1 + exp)


def freq_to_information_content(F):
    """Convert frequency matrix to information content for sequence logos.

    Args:
        F (DataFrame): Frequency matrix (positions × tokens).

    Returns:
        DataFrame: Information content matrix for logo plotting.
    """
    IC = F.copy()
    n_alphabet = F.shape[1]

    for i in F.index:
        freqs = F.loc[i].values
        entropy = -np.sum([p * np.log2(p) for p in freqs if p > 0])
        R = np.log2(n_alphabet) - entropy
        IC.loc[i] = freqs * float(R)

    return IC
