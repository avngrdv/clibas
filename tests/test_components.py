"""
Created on Sat Oct 18 01:57:47 2025
@author: Alex Vinogradov
"""

import pytest


def test_feature_matrix(basic_valid_config):
    
    import clibas as C
    C.initialize(basic_valid_config)

    from clibas.featurization import FeatureMatrix

    F = FeatureMatrix.make(descr="pep_ECFP4", constants=C.fastq_parser.constants)
    F.single_linkage_tree(labels=F.constants.aas)
    F.single_linkage_tree(labels=None)
    assert F.F.shape == (20, 208)

    F = FeatureMatrix.make(descr="varimax", constants=C.fastq_parser.constants)
    assert 1.185 > F.F.sum() > 1.179 

    F = FeatureMatrix.make(descr="pep_ECFP3", constants=C.fastq_parser.constants)
    assert F.F.shape == (20, 195)
    
    # this is actually deprecated but let's test it
    F = FeatureMatrix.make(descr="pep_SMILES", constants=C.fastq_parser.constants)
    assert F.F.shape == (20, 31)    

    F = FeatureMatrix.make(descr=None, constants=C.fastq_parser.constants)
    assert F.F is None


invalid_designs = [
    # topology not the same
    (
        [
            b"GGGGCG112112TGC112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA123",
        ],
        {1: (b"A", b"G", b"T", b"C"), 2: (b"G", b"T")},
        "dna",
        (b"A", b"C", b"G", b"T"),
    ),
    # fails validation (Q monomer)
    (
        [
            b"GGGGCG112112TGC112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATQCCCA",
        ],
        {1: (b"A", b"G", b"T", b"C"), 2: (b"G", b"T")},
        "dna",
        (b"A", b"C", b"G", b"T"),
    ),
    # fails validation v2 (Q monomer)
    (
        [
            b"GGGGCG112112TGC112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA",
        ],
        {1: (b"A", b"G", b"T", b"C"), 2: (b"Q", b"T")},
        "dna",
        (b"A", b"C", b"G", b"T"),
    ),
    # not a byte-string
    (
        [
            "GGGGCG112112TGC112112GGAGGATACCCA",
            "GGGGCG112112TGC112112112GGAGGATACCCA",
            "GGGGCG112112TGC112112112GGAGGATACCCA",
        ],
        {1: (b"A", b"G", b"T", b"C"), 2: (b"G", b"T")},
        "dna",
        (b"A", b"C", b"G", b"T"),
    ),
    # not a byte-string v2
    (
        [
            b"GGGGCG112112TGC112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA123",
        ],
        {1: (b"A", b"G", b"T", b"C"), 2: ("G", "T")},
        "dna",
        (b"A", b"C", b"G", b"T"),
    ),
    # not all random sets specified
    (
        [
            b"GGGGCG112112TGC112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA",
        ],
        {1: (b"A", b"G", b"T", b"C")},
        "dna",
        (b"A", b"C", b"G", b"T"),
    ),
]

@pytest.mark.parametrize("templates,monomers,lib_type,val_monomer", invalid_designs)
def test_invalid_lib_designs(templates, monomers, lib_type, val_monomer):
    from clibas.lib_design import LibraryDesign

    with pytest.raises(Exception):
        LibraryDesign(
            templates=templates,
            monomers=monomers,
            lib_type=lib_type,
            val_monomer=val_monomer,
        )


def test_valid_lib_design():
    from clibas.lib_design import LibraryDesign

    LibraryDesign(
        templates=[
            b"GGGGCG112112TGC112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA",
            b"GGGGCG112112TGC112112112GGAGGATACCCA",
        ],
        monomers={1: (b"A", b"G", b"T", b"C"), 2: (b"G", b"T")},
        lib_type="dna",
        val_monomer=(b"A", b"C", b"G", b"T"),
    )


def test_array_casting_edge_cases(basic_valid_peptides):
    import clibas as C
    import numpy as np
    import pandas as pd
    C.initialize()
    
    from clibas.datatypes import SequencingSample, AnalysisSample
    
    # empty array
    empty = SequencingSample(dna=np.array([]), Q=np.array([]))
    assert empty.size == 0
    
    # single sequence
    single = SequencingSample(
        dna=np.array([b'ACGT']), 
        Q=np.array([b'IIII'])
    )
    assert single.size == 1
    
    # mismatched dimensions (should fail)
    with pytest.raises(AssertionError):
        SequencingSample(dna=np.array([b'ACGT', b'GGGG']), Q=np.array([b'IIII']))

    # empty array: AnalysisSample
    empty = AnalysisSample(X=[], name='!!!')          
    assert empty.size == 0

    # basic from list
    basic = AnalysisSample(X=basic_valid_peptides, name='peptides')          
    assert basic.size == 8
    
    # basic from nd array
    basic = AnalysisSample(X=np.asarray(basic_valid_peptides), name='peptides')          
    assert basic.size == 8    
    
    # basic from pd Series
    basic = AnalysisSample(X=pd.Series(basic_valid_peptides), name='peptides')          
    assert basic.size == 8        

    # unacceptable name
    with pytest.raises(ValueError):
        AnalysisSample(X=basic_valid_peptides, name=basic_valid_peptides)    

    return

def test_translation_edge_cases():
    
    import clibas as C
    C.initialize()
    
    parser = C.fastq_parser
    
    # incomplete codon at end
    assert parser._dna_to_pep_v2(b'ATGAAA', force_at_frame=0) == b'MK'
    assert parser._dna_to_pep_v2(b'ATGAAAT', force_at_frame=0) == b'MK_'
    
    # ambiguous bases (N)
    result = parser._dna_to_pep_v2(b'ATGNNN', force_at_frame=0)
    assert b'+' in result
    
    # stop codon without readthrough
    result = parser._dna_to_pep_v2(b'ATGTAAGGG', force_at_frame=0, stop_readthrough=False)
    assert result == b'M'  # stops at TAA
    
    # stop codon WITH readthrough
    result = parser._dna_to_pep_v2(b'ATGTAAGGG', force_at_frame=0, stop_readthrough=True)
    assert result == b'M*G'  # continues past stop
    
    # different frames
    seq = b'ATACTGAAAGGG'
    assert parser._dna_to_pep_v2(seq, force_at_frame=0) == b'ILKG'
    assert parser._dna_to_pep_v2(seq, force_at_frame=1) == b'Y'
    assert parser._dna_to_pep_v2(seq, force_at_frame=2) == b'TER_'
    
    # empty sequence
    assert parser._dna_to_pep_v2(b'', force_at_frame=0) == b''
    return


def test_internal_state_collapse():
    import clibas as C
    import numpy as np
    C.initialize()
    
    from clibas.datatypes import SequencingSample
    
    # create sample with multi-template matches
    sample = SequencingSample(
        dna=np.array([b'ACGT', b'GGGG', b'TTTT']),
          Q=np.array([b'IIII', b'IIII', b'IIII']),
    )
    
    # simulate multiple template matches
    sample._internal_state = np.array([
        [True, True, False],   # read 0 matches templates 0,1
        [False, True, True],   # read 1 matches templates 1,2
        [True, False, False],  # read 2 matches template 0 only
    ])
    
    assert not sample._is_collapsed()
    sample._collapse_internal_state()
    assert sample._is_collapsed()
    
    # after collapse, each read matches exactly one template
    assert np.sum(sample._internal_state, axis=1).tolist() == [1, 1, 1]
    
    # first matching template should be selected (argmax behavior)
    assert sample._internal_state[0, 0] == True  # selected template 0
    assert sample._internal_state[1, 1] == True  # selected template 1


def test_data_stacking_different_widths():
    import clibas as C
    import numpy as np
    C.initialize()
    
    from clibas.datatypes import AnalysisSample, Data
    
    # different width sequences
    s1 = AnalysisSample(X=np.array([[b'A', b'C'], [b'G', b'T']]), name='s1')
    s2 = AnalysisSample(X=np.array([[b'A', b'C', b'G', b'T']]), name='s2')
    
    data = Data([s1, s2])
    stacked = data.stack(in_place=False)
    
    assert stacked.size == 1
    assert stacked[0].X.shape == (3, 4)
    # first two rows should have pads in columns 2,3
    assert stacked[0].X[0, 2] == b''
    assert stacked[0].X[1, 3] == b''
    assert stacked[0].y.shape == ()


def test_filtering_preserves_internal_state():
    import clibas as C
    import numpy as np
    C.initialize()
    
    from clibas.datatypes import SequencingSample
    
    sample = SequencingSample(
        dna=np.array([b'ACGT', b'GGGG', b'TTTT']),
        Q=np.array([b'IIII', b'IIII', b'IIII']),
        pep=np.array([b'AC', b'GG', b'TT'])
    )
    
    sample._internal_state = np.ones((3, 2), dtype=bool)
    
    # filter out middle sequence
    sample.ind_filter(np.array([True, False, True]))
    
    assert sample.size == 2
    assert sample._internal_state.shape == (2, 2)
    assert sample.dna.shape[0] == 2
    assert sample.pep.shape[0] == 2


def test_q_score_operations():
    import clibas as C
    import numpy as np
    C.initialize()
    
    from clibas.datatypes import SequencingSample
    
    # Phred+33 encoding: '!' = 0, 'I' = 40
    sample = SequencingSample(
        dna=np.array([b'ACGT']),
          Q=np.array([b'!II~'])  # low, high, very high
    )
    
    civilized = sample.civilized_Q
    assert civilized[0, 0] == 0    # '!' - 33 = 0
    assert civilized[0, 1] == 40   # 'I' - 33 = 40
    assert civilized[0, 2] == 40
    assert civilized[0, 3] == 93   # '~' - 33 = 93
    
    # padding should remain 0
    sample2 = SequencingSample(
        dna=np.array([b'AC', b'ACGT']),
          Q=np.array([b'II', b'IIII'])
    )
    
    civil2 = sample2.civilized_Q
    assert civil2[0, 2] == 0  # padding byte
    assert civil2[0, 3] == 0  # padding byte


def test_hamming_distance():
    import numpy as np
    from clibas.misc import hamming_distance

    ref = np.array([b'A', b'C', b'G', b'T'])
    
    seqs = np.array([
        [b'A', b'C', b'G', b'T'],
        [b'A', b'C', b'G', b'G'],  # HD=1 from ref
        [b'A', b'T', b'G', b'T'],  # HD=1 from ref
        [b'T', b'T', b'T', b'T'],  # HD=3 from ref
    ])

    # exact match (h=0)
    result = hamming_distance(seqs, ref, h=0, cum=False)
    assert result.shape[0] == 1  # only first sequence
    
    # within distance 1 (cumulative)
    result = hamming_distance(seqs, ref, h=1, cum=True)
    assert result.shape[0] == 3  # first three sequences
    
    # return distances
    distances = hamming_distance(seqs, ref, return_distance=True)
    assert distances.tolist() == [0, 1, 1, 3]
    
    # return indices
    indices = hamming_distance(seqs, ref, h=1, cum=True, return_index=True)
    assert indices.tolist() == [0, 1, 2]


def test_matrix_featurization_one_hot():
    import numpy as np
    from clibas.featurization import from_matrix_v3
    
    X = np.array([[b'A', b'C'], [b'G', b'T']])
    alphabet = np.array([b'A', b'C', b'G', b'T'])
    
    # should work
    result = from_matrix_v3(X, alphabet=alphabet, F=None)
    assert result.shape == (2, 2, 4)  # one-hot
    
    # alphabet too small (missing 'T')
    bad_alphabet = np.array([b'A', b'C', b'G'])
    with pytest.raises(Exception):
        from_matrix_v3(X, alphabet=bad_alphabet, F=None)

    X = np.array([[b'A', b'C', b''], [b'G', b'T', b''], [b'G', b'T', b'T']])
    result = from_matrix_v3(X, alphabet=alphabet, F=None)
    assert result.shape == (3, 3, 5)  # one-hot    
    assert result[:,:,0].sum() == 2


def test_matrix_featurization_from_F(basic_valid_peptides):
    
    import numpy as np
    from clibas.featurization import from_matrix_v3
    from clibas.featurization import FeatureMatrix
    from clibas.datatypes import AnalysisSample
    import clibas as C
    C.initialize()
    
    F = FeatureMatrix.make(
        descr='pep_ECFP4', constants=C.analysis_tools.constants
    ).F
    
    X = AnalysisSample(X=basic_valid_peptides).X
    
    alphabet = C.analysis_tools._infer_alphabet(alphabet='aa')
    result = from_matrix_v3(X, F=F, alphabet=alphabet)
    assert result.shape == (8, 14, 208)
    assert result.sum() == 1408 
    assert np.all(result[0][7:, :] == 0) #aa pads should be represented as zeros
    assert np.all(result[0,6]  == result[0,0]) #one aa A has a consistent repr

def test_translation_performance():
    
    import time
    import clibas as C
    import numpy as np
    C.initialize()
    
    # generate 10k sequences
    from clibas.datatypes import SequencingSample, Data
    n = 10000
    dna = np.array([b'ATGAAAGGGCCCTTT'] * n)
    Q = np.array([b'I' * 15] * n)
    sample = SequencingSample(dna=dna, Q=Q, pep=None)
    data = Data([sample])
    
    translator = C.fastq_parser.translate()
    
    start = time.time()
    data = translator(data)
    elapsed = time.time() - start
    
    # should translate 10k sequences in under 1 second
    assert elapsed < 1.0, f"Translation took {elapsed:.2f}s, expected <1s"
    assert data[0].pep.shape[0] == n












