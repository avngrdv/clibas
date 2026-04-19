"""
Created on Sun Oct 26 01:08:15 2025
@author: Alex Vinogradov
"""

import pytest


def test_L_basics(tmp_path, basic_valid_peptides):
    import clibas as C

    C.initialize()

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "outputs")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    L = C.analysis_tools.length_analysis()

    data = L(tuple(basic_valid_peptides))
    assert data == tuple(basic_valid_peptides)

    out_dir = tmp_path / "outputs" / "unnamed"
    f = any(out_dir.glob("unnamed_None_L_distribution*"))
    assert f


def test_token_convergence_basics(tmp_path, basic_valid_peptides):
    import numpy as np

    import clibas as C

    C.initialize()

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "outputs")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    # valid alphabet initialization in the absence of 'where'
    aas = [x.decode("ascii") for x in C.analysis_tools.constants.aas]
    T = C.analysis_tools.token_convergence_analysis(alphabet=aas)

    data = T(basic_valid_peptides)
    assert data == basic_valid_peptides

    out_dir = tmp_path / "outputs" / "unnamed"
    f = any(out_dir.glob("unnamed_None_reg_overall_seq_logo*"))
    assert f

    f = any(out_dir.glob("unnamed_None_reg_overall_tokenwise_frequency*"))
    assert f

    # valid alphabet initialization with where
    T = C.analysis_tools.token_convergence_analysis(where="pep")
    data = T(basic_valid_peptides)
    assert data == basic_valid_peptides

    # invalid alphabet should raise
    with pytest.raises(ValueError):
        aa2d = np.reshape(aas, shape=(10, 2))
        T = C.analysis_tools.token_convergence_analysis(alphabet=aa2d)

    return


def test_sequence_convergence_basics(
    tmp_path, basic_valid_peptides, basic_invalid_peptides
):
    import numpy as np

    import clibas as C

    C.initialize()

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "outputs")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    # valid alphabet initialization in the absence of 'where'
    S = C.analysis_tools.sequence_convergence_analysis()

    data = S(basic_valid_peptides)
    assert data == basic_valid_peptides

    out_dir = tmp_path / "outputs" / "unnamed"
    f = any(out_dir.glob("unnamed_None_library_convergence*"))
    assert f

    # this should also work even in the presence of 'where'
    S = C.analysis_tools.sequence_convergence_analysis(where="pep")

    data = S(np.array(basic_valid_peptides))
    assert np.all(data == np.array(basic_valid_peptides))
    return
