"""
Created on Sat Oct 18 01:38:54 2025
@author: Alex Vinogradov
"""

import pytest


@pytest.mark.slow
def test_umap_Heinis(Heinis_config, tmp_path, Heinis_dataset):
    import clibas as C

    C.initialize(Heinis_config)

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "outputs")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    umap = C.analysis_tools.umap_hdbscan_analysis(
        top_n=1000,
        where="pep",
        F="pep_ECFP4",
        single_manifold=True,
        return_modified=True,
    )
    data = umap(Heinis_dataset)
    assert hasattr(data[0], "labels")


@pytest.mark.slow
def test_umap_Walport(Walport_config, tmp_path, Walport_dataset):
    import clibas as C

    C.initialize(Walport_config)

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "outputs")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    umap = C.analysis_tools.umap_hdbscan_analysis(
        top_n=400, where="dna", single_manifold=True, return_modified=True
    )
    data = umap(Walport_dataset)
    assert hasattr(data[0], "C")
    assert hasattr(data[0], "labels")
    assert data.size == 2
    return


def test_umap_basics(tmp_path, basic_valid_peptides, basic_invalid_peptides):
    import numpy as np
    import pandas as pd

    import clibas as C

    C.initialize()

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "outputs")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    # where: specified, alphabet is inferred
    umap = C.analysis_tools.umap_hdbscan_analysis(where="pep", return_modified=True)

    # take lists as X inputs
    data = umap(basic_valid_peptides)
    assert data[0].size == 6
    assert hasattr(data[0], "C")
    assert hasattr(data[0], "labels")

    # where: NOT specified, valid alphabet is passed explicitly
    aas = [x.decode("ascii") for x in C.analysis_tools.constants.aas]
    umap = C.analysis_tools.umap_hdbscan_analysis(return_modified=True, alphabet=aas)

    # take nd arrays as X inputs
    data = umap(np.array(basic_valid_peptides))
    assert data[0].size == 6
    assert hasattr(data[0], "C")
    assert hasattr(data[0], "labels")

    # take pandas series as X inputs
    data = umap(pd.Series(basic_valid_peptides))
    assert data[0].size == 6
    assert hasattr(data[0], "C")
    assert hasattr(data[0], "labels")

    # where: NOT specified, INvalid alphabet is passed explicitly
    with pytest.raises(ValueError):
        aa2d = np.reshape(aas, shape=(10, 2))
        umap = C.analysis_tools.umap_hdbscan_analysis(
            return_modified=True, alphabet=aa2d
        )

    # invalid peptide lists should raise
    with pytest.raises(Exception):
        umap = C.analysis_tools.umap_hdbscan_analysis(
            return_modified=True, alphabet=aas
        )
        data = umap(basic_invalid_peptides)

    return
