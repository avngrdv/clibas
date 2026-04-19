"""
Created on Fri Oct 24 23:41:21 2025
@author: Alex Vinogradov
"""

from pathlib import Path

import pytest

valid_configs = [
    "valid_minimum_config.yaml",
    "valid_no_aa_SMILES_config_non_standard.yaml",
    "valid_no_codon_table_config.yaml",
    "valid_no_complement_config.yaml",
    "valid_no_const_config.yaml",
    "valid_no_dirs_config.yaml",
    "valid_no_lib_design_config.yaml",
    "valid_no_logger_config.yaml",
    "valid_no_utr_config.yaml",
    "valid__no_aa_SMILES_config.yaml",
    None,
]

invalid_configs = [
    "invalid_faulty_aas_config.yaml",
    "invalid_faulty_SMILES_config.yaml",
    "invalid_no_dna_monomers_config.yaml",
    "invalid_no_pep_design_config.yaml",
]


@pytest.mark.parametrize("valid_config_fname", valid_configs)
def test_valid_initialization(valid_config_fname):
    import clibas as C

    if valid_config_fname is not None:
        p = Path(__file__).parent / "data" / valid_config_fname
    else:
        p = None

    C.initialize(config_path=p)
    assert hasattr(C.fastq_parser, "constants")
    assert hasattr(C.fastq_parser, "logger")
    assert hasattr(C.fastq_parser, "dirs")
    assert hasattr(C.fastq_parser.constants, "aas")
    assert hasattr(C.fastq_parser.constants, "bases")
    assert hasattr(C.fastq_parser.constants, "translation_table")
    assert hasattr(C.fastq_parser.constants, "complement_table")
    # SMILES, lib designs may or may not be present


@pytest.mark.parametrize("invalid_config_fname", invalid_configs)
def test_invalid_initialization(invalid_config_fname):
    import clibas as C

    p = Path(__file__).parent / "data" / invalid_config_fname
    with pytest.raises(Exception):
        C.initialize(config_path=p)
