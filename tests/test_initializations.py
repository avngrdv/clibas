"""
Created on Sat Oct 18 01:38:54 2025
@author: Alex Vinogradov
"""

import pytest


def test_valid_initialization(basic_valid_config):
    import clibas as C

    C.initialize(basic_valid_config)
    assert C.fastq_parser


def test_invalid_initialization(basic_invalid_config):
    import clibas as C

    with pytest.raises(Exception):
        C.initialize(basic_invalid_config)
