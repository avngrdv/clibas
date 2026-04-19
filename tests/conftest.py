"""
Created on Sat Oct 18 01:38:54 2025
@author: Alex Vinogradov
"""

from pathlib import Path

import pytest


@pytest.fixture
def basic_valid_config():
    return Path(__file__).parent / "data" / "basic_valid_config.yaml"


@pytest.fixture
def basic_invalid_config():
    return Path(__file__).parent / "data" / "basic_invalid_config.yaml"


@pytest.fixture
def Walport_config():
    return Path(__file__).parent / "data" / "Walport_PADI4_config.yaml"


@pytest.fixture
def Heinis_config():
    return Path(__file__).parent / "data" / "Heinis_FXIa_config.yaml"


@pytest.fixture
def data_dir():
    return Path(__file__).parent / "data"


@pytest.fixture
def Heinis_dataset():
    import numpy as np

    from clibas.datatypes import AnalysisSample, Data

    f1 = Path(__file__).parent / "data" / "FXIa_r3_L7_pep.npy"
    f2 = Path(__file__).parent / "data" / "FXIa_r6_L7_pep.npy"

    data = Data(
        [
            AnalysisSample(X=np.load(f1), name="FXIa_r3_L7_pep"),
            AnalysisSample(X=np.load(f2), name="FXIa_r6_L7_pep"),
        ]
    )

    return data


@pytest.fixture
def Walport_dataset():
    import numpy as np

    from clibas.datatypes import AnalysisSample, Data

    f1 = Path(__file__).parent / "data" / "PADI4_r5_15_min_test_dna.npy"
    f2 = Path(__file__).parent / "data" / "PADI4_r6_15_min_test_dna.npy"

    data = Data(
        [
            AnalysisSample(X=np.load(f1), name="PADI4_r5_15_min_test_dna"),
            AnalysisSample(X=np.load(f2), name="PADI4_r6_15_min_test_dna"),
        ]
    )

    return data

@pytest.fixture
def basic_valid_peptides():
    
    return [
                'AYATAKA',
                'HEEISGYRY',
                'MANADRAIN',
                'AYATAKA',
                'VCHERAVECHERYM',
                'QRIIVA',
                'AYATAKA',
                'MANKQSTWIN'        
    ]

@pytest.fixture
def basic_invalid_peptides():
    
    return [
                'AYATAKa',
                'HEEISGYRY',
                'MANADRAIN',
                'AYATAKA',
                'VCHERAVECHERYM',
                'QRI-VA',
                'AYATAKA',
                '2ANKQSTWIN'        
    ]
























