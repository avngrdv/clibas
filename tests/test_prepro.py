"""
Created on Sat Oct 18 01:38:54 2025
@author: Alex Vinogradov
"""

import pytest


def test_RFA_featurization_Heinis(Heinis_config, tmp_path, Heinis_dataset):
    import clibas as C

    C.initialize(Heinis_config)

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.ml_data = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    to_pop = [
        "ACCHMRSVCYRCP",
        "VCDCCYRCP",
        "KCFDCCYRCP",
        "YCFSCCYRCP",
        "ECSHCERSMCDRCP",
    ]

    C.pipeline.enque(
        [
            C.preprocessor.intrasample_unique(),
            C.preprocessor.intersample_unique(),
            C.preprocessor.merge(),
            C.preprocessor.intrasample_unique(),
            C.preprocessor.filter_external(external=to_pop, max_hd=0),
            C.preprocessor.shuffle(),
            C.preprocessor.featurize_for_RFA(alphabet="aa", order="second"),
        ]
    )

    data = C.pipeline.run(Heinis_dataset)
    assert data[0].X.shape == (817, 36680)
    assert data[0].X.sum() == 85785


def test_ECFP_featurization_Heinis(Heinis_config, tmp_path, Heinis_dataset):
    import clibas as C

    C.initialize(Heinis_config)

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.ml_data = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    C.pipeline.enque(
        [
            C.preprocessor.merge(),
            C.preprocessor.intrasample_unique(),
            C.preprocessor.pad_and_random_shift(new_x_dim=20),
            C.preprocessor.featurize_X(F="pep_ECFP4", alphabet="aa", reshape=True),
        ]
    )

    data = C.pipeline.run(Heinis_dataset)
    assert data[0].X.shape == (846, 20, 208)
    assert data[0].X.sum() == 192501


def test_toh5_featurization_Heinis(Heinis_config, tmp_path, Heinis_dataset):
    import clibas as C

    C.initialize(Heinis_config)

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.ml_data = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    C.pipeline.enque(
        [
            C.preprocessor.merge(),
            C.preprocessor.intrasample_unique(),
            C.preprocessor.tt_split(test_fraction=0.1),
            C.preprocessor.to_h5(alphabet="aa", reshape=True, chunks=2),
        ]
    )

    C.pipeline.run(Heinis_dataset)

    fname = str(tmp_path / "outputs" / "train_data.hdf5")

    import h5py

    with h5py.File(fname, "r") as f:
        X = f["X"][:]
        assert X.shape == (762, 14, 21)
        assert X.sum() == 10668

        with pytest.raises(Exception):
            y = f["y"][:]
            assert y.shape == (762,)
