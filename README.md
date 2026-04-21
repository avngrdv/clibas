<div align="center">
  <picture>
    <source srcset="/docs/source/_static/clibas_logo_dark_v1.png" media="(prefers-color-scheme: dark)">
    <img src="/docs/source/_static/clibas_logo_v1.png" alt="clibas logo" width="50%">
  </picture>
</div>

<div align="center">

[![PyPI version](https://img.shields.io/pypi/v/clibas?color=blue)](https://pypi.org/project/clibas/)
[![Python versions](https://img.shields.io/pypi/pyversions/clibas)](https://pypi.org/project/clibas/)
[![Downloads](https://img.shields.io/pypi/dm/clibas)](https://pypi.org/project/clibas/)
[![Docs](https://readthedocs.org/projects/clibas/badge/?version=latest)](https://clibas.readthedocs.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

</div>

# clibas

Welcome to clibas – a Python package for analyzing NGS data from combinatorial genetically encoded libraries, including techniques like mRNA/phage/yeast display and SELEX selections. The library provides fast and scalable tools for parsing, filtering, and analyzing .fastq files at both DNA and translated peptide levels, with a high-level API to build sophisticated analysis pipelines in just a few lines of code.

## Documentation & examples
For full documentation, including a tutorial, API reference, and jupyter notebook examples, please visit our <a href="https://clibas.readthedocs.io/en/latest/index.html#" target="_blank" rel="noopener noreferrer">Read the Docs</a> page. Example .ipynb notebooks are also available in this repo in `docs\source\examples`

## Quick start

It is recommended that clibas is installed in a dedicated virtual environment to avoid potential version conflicts with existing packages. Any virtual environment (e.g., `conda` or `pipenv`) will work for this purpose. 

The library can be installed from PyPI:

```bash
pip install clibas[ml]
```

This will also install `scikit-learn`, `rdkit`, `umap-learn`, `hdbscan`, `plotly`, and `h5py` packages. These libraries are used to run UMAP–HDBSCAN sequence embedding and clustering analyses. If these capabilities are not required, a lightweight package can be installed like this:

```bash
pip install clibas
```

## Contact

To report bugs, seek technical assistance, and general correspondence, please contact Alex Vinogradov at <vngrdv@nus.edu.sg>, or here on github. 

Any contributions – code, feature requests, jupyter notebooks – are welcome! 

## Links

[Project on PyPI](https://pypi.org/project/clibas/) \
[Documentation](https://clibas.readthedocs.io/en/latest/index.html#) \
Accompanying paper: link tbd \
[Vinogradov Lab at NUS](https://vinogradov.science/)

## Citation

If you use clibas, please cite the accompanying publication: *tbd / link tbd*