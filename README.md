<div align="center">
  <picture>
    <source srcset="/docs/source/_static/clibas_logo_dark_v1.png" media="(prefers-color-scheme: dark)">
    <img src="/docs/source/_static/clibas_logo_v1.png" alt="clibas logo" width="50%">
  </picture>
</div>

# clibas

### CAREFUL: WORK IN PROGRESS

Welcome to clibas – a Python package for analyzing NGS data from combinatorial genetically encoded libraries, including techniques like mRNA/phage/yeast display and SELEX selections. The library provides fast and scalable tools for parsing, filtering, and analyzing .fastq files at both DNA and translated peptide levels, with a high-level API to build sophisticated analysis pipelines in just a few lines of code.

## Quick start

The library can be installed from PyPI:

```bash
pip install clibas[ml]
```

This will also install `scikit-learn`, `rdkit`, `umap-learn`, `hdbscan`, `plotly`, and `h5py` packages. These libraries are used to run UMAP–HDBSCAN sequence embedding and clustering analyses. If these capabilities are not required, a lightweight package can be installed like this:

```bash
pip install clibas
```

## Documentation & examples
For full documentation, including a tutorial, API reference, and jupyter notebook examples, please visit our *documentation site*. Example .ipynb notebooks are also available in this repo in `docs\source\examples`

## Contact

To report bugs, seek technical assistance, and general correspondence, please contact Alex Vinogradov at <vngrdv@nus.edu.sg>, or here on github. 

Any contributions – code, feature requests, jupyter notebooks – are welcome! 

## Links

Project on PyPI: link tbd \
Documentation: link tbd \
Accompanying paper: link tbd \
<a href="https://vinogradov.science/" target="_blank" rel="noopener noreferrer">Vinogradov Lab at NUS</a>


## Citation

If find this code useful, please cite the accompanying publication: _tbd_ / _link tbd_