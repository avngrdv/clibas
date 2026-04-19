# About

```{toctree}
:hidden:
:maxdepth: 2
:caption: Getting Started

self
tutorial
documentation
```
```{toctree}
:hidden:
:maxdepth: 1
:caption: Examples
:glob:

examples/*
```

Welcome to clibas – a python package designed for the analysis of next-generation sequencing (NGS) data derived from combinatorial genetically encoded libraries. It is primarily intended for workflows involving mRNA, phage, or yeast display, as well as SELEX affinity selections, but can be readily adapted to related experimental setups where analysis of `.fastq` sequencing data is necessary. 

The package provides fast and efficient tools for parsing, filtering, counting, and analyzing `.fastq` and `.fastq.gz` files containing DNA sequencing reads. Analyses can be performed at both the DNA level and the translated peptide/protein level.

A high-level API makes it straightforward to construct flexible and fully customizable NGS analysis pipelines. Going from a compressed `.fastq.gz` file to clustered peptide results can be achieved in just a few lines of code. Pipelines are concise and human-readable but support sophisticated data analysis workflows.

clibas is performant and scalable. A typical pipeline – including loading `.fastq` files, performing in-silico translation, applying quality filters, collecting statistics, and writing results – takes roughly ~20&nbsp;s per million reads on a single CPU (Intel Core i9-10900K). The architecture supports datasets containing billions of reads.

We have used clibas to make sense of mRNA display affinity selection data, saturation mutagenesis experiments, enzyme substrate profiling workflows, and quantitative enzymatic display pipelines. Some of these and similar pipelines, as well as examples of using clibas to parse phage/yeast display and SELEX selection data, are provided here as Jupyter notebook examples.

For step-by-step instructions on building analysis pipelines, see the Tutorial section. Additional details on the design and capabilities of the software will be described in the accompanying paper (_tbd_).

---

(Installation)=
## Installation

The library can be installed from PyPI:

```bash
pip install clibas[ml]
```

This will also install `scikit-learn`, `rdkit`, `umap-learn`, `hdbscan`, `plotly`, and `h5py` packages. These libraries are used to run UMAP–HDBSCAN sequence embedding and clustering analyses. If these capabilities are not required, a lightweight package can be installed like this:

```bash
pip install clibas
```

Source project is available on GitHub: [clibas](https://github.com/avngrdv/clibas){target=_blank}.

(Contact)=
## Contact

To report bugs, seek technical assistance, and general correspondence, please contact Alex Vinogradov at <vngrdv@nus.edu.sg>, or on [GitHub](https://github.com/avngrdv/clibas){target=_blank}. Please feel free also to request new features. 

Any contributions – code, jupyter notebooks, etc – are welcome! 

(Links)=
## Links

Project on PyPI: link tbd \
Accompanying paper: link tbd \
[Project on GitHub](https://github.com/avngrdv/clibas){target=_blank} \
[Vinogradov Lab at NUS](https://vinogradov.science/){target=_blank}

(Citation)=
## Citation

If find this code useful, please cite the accompanying publication: tbd / link tbd