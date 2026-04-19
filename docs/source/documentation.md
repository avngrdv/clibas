# Documentation

## Initialization

All functionality is accessed through the `clibas` facade after initialization:

```python
import clibas as C
C.initialize('config.yaml')
```

`initialize()` reads the config file and sets up all tools. A no-config initialization is also possible:

```python
import clibas as C
C.initialize()
```

After this call, the following attributes are available:

| Attribute | Description |
|---|---|
| `C.fastq_parser` | Parsing, filtering, and manipulations with .fastq data |
| `C.analysis_tools` | Data analysis and statistics |
| `C.pipeline` | Pipeline construction and execution |
| `C.preprocessor` | Data preprocessing tools for ML applications |
| `C.data_loader` | Loading .fastq and .fastq.gz data from disk |

## Basic usage
```python
import clibas as C
C.initialize('config.yaml')

C.pipeline.enque(
[
    C.fastq_parser.translate(stop_readthrough=False),
    C.analysis_tools.length_analysis(where='pep'),
    C.fastq_parser.len_filter(where='pep'),
	...
]
)

loader = C.data_loader.fetch_gz_from_dir(data_dir='./sequencing_data/')
C.pipeline.load_and_run(loader=loader, save_summary=True)
```

Refer to each module's API reference below for full method documentation.

## Module API reference

- {doc}`C.fastq_parser <parsers>`
- {doc}`C.analysis_tools <dataanalysis>`
- {doc}`C.pipeline <pipelines>`
- {doc}`C.preprocessor <datapreprocessors>`
- {doc}`C.data_loader <dataloaders>`

```{toctree}
:hidden:

parsers
dataanalysis
pipelines
datapreprocessors
dataloaders
```

