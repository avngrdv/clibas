from .dataanalysis import DataAnalysisTools
from .dataloaders import FastqLoader
from .datapreprocessors import DataPreprocessor
from .dispatchers import Dispatcher
from .parsers import FastqParser
from .pipelines import Pipeline

__version__ = "0.4.3"

# initialize later - need to pass the location of the config file
_c_instance = None
    
class _clibas_facade:
    def __init__(self, config_path):
        self._config_path = config_path

        handlers = (
            Pipeline,
            FastqParser,
            DataAnalysisTools,
            DataPreprocessor,
            FastqLoader,
        )
        (
            self._pip,
            self._par,
            self._dta,
            self._pre,
            self._loa,
        ) = Dispatcher.dispatch(config_path, handlers)

        msg = "<clibas> succesfully loaded config and is ready. . ."
        self._pip.logger.info(msg)
        return

    @property
    def pipeline(self):
        return self._pip

    @property
    def fastq_parser(self):
        return self._par

    @property
    def preprocessor(self):
        return self._pre

    @property
    def analysis_tools(self):
        return self._dta

    @property
    def data_loader(self):
        return self._loa


def initialize(config_path=None):
    global _c_instance
    _c_instance = _clibas_facade(config_path)


_PUBLIC_ATTRS = [
    "initialize",
    "pipeline",
    "fastq_parser",
    "preprocessor",
    "analysis_tools",
    "data_loader",
]


def __getattr__(name):
    if _c_instance is None:
        raise RuntimeError(
            "clibas not initialized. Call clibas.initialize() first."
        )

    if name not in _PUBLIC_ATTRS:
        raise AttributeError(f"Module 'clibas' has no attribute '{name}'")

    return getattr(_c_instance, name)


def __dir__():
    return sorted(_PUBLIC_ATTRS)


__all__ = ["initialize"]
