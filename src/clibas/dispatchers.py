"""
Config file parsing and handler dispatching.

Provides Dispatcher for parsing YAML config files and initializing
clibas handlers with coordinated settings. Internal infrastructure:
not typically used directly by end users.
"""

import datetime
import importlib
import os

import yaml

import clibas.config_defaults
from clibas.baseclasses import DirectoryTracker, Logger

_reserved_names = ("_", "+", "*", "1", "2", "3", "4", "5", "6", "7", "8", "9", "0")


class Dispatcher:
    """
    Config file parsing and handler dispatching.

    Provides Dispatcher for parsing YAML config files and initializing
    clibas handlers with coordinated settings. Internal infrastructure:
    not typically used directly by end users.
    """

    def __init__(self, config_fname):
        self._preparse_yaml_config(config_fname)
        self._parse_config()
        return

    def __repr__(self):
        return f"<Dispatcher object at {hex(id(self))}>"

    def _preparse_yaml_config(self, config_fname):
        # if nothing at all is passed - it's fine
        if not config_fname:
            self.config = {}
            return

        # but a faulty file should raise
        try:
            with open(config_fname) as file:
                raw = yaml.safe_load(file)
        except:
            msg = "<clibas.Dispatcher>: could not load the config file. . ."
            raise OSError(msg)

        config_dict = {}
        for config_type, attribs in raw.items():
            if not isinstance(attribs, dict):
                attribs = {"value": attribs}

            config_dict[config_type] = type(config_type, (object,), attribs)

        self.config = config_dict
        return

    def _parse_dir_config(self):
        if "TrackerConfig" not in self.config.keys():
            self.dirs = DirectoryTracker()

        else:
            # don't have to check individual attribs becase DirectoryTracker
            # can handle it by itself
            self.dirs = DirectoryTracker(self.config["TrackerConfig"])
        return

    def _parse_logger_config(self):
        # logger: just need to set up the log file location, if needed
        if "LoggerConfig" in self.config.keys():
            if hasattr(self.config["LoggerConfig"], "log_to_file"):
                if self.config["LoggerConfig"].log_to_file:
                    if not hasattr(self.config["LoggerConfig"], "log_fname"):
                        t = datetime.datetime.now()
                        timestamp = f"_{t.year}_{t.month}_{t.day}_{t.hour}_{t.minute}_{t.second}"

                        if not os.path.isdir(self.dirs.logs):
                            os.makedirs(self.dirs.logs)

                        self.config["LoggerConfig"].log_fname = os.path.join(
                            self.dirs.logs,
                            self.config["experiment"].value + timestamp + "_logs.txt",
                        )

            self.L = Logger(config=self.config["LoggerConfig"]).logger

        else:
            self.L = Logger().logger
        return

    def _validate_translation_table(self, table):
        inter = set(table.values()) & set(_reserved_names)
        if len(inter) == 1:
            if inter.pop() == "*":
                return

        msg = "<Dispatcher>: Amino acids in the translation table can not be encoded by reserved names. . . "
        self.L.error(msg)
        raise ValueError(msg)
        return

    def _parse_constants_config(self):
        # setup constants: if no constants are given, fallback to defaults right away
        if "constants" in self.config:
            self.constants = self.config["constants"]
            if not all(
                (
                    hasattr(self, "aas"),
                    hasattr(self, "translation_table"),
                    hasattr(self, "bases"),
                )
            ):
                msg = "<Dispatcher>: config is missing some basic definitions (bases, amino acids, translation tables, etc); will use defaults for the missing values. . . "
                self.L.debug(msg)

        else:
            msg = "<Dispatcher>: config is missing some basic definitions (bases, amino acids, translation tables, etc); will use defaults for the missing values. . . "
            self.L.debug(msg)

            import copy
            import types

            self.constants = types.SimpleNamespace(
                **{
                    k: copy.deepcopy(v)
                    for k, v in vars(clibas.config_defaults).items()
                    if not k.startswith("_")
                }
            )

        b_def_aas = tuple(s.encode("ascii") for s in clibas.config_defaults.aas)
        b_def_base = tuple(s.encode("ascii") for s in clibas.config_defaults.bases)

        # deal with amino acids
        if hasattr(self.constants, "aas"):
            if set(self.constants.aas) & set(_reserved_names):
                msg = "<Dispatcher>: Amino acid monomers can not be encoded by reserved names. . . "
                self.L.error(msg)
                raise ValueError(msg)

            self.constants.aas = tuple(
                sorted(s.encode("ascii") for s in self.constants.aas)
            )

        # amino acids can also be inferred from the translation table
        elif hasattr(self.constants, "translation_table"):
            self._validate_translation_table(self.constants.translation_table)
            self.constants.aas = tuple(
                sorted(
                    set(
                        x.encode("ascii")
                        for x in self.constants.translation_table.values()
                        if x not in _reserved_names
                    )
                )
            )

        # if no amino acid info is passed, fallback to defaults
        else:
            self.constants.aas = b_def_aas

        # set up translation table
        if hasattr(self.constants, "translation_table"):
            table = self.constants.translation_table
        else:
            table = clibas.config_defaults.translation_table

        self._validate_translation_table(table)
        self.constants.translation_table = {
            k.encode("ascii"): v.encode("ascii") for k, v in table.items()
        }

        # set up custom initiator amino acid: if any; there is no fallback for this.
        if hasattr(self.constants, "custom_ini_aa"):
            if self.constants.custom_ini_aa in _reserved_names:
                msg = "<Dispatcher>: custom initiator amino acid can not be encoded by a reserved name. . . "
                self.L.error(msg)
                raise ValueError(msg)

            self.constants.custom_ini_aa = self.constants.custom_ini_aa.encode("ascii")

            if self.constants.custom_ini_aa not in self.constants.aas:
                self.constants.aas = tuple(
                    sorted(self.constants.aas + (self.constants.custom_ini_aa,))
                )

        # deal with aa_SMILES: make sure that the aa_SMILES and aas alphabets match.
        if hasattr(self.constants, "aa_SMILES"):
            if not isinstance(self.constants.aa_SMILES, dict):
                msg = f"<Dispatcher>: config must supply aa_SMILES as a dictionary, found {type(self.constants.aa_SMILES)} instead. . . "
                self.L.error(msg)
                raise ValueError(msg)

            sorted_keys = tuple(
                sorted(list(s.encode("ascii") for s in self.constants.aa_SMILES.keys()))
            )
            if sorted_keys != self.constants.aas:
                msg = "<Dispatcher>: Amino acid alphabet does not match the aa_SMILES dictionary. . . "
                self.L.error(msg)
                raise ValueError(msg)

            self.constants.aa_SMILES = tuple(
                self.constants.aa_SMILES[k.decode("ascii")] for k in sorted_keys
            )

        # if aa_SMILES is not explicitly specified, but aas are all standard, we
        # can still rescue this:
        elif self.constants.aas == b_def_aas:
            self.constants.aa_SMILES = tuple(
                clibas.config_defaults.aa_SMILES[k.decode("ascii")]
                for k in self.constants.aas
            )

        # parse the barcode_table if present
        if hasattr(self.constants, "barcode_table"):
            self.constants.barcode_table = {
                tuple(s.encode("ascii") for s in k): {
                    i.encode("ascii"): j.encode("ascii") for i, j in v.items()
                }
                for k, v in self.constants.barcode_table.items()
            }

        # now deal with nucleic acids
        if hasattr(self.constants, "bases"):
            if set(self.constants.bases) & set(_reserved_names):
                msg = "<Dispatcher>: Nucleotide monomers can not be encoded by reserved names. . . "
                self.L.error(msg)
                raise ValueError(msg)

            self.constants.bases = tuple(
                sorted(s.encode("ascii") for s in self.constants.bases)
            )

        else:
            self.constants.bases = b_def_base

        # deal with base_SMILES: make sure that the base_SMILES and base alphabets match.
        if hasattr(self.constants, "base_SMILES"):
            if not isinstance(self.constants.base_SMILES, dict):
                msg = f"<Dispatcher>: config must supply base_SMILES as a dictionary, found {type(self.constants.aa_SMILES)} instead. . . "
                self.L.error(msg)
                raise ValueError(msg)

            sorted_keys = tuple(
                sorted(
                    list(s.encode("ascii") for s in self.constants.base_SMILES.keys())
                )
            )
            if sorted_keys != self.constants.bases:
                msg = "<Dispatcher>: Nucleotide alphabet does not match the base_SMILES dictionary. . . "
                self.L.error(msg)
                raise ValueError(msg)

            self.constants.base_SMILES = tuple(
                self.constants.base_SMILES[k.decode("ascii")] for k in sorted_keys
            )

        # if base_SMILES is not explicitly specified, but bases are all standard,
        # we can still rescue this:
        elif self.constants.bases == b_def_base:
            self.constants.base_SMILES = tuple(
                clibas.config_defaults.base_SMILES[k.decode("ascii")]
                for k in self.constants.bases
            )

        # setup the complement table: only really needed if revcom is used
        if hasattr(self.constants, "complement_table"):
            t = self.constants.complement_table
            self.constants.complement_table = bytes.maketrans(
                t[0].encode("ascii"), t[1].encode("ascii")
            )
        else:
            self.constants.complement_table = clibas.config_defaults.complement_table
        return

    def _parse_lib_design_config(self):
        from clibas.lib_design import LibraryDesign

        if "LibraryDesigns" not in self.config:
            msg = "<Dispatcher>: config has not specified any library designs. . ."
            self.L.warning(msg)
            return

        # deal with peptide libraries first
        if hasattr(self.config["LibraryDesigns"], "pep_templates") and hasattr(
            self.config["LibraryDesigns"], "pep_monomers"
        ):
            if hasattr(self.constants, "aas"):
                pep_templates = list(
                    s.encode("ascii")
                    for s in self.config["LibraryDesigns"].pep_templates
                )
                pep_monomers = {
                    k: [s.encode("ascii") for s in v]
                    for k, v in self.config["LibraryDesigns"].pep_monomers.items()
                }

                self.P_design = LibraryDesign(
                    templates=pep_templates,
                    monomers=pep_monomers,
                    lib_type="pep",
                    val_monomer=self.constants.aas,
                )

            # note that the condition below is very difficult to violate at the
            # moment, but let's just keep it for now.
            else:
                msg = "Dispatcher: cannot setup a peptide library design without the amino acid alphabet specifications. . ."
                self.L.error(msg)
                raise ValueError(msg)

        # if only of the two necessary atrributes is present, raise
        # note that if both are absent, it's OK
        elif hasattr(self.config["LibraryDesigns"], "pep_templates") != hasattr(
            self.config["LibraryDesigns"], "pep_monomers"
        ):
            msg = "Dispatcher: peptide library design specifications parsed from config are incomplete!"
            self.L.error(msg)
            raise ValueError(msg)

        # same, but now for DNA libs
        if hasattr(self.config["LibraryDesigns"], "dna_templates") and hasattr(
            self.config["LibraryDesigns"], "dna_monomers"
        ):
            if hasattr(self.constants, "bases"):
                dna_templates = list(
                    s.encode("ascii")
                    for s in self.config["LibraryDesigns"].dna_templates
                )
                dna_monomers = {
                    k: [s.encode("ascii") for s in v]
                    for k, v in self.config["LibraryDesigns"].dna_monomers.items()
                }

                self.D_design = LibraryDesign(
                    templates=dna_templates,
                    monomers=dna_monomers,
                    lib_type="dna",
                    val_monomer=self.constants.bases,
                )
            else:
                msg = "Dispatcher: cannot setup a DNA library design without the nucleotide alphabet specification. . ."
                self.L.error(msg)
                raise ValueError(msg)

        # if only of the two necesasry atrributes is present, raise
        elif hasattr(self.config["LibraryDesigns"], "dna_templates") != hasattr(
            self.config["LibraryDesigns"], "dna_monomers"
        ):
            msg = "Dispatcher: config is missing necessary parameters to setup DNA library designs. . ."
            self.L.error(msg)
            raise ValueError(msg)

        return

    def _parse_config(self):
        # check that all necessary config params are in place, one by one
        # if not, fall back to some innocuous defaults where possible
        if "experiment" in self.config:
            self.experiment = self.config["experiment"].value
        else:
            self.experiment = "untitled_experiment"

        self._parse_dir_config()
        self._parse_logger_config()
        self._parse_constants_config()
        self._parse_lib_design_config()

        self.common = {
            key: getattr(self, attr)
            for key, attr in {
                "logger": "L",
                "dirs": "dirs",
                "exp_name": "experiment",
                "constants": "constants",
                "P_design": "P_design",
                "D_design": "D_design",
            }.items()
            if hasattr(self, attr)
        }

        return

    def _config_to_dict(self, conf):
        return dict(
            (name, getattr(conf, name))
            for name in dir(conf)
            if not name.startswith("__")
        )

    def _dispatch_handler(self, handler):
        lookup_modules = [
            "parsers",
            "dataanalysis",
            "datapreprocessors",
            "pipelines",
            "dataloaders",
        ]

        config_mapping = {
            "Pipeline": "PipelineConfig",
            "FastqParser": "FastqParserConfig",
            "DataAnalysisTools": "DataAnalysisConfig",
            "DataPreprocessor": "PreproConfig",
            "FastqLoader": "LoaderConfig",
        }

        for module in lookup_modules:
            try:
                m = importlib.import_module(f"clibas.{module}")
                obj = getattr(m, handler.__name__)
                break
            except:
                obj = None

        if not obj:
            msg = f"<Dispatcher> could not set up {handler}"
            self.L.error(msg)
            raise ImportError(msg)

        params = dict()
        if config_mapping[handler.__name__] in self.config:
            params = self._config_to_dict(self.config[config_mapping[handler.__name__]])

        params.update(self.common)
        return obj(params)

    @classmethod
    def dispatch(cls, config_fname, handlers):
        """
        Load configuration and instantiate handlers.

        The primary/only public method. Parses a config file and initializes
        specified handlers with consistent settings. All handlers share logger,
        directory tracker, and other common config parameters.

        Args:
            config_fname (str): Path to YAML configuration file.

            handlers (class or tuple): Handler class(es) to instantiate.
                Can be a single handler (e.g., FastqParser) or tuple of
                handlers (e.g., (Pipeline, FastqParser, DataPreprocessor)).

        Returns:
            tuple or object: Tuple of initialized handler instances if multiple
                handlers provided, otherwise single handler instance.

        Example:
            >>> from clibas.dispatchers import Dispatcher
            >>> from clibas.parsers import FastqParser
            >>> parser = Dispatcher.dispatch('config.yaml', FastqParser)
            >>>
            >>> #or multiple handlers
            >>> handlers = (Pipeline, FastqParser, DataPreprocessor)
            >>> pip, par, pre = Dispatcher.dispatch('config.yaml', handlers)
        """
        self = cls(config_fname)

        # make sure that either individual handlers or tuples thereof work
        handlers = [handlers]
        try:
            handlers = [item for sublist in handlers for item in sublist]
        except:
            pass

        h = list()
        for handler in handlers:
            dispatched_handler = self._dispatch_handler(handler)
            h.append(dispatched_handler)

        return tuple(h)
