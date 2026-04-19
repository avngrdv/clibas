"""
Base classes and utilities for clibas handlers.

Provides foundational classes for logging, directory management, and data
handling. These classes are not never invoked directly!
"""

import datetime
import logging
import os

import numpy as np


class Logger:
    """
    Configured logging wrapper for clibas operations.

    Provides customizable logging with optional console output and file logging.
    Typically configured through the clibas configuration system rather than
    instantiated directly.

    Args:
        config: Configuration object with logger settings (name, verbose,
            log_to_file, log_fname, level).
    """

    def __init__(self, config=None):
        self.conf = config

        self.__fallback()
        self.__configure_logger()
        return

    def __repr__(self):
        return f"<Logger {self.name}; verbose: {self.verbose}; log_to_file: {self.log_to_file}; level: {self.level}>"

    def __fallback(self):
        """
        If no config is passed fallback to some innocuous defaults, which is
        basically an info console logger (no file logging).
        """

        attribs = {
            "name": "unnamed",
            "verbose": False,
            "log_to_file": False,
            "log_fname": None,
            "level": "INFO",
        }

        for attr in attribs:
            if not hasattr(self.conf, attr):
                setattr(self, attr, attribs[attr])
            else:
                setattr(self, attr, getattr(self.conf, attr))

        return

    def __configure_logger(self):
        levels = {
            "DEBUG": 10,
            "INFO": 20,
            "WARNING": 30,
            "ERROR": 40,
        }

        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(levels[self.level])
        formatter = logging.Formatter("[%(levelname)s]: %(message)s")

        # clear any preexisting handlers to avoid stream duplication
        self.logger.handlers.clear()

        if self.verbose:
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            console_handler.setLevel(levels[self.level])

            self.logger.addHandler(console_handler)

        if self.log_to_file:
            if self.log_fname is None:
                cwd = os.getcwd()
                self.log_fname = os.path.join(cwd, "unnamed_log_file.txt")

            filehandler = logging.FileHandler(self.log_fname)
            filehandler.setFormatter(formatter)
            filehandler.setLevel(levels[self.level])

            self.logger.addHandler(filehandler)

        return


class DirectoryTracker:
    """
    Directory management for clibas file operations.

    Tracks and creates directories for data input/output, logs, and analysis results.
    Automatically creates directories if they don't exist. Typically configured
    through the clibas configuration system.

    Args:
        config: Configuration object specifying directory paths (seq_data, logs,
            parser_out, analysis_out).
    """

    def __init__(self, config=None):
        self._conf = config
        self.__fallback()
        # self.__setup_dirs()

    def __repr__(self):
        return "<DirectoryTracker object>"

    def __fallback(self):
        # any directories left unspecified by config will be set to cwd
        cwd = os.getcwd()

        # all concerned directories
        attribs = ["seq_data", "logs", "parser_out", "analysis_out"]

        for attr in attribs:
            if not hasattr(self._conf, attr):
                setattr(self, attr, cwd)
            else:
                setattr(self, attr, getattr(self._conf, attr))
        return
    
    # this doesn't actually seem necessary -> rely on a lazy folder creation instead
    def __setup_dirs(self):
        for d in [x for x in dir(self) if not x.startswith("_")]:
            if not os.path.isdir(getattr(self, d)):
                os.makedirs(getattr(self, d))
        return


class Handler:
    """
    Base class for clibas data handlers.

    Provides common functionality for data processing operations including logging,
    directory tracking, array manipulation, and validation. Should not be instantiated
    directly - use subclasses like FastqParser instead.

    Note:
        This is an internal base class. Users should interact with subclasses
        like FastqParser, DataPreprocessor, etc.
    """

    def __init__(self, *args):
        self.__dict__.update(*args)
        self.__logger_fallback()
        self.__tracker_fallback()
        return

    def __logger_fallback(self):
        # if no logger is passed to a data handler, a default Logger
        # object will be invoked. The default logger is silent.
        if not hasattr(self, "logger"):
            self.logger = Logger().logger

        if self.logger is None:
            self.logger = Logger().logger
        return

    def __tracker_fallback(self):
        # If no DirectorTracker object was passed to a handler, a default
        # tracker will be invoked (everything in the cwd)
        if not hasattr(self, "dirs"):
            self.dirs = DirectoryTracker()
        return

    def _validate_designs(self):
        if hasattr(self, "P_design") and hasattr(self, "D_design"):
            if not len(self.P_design) == len(self.D_design):
                msg = "Peptide and DNA library designs must contain the same number of templates; cannot inialize FastqProcessor. . ."
                self.logger.error(msg)
                raise ValueError(msg)
        return

    def _validate_constants(self):
        if not hasattr(self, "constants"):
            msg = "FastqParser requires constants for setup. . ."
            self.logger.error(msg)
            raise ValueError(msg)
        return

    def _where_check(self, where):
        if where not in ("pep", "dna"):
            msg = f"The parser did not understand which dataset it should operate on. Passed value: {where}; allowed values: pep/dna."
            self.logger.error(msg)
            raise ValueError(msg)
        return

    def _infer_design(self, where):
        if where == "pep":
            if hasattr(self, "P_design"):
                return self.P_design
            else:
                msg = "Peptide library design not set: cannot analyze peptide datasets without specifying a library design."
                self.logger.error(msg)
                raise ValueError(msg)                

        if where == "dna":
            if hasattr(self, "D_design"):
                return self.D_design
            else:
                msg = "DNA library design not set: cannot analyze DNA datasets without specifying a library design."
                self.logger.error(msg)
                raise ValueError(msg)
 
    
    def _infer_alphabet(self, where=None, alphabet=None):
        
        # if alphabet is specified: override any inferences
        if alphabet is not None:
            if isinstance(alphabet, np.ndarray):
                if alphabet.ndim != 1:
                    msg = f"Token alphabet must be supplied as a 1-dimensional array; recived n_dim={alphabet.ndim}. . ."
                    self.logger.error(msg)
                    raise ValueError(msg)

            if isinstance(alphabet, (np.ndarray, list, tuple)):
                return tuple(i if isinstance(i, bytes) else str(i).encode()
                             for i in alphabet
                            )
            else:
                if alphabet == "aa":
                    return self.constants.aas

                elif alphabet == "base":
                    return self.constants.bases

                else:
                    msg = 'Parameter "alphabet" was expected as type=(list, tuple, np.ndarray) or a keyword ("aa", "base"); received: {type(alphabet)}'
                    self.logger.error(msg)
                    raise ValueError(msg)
                
        # if no alphabet is supplied, infer
        if where == "pep":
            return self.constants.aas

        elif where == "dna":
            return self.constants.bases

        msg = "Token alphabet was not supplied or was not understood. . ."
        self.logger.error(msg)
        raise ValueError(msg)        
    
    
    def _loc_check(self, loc, design):
        if not isinstance(loc, (list, tuple)):
            msg = f"The Parser expected to receive a list of region indexes to parse; received: {type(loc)}"
            self.logger.error(msg)
            raise ValueError(msg)

        if max(loc) > design.loc.max():
            msg = f"{design.lib_type} library design does not contain enough regions. Library design contains {design.loc.max() + 1} regions; specified: up to {max(loc) + 1}"
            self.logger.error(msg)
            raise AssertionError(msg)
        return

    def _prepare_destinations(self, root=None, data=None):
        if data is not None:
            for sample in data:
                destination = os.path.join(root, sample.name)
                if not os.path.isdir(destination):
                    os.makedirs(destination)
        else:
            if not os.path.isdir(root):
                os.makedirs(root)            
        return

    def _L_summary(self, arr):
        # infer what the pad token is
        pad = arr.dtype.type()

        # fetch the indexes where dna/pep length == designed
        return np.sum(arr != pad, axis=1)

    def _find_max_len(self, design, loc):
        """
        When trying to get a column-wise slice of the array, the slices for
        different designs can have a different shape (for example, different
        variable region size). This op will find the largest possible
        column-wise view.  Output m is used to as a shape parameter during
        array creation.
        """
        m = 0
        for template in design:
            if len(template(loc)) > m:
                m = len(template(loc))
        return m

    # helper function to quickly create 1D-views of sequence arrays
    def _arr_as_1d(self, arr):
        if arr.ndim == 1:
            return arr

        if arr.ndim == 2:
            arr = np.ascontiguousarray(arr)
            kind = arr.dtype.kind

            if arr.dtype.kind in ("S", "U"):
                s_len = arr.shape[1] * arr.dtype.itemsize // (1 if kind == "S" else 4)
                return arr.view(f"{kind}{s_len}").reshape(arr.shape[0])

            # there is just one special case: viewing uint8 arrays as SX
            if arr.dtype == np.uint8:
                return arr.view(f"S{arr.shape[1]}").reshape(arr.shape[0])

        msg = f"{self} could not create a 1D array view. . ."
        self.logger.error(msg)
        raise ValueError(msg)

    # helper function to quickly create 2D-views of sequence arrays
    def _arr_as_2d(self, arr):
        if arr.ndim == 2:
            return arr

        if arr.ndim == 1:
            arr = np.ascontiguousarray(arr)
            kind = arr.dtype.kind

            # currently, we can recast only SX, UX and object-type arrays
            if arr.dtype.kind == "O":
                arr = arr.astype(f"U{max(len(s) for s in arr)}")

            if arr.dtype.kind in ("S", "U"):
                return arr.view(f"{kind}1").reshape((arr.shape[0], -1))

        msg = f"{self} could not create a 2D array view. . ."
        self.logger.error(msg)
        raise ValueError(msg)

    def _cast(self, arr, dim):
        if dim == "1d":
            return self._arr_as_1d(arr)

        elif dim == "2d":
            return self._arr_as_2d(arr)

        else:
            msg = f"{self} did not understand array dimensions to cast to"
            self.logger.error(msg)
            raise ValueError(msg)

    def _empty_array_check(self, arr, func):
        if arr.ndim > 0:
            if arr.size > 0:
                return

        msg = f"{func} received an empty dataset to operate on! Aborting. . ."
        self.logger.error(msg)
        raise ValueError(msg)

    def _joint_alphabet_X_check(self, arr, alphabet):
        arr = self._cast(arr, "2d")

        try:
            alphabet = np.asarray(alphabet).astype(arr.dtype)
        except:
            msg = "<DataPreprocessor> could not cast the alphabet and the X-arrays in data to the same datatype. . ."
            self.logger.error(msg)
            raise ValueError(msg)

        exp = np.hstack((alphabet, arr.dtype.type()))  # we have to allow pads in X
        if not np.isin(arr, exp).all():
            msg = "<DataPreprocessor>: Sample's X-array contains elements not encoded in the alphabet!"
            self.logger.error(msg)
            raise ValueError(msg)

        return arr

    @property
    def timestamp(self):
        t = datetime.datetime.now()
        return f"_{t.year}_{t.month}_{t.day}_{t.hour}_{t.minute}_{t.second}"
