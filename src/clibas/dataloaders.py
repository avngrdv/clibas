"""
FASTQ file loading and streaming utilities.

Provides FastqLoader for reading sequencing data from FASTQ files, supporting
both uncompressed (.fastq) and gzipped (.fastq.gz) formats with optional
streaming for memory-efficient processing of large files.
"""

import gzip
import os

import numpy as np

from clibas.baseclasses import Handler
from clibas.datatypes import Data, SequencingSample


class FastqLoader(Handler):
    """
    FASTQ file loading and streaming utilities.

    Provides FastqLoader for reading sequencing data from FASTQ files, supporting
    both uncompressed (.fastq) and gzipped (.fastq.gz) formats with optional
    streaming for memory-efficient processing of large files. Typically accessed
    through the clibas facade after initialization.

    Example:
        >>> import clibas as C
        >>> C.initialize_from_config('config.yaml')
        >>> #loader tools are now ready to use
        >>> #as C.data_loader

    Note:
        This class is not typically instantiated directly. Use the clibas
        initialization system to access analysis functionality.
    """

    def __init__(self, *args):
        super(FastqLoader, self).__init__(*args)

        self._validate_designs()
        self._validate_constants()
        return

    def _parse_fastq_lines(self, content, sample_name):
        DNA = content[1::4]
        DNA = np.array([x for x in DNA])

        Q = content[3::4]
        Q = np.array([x for x in Q])

        sample = SequencingSample(name=sample_name, dna=DNA, Q=Q, pep=None)
        # it doesn't matter which design (P or D) is used here
        # because their size should be equivalent either way
        if hasattr(self, "P_design"):
            y_dim = len(self.P_design)
        elif hasattr(self, "D_design"):
            y_dim = len(self.D_design)
        else:
            y_dim = 1

        shape = (sample.size, y_dim)
        sample._internal_state = np.ones(shape, dtype=bool)

        return sample

    def _fetch_fastq_file(self, reader):
        """
        Read and parse a FASTQ file from a file reader object.

        Args:
            reader: Open file reader for a FASTQ file.

        Returns:
            SequencingSample: Parsed sequencing data as a SequencingSample object.

        Note:
            Expects Illumina MiSeq single pair read format. Automatically extracts
            sample name from filename.
        """
        basename = os.path.basename(reader.name)
        sample_name = os.path.splitext(basename)[0]

        with reader as f:
            msg = f"Fetching {basename}. . ."
            self.logger.info(msg)
            content = []
            for line in f:
                content.append(line.rstrip("\n").encode("ascii"))

        return self._parse_fastq_lines(content, sample_name)

    def stream_from_fastq_dir(self, data_dir=None):
        """
        Stream sequencing samples from a directory containing FASTQ files.

        Generator that yields one SequencingSample per FASTQ file, enabling
        memory-efficient processing of multiple files without loading all data
        at once.

        Args:
            data_dir (str, optional): Directory containing FASTQ files. If None,
                uses the sequencing_data directory from the config file.

        Yields:
            SequencingSample: One sample per FASTQ file in the directory.

        Raises:
            IOError: If directory is invalid or contains no .fastq files.

        Example:
            >>> streamer = C.data_loader.stream_from_fastq_dir(data_dir='../fastq')
            >>> #process samples one by one
            >>> C.pipeline.stream(streamer=streamer, save_summary=True)
        """
        d = self.dirs.seq_data
        if data_dir is not None:
            if not os.path.isdir(data_dir):
                msg = '<fetch_gz_from_dir>: specified "data_dir" does not point to a valid directory!'
                self.logger.error(msg)
                raise OSError(msg)
            d = data_dir

        fnames = [os.path.join(d, x) for x in os.listdir(d) if x.endswith(".fastq")]

        if not fnames:
            msg = f"No .fastq files were found in {d}! Aborting."
            self.logger.error(msg)
            raise OSError(msg)

        for f in fnames:
            reader = open(f)
            sample = self._fetch_fastq_file(reader)
            yield sample

    def stream_from_gz_dir(self, data_dir=None):
        """
        Stream sequencing samples from a directory of gzipped FASTQ files.

        Generator that yields one SequencingSample per .fastq.gz file, enabling
        memory-efficient processing of multiple compressed files without loading
        all data at once.

        Args:
            data_dir (str, optional): Directory containing gzipped FASTQ files.
                If None, uses the sequencing_data directory from the config file.

        Yields:
            SequencingSample: One sample per .fastq.gz file in the directory.

        Raises:
            IOError: If directory is invalid or contains no .fastq.gz files.

        Example:
            >>> streamer = C.data_loader.stream_from_gz_dir(data_dir='../fastq')
            >>> #process samples one by one
            >>> C.pipeline.stream(streamer=streamer, save_summary=True)
        """
        d = self.dirs.seq_data
        if data_dir is not None:
            if not os.path.isdir(data_dir):
                msg = '<fetch_gz_from_dir>: specified "data_dir" does not point to a valid directory!'
                self.logger.error(msg)
                raise OSError(msg)
            d = data_dir

        fnames = [os.path.join(d, x) for x in os.listdir(d) if x.endswith(".gz")]

        if not fnames:
            msg = f"No .fastq.gz files were found in {d}! Aborting."
            self.logger.error(msg)
            raise OSError(msg)

        for f in fnames:
            reader = gzip.open(f, "rt")
            sample = self._fetch_fastq_file(reader)
            yield sample

    def stream_from_gz_file(self, fname=None, reads_per_chunk=None):
        """
        Stream a gzipped FASTQ file in chunks for memory-efficient processing.

        Reads a large .fastq.gz file in manageable chunks, yielding SequencingSample
        objects containing the specified number of reads. Useful for processing
        files too large to fit in memory.

        Args:
            fname (str): Path to the .fastq.gz file to stream.
            reads_per_chunk (int): Number of reads to process in each chunk.

        Yields:
            SequencingSample: Sample objects containing reads_per_chunk sequences. Sample names are suffixed with chunk numbers (e.g., sample_001, sample_002).

        Raises:
            ValueError: If reads_per_chunk is not an integer.
            IOError: If file cannot be opened or read.

        Example:
            >>> streamer = C.data_loader.stream_from_gz_file(fname='example.fastq.gz', reads_per_chunk=int(5e6))
            >>> #process file chunks one by one
            >>> C.pipeline.stream(streamer=streamer, save_summary=True)
        """
        if not isinstance(reads_per_chunk, int):
            msg = f'Invalid argument passed to <stream_from_gz_file> op. Expected "reads_per_chunk" as dtype=int, received: {type(reads_per_chunk)}'
            self.logger.error(msg)
            raise ValueError(msg)

        try:
            reader = gzip.open(fname, "rt")
            basename = os.path.basename(reader.name)
            msg = f"Streaming from {basename} file. . ."
            self.logger.info(msg)
        except:
            msg = f"<stream_from_gz_file> op could load {fname}. . ."
            self.logger.error(msg)
            raise OSError(msg)

        n_chunk_lines = reads_per_chunk * 4
        chunk = []
        chunk_counter = 1
        for line in reader:
            chunk.append(line.rstrip().encode("ascii"))
            if len(chunk) == n_chunk_lines:
                sample_name = os.path.splitext(basename)[0] + f"_{chunk_counter:03d}"
                msg = f"Fetching {sample_name}. . ."
                self.logger.info(msg)
                yield self._parse_fastq_lines(chunk, sample_name)

                chunk = []
                chunk_counter += 1

        # return leftovers if any
        if chunk:
            sample_name = os.path.splitext(basename)[0] + f"_{chunk_counter:03d}"
            yield self._parse_fastq_lines(chunk, sample_name)

    def fetch_fastq_from_dir(self, data_dir=None):
        """
        Load all FASTQ files from a directory into memory.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Reads all .fastq files in the specified directory and combines them into
        a single Data object. Use for small to medium datasets that fit in memory.

        Args:
            data_dir (str, optional): Directory containing FASTQ files. If None,
                uses the sequencing_data directory from configuration.

        Returns:
            callable: Operation that when called returns a Data object containing all samples from the directory.

        Raises:
            IOError: If directory is invalid or contains no .fastq files.

        Example:
            >>> fetch_op = C.data_loader.fetch_fastq_from_dir('path/to/fastq/files')
            >>> data = fetch_op()
        """
        if data_dir is not None:
            if not os.path.isdir(data_dir):
                msg = '<fetch_gz_from_dir>: specified "data_dir" does not point to a valid directory!'
                self.logger.error(msg)
                raise OSError(msg)

        def fetch_dir_fastq(*args):
            samples = list()
            for sample in self.stream_from_fastq_dir(data_dir=data_dir):
                samples.append(sample)

            return Data(samples=samples)

        return fetch_dir_fastq

    def fetch_gz_from_dir(self, data_dir=None):
        """
        Load all gzipped FASTQ files from a directory into memory.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Reads all .fastq.gz files in the specified directory and combines them into
        a single Data object. Use for small to medium datasets that fit in memory.

        Args:
            data_dir (str, optional): Directory containing gzipped FASTQ files.
                If None, uses the sequencing_data directory from configuration.

        Returns:
            callable: Operation that when called returns a Data object containing all samples from the directory.

        Raises:
            IOError: If directory is invalid or contains no .fastq.gz files.

        Example:
            >>> fetch_op = C.data_loader.fetch_gz_from_dir('path/to/fastq/files')
            >>> data = fetch_op()
        """
        if data_dir is not None:
            if not os.path.isdir(data_dir):
                msg = '<fetch_gz_from_dir>: specified "data_dir" does not point to a valid directory!'
                self.logger.error(msg)
                raise OSError(msg)

        def fetch_dir_gz(*args):
            samples = list()
            for sample in self.stream_from_gz_dir(data_dir=data_dir):
                samples.append(sample)

            return Data(samples=samples)

        return fetch_dir_gz
