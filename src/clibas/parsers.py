"""
Parser objects for data processing.

At present, the module implements FastqParser for processing NGS data with
operations for in silico translation, filtering, quality control, and
basic statistics.

FastqParser's key capabilities:

    - DNA to peptide translation with ORF detection
    - Multiple filtration operations
    - Demultiplexing and barcode translation
    - Dataset summary operations
"""

import inspect
import os
import re

import numpy as np
import pandas as pd

from clibas.baseclasses import Handler
from clibas.datatypes import Data


class FastqParser(Handler):
    """
    Processor for FASTQ sequencing data with filtering and translation capabilities.

    FastqParser provides methods for processing NGS data, including DNA-to-peptide
    translation, quality filtering, and sequence validation against library designs.
    Most operations work on Data objects and return transformed Data instances.

    Most methods are pipeline operations: they return a callable for use in
    processing pipelines.

    The parser is typically accessed through the clibas facade after initialization:

    Example:
        >>> import clibas as C
        >>> C.initialize_from_config('config.yaml')
        >>> #parser is now ready to use
        >>> #as C.fastq_parser

    Note:
        This class is not typically invoked directly. Use the clibas
        initialization system to access parser functionality.
    """

    def __init__(self, *args):
        super(FastqParser, self).__init__(*args)

        self._validate_designs()
        self._validate_constants()

        if hasattr(self, "orf_locator"):
            self.orf_locator = self.orf_locator.encode("ascii")
        
        if hasattr(self, "sample_barcodes"):
            self._validate_barcodes()
        return

    def __repr__(self):
        return "<FastqParser object>"

    def _validate_barcodes(self):
        if not isinstance(self.sample_barcodes, dict):
            msg = f"<FastqParser> received faulty sample barcodes; expected type=dict, received: {type(self.sample_barcodes)}"
            self.logger.error(msg)
            raise ValueError(msg)

        validated_dict = {}
        barcodes = list(self.sample_barcodes.values())
        names = list(self.sample_barcodes.keys())
        for i in range(len(barcodes)):
            bc = barcodes[i]
            
            if not isinstance(bc, str):
                msg = f"<FastqParser> received faulty sample barcodes; expected individual barcodes of str type, received: {type(bc)}"
                self.logger.error(msg)
                raise ValueError(msg)    
                
            bc = tuple(s.encode('ascii') for s in list(bc))
            if not all(x in self.constants.bases for x in bc):
                msg = "<FastqParser> received barcodes containing unrecognized symbols. . ."
                self.logger.error(msg)
                raise ValueError(msg)      
                
            validated_dict[names[i]] = bc
            
        self.sample_barcodes = validated_dict
        return 

    def _dna_to_pep_v2(
        self, seq, first_aa=None, force_at_frame=None, stop_readthrough=False
    ):
        def find_orf(seq):
            loc = re.search(self.orf_locator, seq)
            if loc is not None:
                return seq[loc.end() :]
            else:
                return None

        pep = b""
        # figure out what to use as the orf
        if force_at_frame is None:
            orf = find_orf(seq)
            if orf:
                pep += first_aa
        else:
            orf = seq[force_at_frame:]

        # throughout, '+' is a reserved symbol to denote failed sequences
        # mostly for ambiguous peptides stemming from the occasional "N" base calls
        if orf is not None:
            for i in range(0, len(orf), 3):
                try:
                    next_aa = self.constants.translation_table[orf[i : i + 3]]

                    if not stop_readthrough and next_aa == b"*":
                        break

                    pep += next_aa

                except KeyError:
                    if len(orf[i : i + 3]) != 3:
                        pep += b"_"
                    else:
                        pep += b"+"
        return pep

    def trim_reads(self, left="", right="", tol=None):
        """
        Trim DNA reads based on constant flanking sequences.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Removes adapter overhangs and non-target sequences by anchoring to specified
        5' and 3' constant regions. Useful for processing merged paired-end reads from
        tools like FLASH that may leave adapter sequences at read termini.

        The operation locks onto the specified left (5') and right (3') sequences,
        keeping only the region between them. Matches must be within the specified
        tolerance. If no match is found on either side, that side remains untrimmed.

        Args:
            left (str): 5' anchor sequence. Everything upstream of this sequence
                is trimmed. Use empty string '' to skip 5' trimming.

            right (str): 3' anchor sequence. Everything downstream of this sequence
                is trimmed. Use empty string '' to skip 3' trimming.

            tol (int): Maximum number of mismatches allowed when matching anchor
                sequences (Hamming distance tolerance).

        Returns:
            callable: Operation that accepts a Data object, applies trimming transformation, and returns the modified Data object.

        Example:
            >>> trimmer = C.fastq_parser.trim_reads(left='ATCG', right='GCTA', tol=1)
            >>> data = trimmer(data)
        """

        if not isinstance(tol, int):
            msg = f"<trim_reads> expected to receive parameter tol as as int; received: {type(tol)}"
            self.logger.error(msg)
            raise ValueError(msg)

        if not isinstance(left, str):
            msg = f'<trim_reads> expected to receive parameters "left" as a string; received: {type(left)}'
            self.logger.error(msg)
            raise ValueError(msg)

        if not isinstance(right, str):
            msg = f'<trim_reads> expected to receive parameters "right" as a string; received: {type(right)}'
            self.logger.error(msg)
            raise ValueError(msg)

        # we need to turn the arguments to byte-strings
        left = np.array(list(left), dtype="S1").view(np.uint8)
        right = np.array(list(right), dtype="S1").view(np.uint8)

        from clibas.cython_hacks import idx_finder

        def trim_exp(data):
            for sample in data:
                self._empty_array_check(sample.dna, inspect.stack()[0][3])
                for i, read in enumerate(self._cast(sample.dna, "2d")):
                    idx1, idx2 = idx_finder(
                        np.frombuffer(read, dtype=np.uint8), left, right, tol
                    )

                    sample.dna[i] = self._cast(sample.dna, "1d")[i][idx1:idx2]
                    sample.Q[i] = self._cast(sample.Q, "1d")[i][idx1:idx2]

            return data

        return trim_exp

    def translate(self, force_at_frame=None, stop_readthrough=False):
        """
        Translate DNA sequences to peptides/proteins in silico.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Performs translation of DNA sequencing data to peptides/proteins. Supports
        ORF detection via a configurable sequence pattern or forced translation at
        a specified frame. Designed for one-ORF-per-read NGS data rather than long
        reads with multiple ORFs.
        
        An underscore '_' is appended to the C-terminus if the final codon is 
        incomplete and no stop codon is present.

        Args:
            force_at_frame (int, optional): Translation frame (0, 1, or 2). If None,
                performs automatic ORF search using the pattern specified in the
                config (orf_locator), in which case translation begins downstream
                of the matched pattern; the search is conducted in the 5'-to-3' direction,
                and the first match is used to initiate translation. If no orf_locator is 
                specified, the operation defaults to "ATG" as orf_locator.
                
                If force_at_frame is specified, it forces translation to start at 
                the given frame regardless of pattern matching.
               
                For example:
                    >>>                DNA: TACGACTCACTATAGGGTTAACTTTAAGAAGGA
                    >>>   force_at_frame=0  ---------->
                    >>>     force_at_frame=1  ---------->
                    >>>      force_at_frame=2  ---------->

            stop_readthrough (bool): If True, translation continues past stop codons
                until the 3' end of the read. If False, translation terminates at the
                first stop codon. Default is False.

        Returns:
            callable: Operation that accepts a Data object, translates DNA to peptides, and returns the modified Data object with peptide sequences.

        Example:
            >>> translator = C.fastq_parser.translate(force_at_frame=0)
            >>> data = translator(data)
        """
        if type(stop_readthrough) != bool:
            msg = f'<translate> op expected to receive param "stop_readthrough" as type=bool; received: {type(stop_readthrough)}'
            self.logger.error(msg)
            raise ValueError(msg)
            
        if force_at_frame is not None:
            if not isinstance(force_at_frame, int):
                msg = f'<translate> op expected to receive param "force_at_frame" as dtype=int; received: {type(force_at_frame)}'
                self.logger.error(msg)
                raise ValueError(msg)
        else:
            if not hasattr(self, "orf_locator"):
                msg = "The ORF locator sequence is not set for the <translate> op; will default to 'ATG'. . ."
                self.logger.warning(msg)
                self.orf_locator = b"ATG"

        # this overhead is introduced to handle custom translation initiation
        if hasattr(self.constants, "custom_ini_aa"):
            first_aa = self.constants.custom_ini_aa
        else:
            if hasattr(self, "orf_locator"):
                first_aa = self.constants.translation_table[self.orf_locator[-3:]]
            else:
                # because we default to 'ATG' above, this never actually happens, but keep for now
                msg = "<translate> op did not receive information about the amino acid used for translation initiation; will default to 'M'. . ."
                self.logger.debug(msg)
                first_aa = b"M"

        def translate_dna(data):
            for sample in data:
                self._empty_array_check(sample.dna, inspect.stack()[0][3])

                sample.pep = np.array(
                    [
                        self._dna_to_pep_v2(
                            x,
                            first_aa=first_aa,
                            force_at_frame=force_at_frame,
                            stop_readthrough=stop_readthrough,
                        )
                        for x in self._cast(sample.dna, "1d")
                    ]
                )
            return data

        return translate_dna

    def revcom(self):
        """
        Reverse complement DNA sequences.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Computes the reverse complement of DNA sequences and reverses the corresponding
        quality scores for each sample in the Data object. Useful for processing reads
        sequenced from the opposite strand.

        Returns:
            callable: Operation that accepts a Data object, applies reverse complement transformation to DNA sequences and reverses quality scores, and returns the modified Data object.

        Note:
            Calling on Data containing peptide sequences will raise a warning;
            The peptide data will be left as is (no transform will be applied)

        Example:
            >>> revcom_op = C.fastq_parser.revcom()
            >>> data = revcom_op(data)
        """

        if not hasattr(self.constants, "complement_table"):
            msg = "Complement table is not specified for <revcom> op. Aborting. . ."
            self.logger.error(msg)
            raise ValueError(msg)

        @np.vectorize
        def _rc(seq):
            return seq.translate(self.constants.complement_table)[::-1]

        @np.vectorize
        def _r(seq):
            return seq[::-1]

        def revcom_data(data):
            for sample in data:
                if sample.pep:
                    msg = "Attempting to to revcom a sample holding a pep dataset. Pep dataset will be ignored. . ."
                    self.logger.warning(msg)

                self._empty_array_check(sample.dna, inspect.stack()[0][3])
                self._empty_array_check(sample.Q, inspect.stack()[0][3])

                sample.dna = _rc(self._cast(sample.dna, "1d"))
                sample.Q = _r(self._cast(sample.Q, "1d"))

            return data

        return revcom_data

    def len_filter(self, where=None, len_range=None):
        """
        Filter dataset sequences by length.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Filters sequences based on library design specifications or a custom length range.
        Sequences that don't meet length criteria are discarded from all dataset in
        each sample.

        Args:
            where (str): Dataset to evaluate for filtering. Must be 'dna' or 'pep'.
                Length is checked against this dataset, but filtering is applied across
                all associated data in each sample.

            len_range (tuple or list, optional): Two-element sequence specifying
                (min_length, max_length) range. If None, filtering is performed
                according to library design specifications.

        Returns:
            callable: Operation that accepts a Data object, applies length filtering, and returns the modified Data object.

        Example:
            >>> #filter by design specifications
            >>> length_filt = C.fastq_parser.len_filter(where='pep')
            >>> data = length_filt(data)
            >>> #filter by custom range
            >>> length_filt = C.fastq_parser.len_filter(where='pep', len_range=(10, 50))
            >>> data = length_filt(data)
        """
        self._where_check(where)
        design = self._infer_design(where)

        if len_range is not None:
            if not isinstance(len_range, (list, tuple)):
                msg = f"<len_filter> op expected to receive len_range argument as a list; received: {type(len_range)}"
                self.logger.error(msg)
                raise ValueError(msg)

            if len(len_range) != 2:
                msg = f"<len_filter> op expected to receive len_range as a list with two values; received: len={len(len_range)}"
                self.logger.error(msg)
                raise ValueError(msg)

        def length_filter(data):
            for sample in data:
                arr = self._cast(sample[where], "2d")
                self._empty_array_check(arr, inspect.stack()[0][3])

                # L is a length summary array
                L = self._L_summary(arr)

                # change the sample internal state
                for i, template in enumerate(design):
                    row_mask = sample._internal_state[:, i]

                    if len_range is None:
                        sample._internal_state[row_mask, i] = L[row_mask] == template.L
                    else:
                        sample._internal_state[row_mask, i] = (
                            L[row_mask] > len_range[0]
                        ) & (L[row_mask] < len_range[1])

                # keep every entry that has at least one positive
                # value in the internal state array
                ind = np.any(sample._internal_state, axis=-1)
                sample.ind_filter(ind)

            return data

        return length_filter

    def cr_filter(self, where=None, loc=None, tol=1):
        """
        Filter datasets by constant region integrity.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Filters sequences based on constant region integrity. Sequences with mutations
        in specified constant regions exceeding the tolerance threshold are discarded
        from all datasets in each sample.

        Args:
            where (str): Dataset to evaluate for filtering. Must be 'dna' or 'pep'.
                Constant regions are checked in this dataset, but filtering is applied
                across all associated data in each sample.

            loc (list): List of integers specifying which constant region indices to
                operate on. Must reference constant (non-variable) regions in the
                library design.

            tol (int): Maximum number of mismatches allowed in the constant regions
                (Hamming distance). Sequences exceeding this threshold are discarded.
                Default is 1.

        Returns:
            callable: Operation that accepts a Data object, applies constant region filtering, and returns the modified Data object.

        Note:
            Insertions and deletions in constant regions are not validated by this
            operation, only substitutions.

        Example:
            Consider a library design with the following structure::

                        seq:  ACDEF11133211AWVFRTQ12345YTPPK
                     region:  [-0-][---1--][--2--][-3-][-4-]
                is_variable:  False  True   False True False

            Regions 0, 2, and 4 are constant regions with expected sequences.

            >>> #allow up to 1 mutation in region 2 (AWVFRTQ)
            >>> cr_filt = C.fastq_parser.cr_filter(where='pep', loc=[2], tol=1)
            >>> data = cr_filt(data)
            >>> #sequences with >1 mutation in AWVFRTQ are discarded
        """
        self._where_check(where)
        design = self._infer_design(where)
        self._loc_check(loc, design)

        if not isinstance(tol, int):
            msg = f"<constant_region_filter> expected to receive parameter tol as as int; received: {type(tol)}"
            self.logger.error(msg)
            raise ValueError(msg)

        if np.any(design.is_vr[loc]):
            msg = "<constant_region_filter> expected a list of contant regions to operate on; some of the specified locations point to variable regions."
            self.logger.error(msg)
            raise AssertionError(msg)

        def constant_region_filter(data):
            from clibas.misc import hamming_distance

            for sample in data:
                arr = self._cast(sample[where], "2d")
                self._empty_array_check(arr, inspect.stack()[0][3])

                # iterativelt fill in the indexing array
                for i, template in enumerate(design):
                    cr = np.array(template(loc))
                    cr_mask = template(loc, return_mask=True)

                    row_mask = sample._internal_state[:, i]
                    if np.sum(row_mask) > 0:
                        dist = hamming_distance(
                            arr[row_mask][:, cr_mask], cr, return_distance=True
                        )
                        sample._internal_state[row_mask, i] = dist <= tol
                    else:
                        continue

                # keep every entry that has at least one positive
                # value in the internal state array
                ind = np.any(sample._internal_state, axis=-1)
                sample.ind_filter(ind)

            return data

        return constant_region_filter

    def vr_filter(self, where=None, loc=None, sets=None):
        """
        Filter by variable region composition.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Filters sequences based on variable region composition. Sequences containing
        amino acids or nucleotides outside the allowed monomer sets are discarded from
        all datasets in each sample.

        Args:
            where (str): Dataset to evaluate for filtering. Must be 'dna' or 'pep'.
                Variable regions are checked in this dataset, but filtering is applied
                across all associated data in each sample.

            loc (list): List of integers specifying which variable region indices to
                validate. Must reference variable regions in the library design.

            sets (list): List of integers specifying which monomer subsets to validate.
                Each integer corresponds to a monomer set defined in the library design
                configuration. Only specified sets are checked; others are ignored.

        Returns:
            callable: Operation that accepts a Data object, applies variable region filtering, and returns the modified Data object.

        Example:
            Consider a library design with the following structure::

                        seq:  ACDEF11133211AWVFRTQ12345YTPPK
                     region:  [-0-][---1--][--2--][-3-][-4-]
                is_variable:  False  True   False True False

            Region 1 contains variable positions with monomer sets 1, 2, and 3.
            The library design defines which amino acids are allowed for each set.

            >>> #validate only sets 1 and 3 in region 1 (set 2 is not checked)
            >>> vr_filt = C.fastq_parser.vr_filter(where='pep', loc=[1], sets=[1, 3])
            >>> data = vr_filt(data)
            >>>
            >>> #this would raise an error (region 2 is constant, not variable)
            >>> vr_filt = C.fastq_parser.vr_filter(where='pep', loc=[2], sets=[1])
        """
        self._where_check(where)
        design = self._infer_design(where)
        self._loc_check(loc, design)

        if not isinstance(sets, list):
            msg = f"variable_region_filter op expected to receive a list of monomer subsets to parse; received: {type(sets)}"
            self.logger.error(msg)
            raise ValueError(msg)

        allowed = set(design.monomers.keys())
        passed = set(sets)
        if not passed.issubset(allowed):
            msg = "Specified variable region sets for <variable_region_filter> op must constitute a subset of library design monomers."
            self.logger.error(msg)
            raise AssertionError(msg)

        if not np.all(design.is_vr[loc]):
            msg = "<variable_region_filter> expected a list of variable regions to operate on; some of the specified locations point to constant regions."
            self.logger.error(msg)
            raise AssertionError(msg)

        def variable_region_filter(data):
            for sample in data:
                arr = self._cast(sample[where], "2d")
                self._empty_array_check(arr, inspect.stack()[0][3])

                # first things first: temporarily expand the internal
                # state array by one dimension; will collapse back at the end
                sample._internal_state = np.repeat(
                    sample._internal_state[:, :, np.newaxis], len(sets), axis=-1
                )

                for i, template in enumerate(design):
                    # use internal state to figure out which
                    # entries are worth focusing on
                    row_mask = sample._internal_state[:, i, 0]

                    for j, subset in enumerate(sets):
                        # work out column-wise mask
                        col_mask = np.array(template(loc, return_mask=True))
                        col_mask = col_mask[np.array(template(loc)) == subset]

                        # get the matching array: check whether entries are all in the corresponding monomer subset
                        match = np.isin(
                            arr[row_mask][:, col_mask], design.monomers[subset]
                        )

                        # np.isin flattens the array, so it needs to be reshaped back
                        match = match.reshape(arr[row_mask][:, col_mask].shape)

                        # the entry is taken only if everything matches
                        sample._internal_state[row_mask, i, j] = np.all(match, axis=1)

                # reduce along the subset axis to return
                # internal state array in its original form
                sample._internal_state = np.all(sample._internal_state, axis=-1)

                # keep every entry that has at least one positive
                # value in the internal state array
                ind = np.any(sample._internal_state, axis=-1)
                sample.ind_filter(ind)

            return data

        return variable_region_filter

    def filt_ambiguous(self, where=None):
        """
        Filter sequences containing ambiguous tokens.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Removes sequences containing ambiguous characters. For DNA, this includes 'N'
        nucleotides from uncertain base calls. For peptides, this includes any amino
        acids outside the translation table specification (such as '+' tokens stemming
        from ambiguous codons).

        Args:
            where (str): Dataset to evaluate for filtering. Must be 'dna' or 'pep'.
                Ambiguous tokens are checked in this dataset, but filtering is applied
                across all associated data in each sample.

        Returns:
            callable: Operation that accepts a Data object, applies ambiguous token filtering, and returns the modified Data object.

        Example:
            >>> ambig_filt = C.fastq_parser.filt_ambiguous(where='pep')
            >>> data = ambig_filt(data)
        """
        self._where_check(where)
        allowed_monomers = self._infer_alphabet(where, alphabet=None)

        def filter_ambiguous(data):
            for sample in data:
                arr = self._cast(sample[where], "2d")
                self._empty_array_check(arr, inspect.stack()[0][3])

                # perform the check; a little annoying because pads are also technically not allowed
                ind = np.isin(arr, allowed_monomers).reshape(arr.shape)
                ind = np.sum(ind, axis=1) == self._L_summary(arr)

                # filter the sample
                sample.ind_filter(ind)

            return data

        return filter_ambiguous

    def drop_data(self, where=None):
        """
        Drop specified datasets from samples.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Removes the specified dataset from all samples to free memory or simplify
        downstream processing. Useful when certain data types are no longer needed
        in the pipeline.

        Args:
            where (str): Dataset to drop. Must be 'dna', 'pep', or 'Q' (quality scores).

        Returns:
            callable: Operation that accepts a Data object, drops the specified dataset, and returns the modified Data object.

        Example:
            >>> drop_op = C.fastq_parser.drop_data(where='dna')
            >>> data = drop_op(data)
        """
        if where not in ("pep", "dna", "Q"):
            msg = f"Invalid argument passed to <drop_dataset> op. Expected where = any of ('pep', 'dna', 'Q'); received: {where}"
            self.logger.error(msg)
            raise ValueError(msg)

        def drop_dataset(data):
            for sample in data:
                sample.drop(where)
            return data

        return drop_dataset

    def q_score_filt(self, minQ=None, avgQ=None, loc=None):
        """
        Filter data by Phred quality scores.

        .. note::
           **Pipeline Operation** — Returns a callable for use in processing pipelines.

        Filters sequences based on Phred quality scores (Q). Reads that fail the
        specified quality criteria are discarded from all datasets in each sample.

        Exactly one of ``minQ`` or ``avgQ`` must be provided:

        - If ``minQ`` is given, a read passes only if **all** relevant Q scores are
          greater than or equal to ``minQ``.
        - If ``avgQ`` is given, a read passes only if the **average** Q score across
          the relevant regions is greater than or equal to ``avgQ``.

        The optional ``loc`` argument restricts filtering to specific regions
        (e.g., primer sites, amplicons). If omitted, filtering is applied across all
        regions.

        Args:
            minQ (int, optional): Minimum quality score threshold. All Q scores in
                the specified regions must be >= this value. Mutually exclusive
                with ``avgQ``.
                
            avgQ (int, optional): Average quality score threshold. The mean Q score
                across the specified regions must be >= this value. Mutually
                exclusive with ``minQ``.
                
            loc (list[int], optional): List of region indices to evaluate for
                quality filtering. If ``None``, all regions are considered.

        Returns:
            callable: A pipeline operation that accepts a ``Data`` object, applies quality score filtering according to the specified criteria, and returns the modified ``Data`` object.

        Example:
            >>> # require all Q scores in regions 0, 1, 2 to be >= 30
            >>> q_filt = C.fastq_parser.q_score_filt(minQ=30, loc=[0, 1, 2])
            >>> data = q_filt(data)

            >>> # alternatively, require average Q score across all regions >= 35
            >>> q_filt = C.fastq_parser.q_score_filt(avgQ=35)
            >>> data = q_filt(data)
        """
        if (minQ is None) == (avgQ is None):
            msg = "<Q_score_filter>: exactly one of 'minQ' or 'avgQ' argument must be specified!"
            self.logger.error(msg)
            raise ValueError(msg)

        thresh, kind = (
            (minQ, "min") if minQ is not None else (avgQ, "avg")
        )

        if not isinstance(thresh, int):
            msg = f"<Q_score_filter> op expected to receive argument minQ/avgQ as as int; received: {type(thresh)}"
            self.logger.error(msg)
            raise ValueError(msg)

        if loc is not None:    
            self._infer_design("dna")
            self._loc_check(loc, self.D_design)

        def q_score_filter(data):
            import numpy.ma as ma
            for sample in data:
                arr = sample.civilized_Q
                self._empty_array_check(arr, inspect.stack()[0][3])

                if loc is None:  
                    if kind == 'min':
                        ind = np.all(ma.array(arr, mask=arr==0) >= thresh, axis=1).data
                        
                    elif kind == 'avg':
                        ind = (ma.array(arr, mask=arr==0).mean(axis=1) >= thresh).data
                        
                    sample.ind_filter(ind)
                    
                else:
                    for i, template in enumerate(self.D_design):
                        row_mask = sample._internal_state[:, i]
                        col_mask = template(loc, return_mask=True)
                        
                        if kind == 'min':
                            sample._internal_state[row_mask, i] = np.all(
                                arr[row_mask][:, col_mask] >= thresh, axis=1
                            )
                            
                        elif kind == 'avg':
                            sample._internal_state[row_mask, i] = np.mean(
                                arr[row_mask][:, col_mask], axis=1) >= thresh                         

                    # keep every entry that has at least one positive
                    # value in the internal state array
                    ind = np.any(sample._internal_state, axis=-1)
                    sample.ind_filter(ind)

            return data

        return q_score_filter

    def fetch_at(self, where=None, loc=None):
        """
        Extract specified regions from sequences.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Extracts and retains only the specified regions from sequences in a dataset
        specified by 'where', discarding other regions. For DNA datasets, quality
        scores are also truncated to match. Collapses the sample's internal state
        after extraction.

        Args:
            where (str): Dataset to extract from. Must be 'dna' or 'pep'.

            loc (list): List of integers specifying which regions to extract and retain.

        Returns:
            callable: Operation that accepts a Data object, extracts specified regions, and returns the modified Data object with truncated sequences.

        Example:
            >>> #truncate peptide sequences to contain only regions 1 and 3:
            >>> fetch_op = C.fastq_parser.fetch_at(where='pep', loc=[1, 3])
            >>> data = fetch_op(data)
        """
        self._where_check(where)
        design = self._infer_design(where)
        self._loc_check(loc, design)

        def _fetch_region(arr, sample, design, loc):
            # initialize the array to hold the results
            max_len = self._find_max_len(design, loc)
            result = np.zeros((arr.shape[0], max_len), dtype=arr.dtype)

            for i, template in enumerate(design):
                col_mask = template(loc, return_mask=True)
                row_mask = sample._internal_state[:, i]

                result[row_mask, : len(col_mask)] = arr[row_mask][:, col_mask]
            return result

        def fetch_region(data):
            for sample in data:
                if not sample._is_collapsed():
                    msg = f"<fetch_region> op will collapse sample {sample.name}'s internal state"
                    self.logger.debug(msg)
                    sample._collapse_internal_state()

                # has to be done manually, because in this case dna and pep args
                # are not equivalent: where='dna' also means Q scores need to be
                # truncated
                if where == "pep":
                    arr = self._cast(sample.pep, "2d")
                    self._empty_array_check(arr, inspect.stack()[0][3])
                    result = _fetch_region(arr, sample, design, loc)
                    sample.pep = result

                if where == "dna":
                    arr = self._cast(sample.dna, "2d")
                    self._empty_array_check(arr, inspect.stack()[0][3])
                    result = _fetch_region(arr, sample, design, loc)
                    sample.dna = result

                    arr = self._cast(sample.Q, "2d")
                    self._empty_array_check(arr, inspect.stack()[0][3])
                    result = _fetch_region(arr, sample, design, loc)
                    sample.Q = result

            # reindex the library design accordingly so that the downstream ops
            # can still be called with originally defined loc pointers
            design.truncate_and_reindex(loc)

            return data

        return fetch_region

    def unpad(self):
        """
        Remove padding from sequence arrays.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Removes padding columns from DNA, peptide, and quality score arrays. Columns
        where every value is a padding token are removed, reducing memory usage and
        simplifying downstream analysis.

        Returns:
            callable: Operation that accepts a Data object, removes padding, and returns the modified Data object.

        Example:
            >>> unpad_op = C.fastq_parser.unpad()
            >>> data = unpad_op(data)
        """

        def unpad_data(data):
            for sample in data:
                sample.unpad_arrays()

            return data

        return unpad_data

    def demultiplex_sample_barcodes(self, barcode_loc=None, barcode_tol=None):
        """
        Demultiplex samples based on DNA barcode sequences.
        
        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.
        
        Splits a multiplexed sample into separate sub-samples according to DNA barcode
        sequences present in specified regions. Each read is assigned to the sample
        whose barcode it most closely matches (within tolerance). Reads that match
        multiple barcodes are assigned to the first matching barcode with a warning.
        
        The operation requires a sample barcode mapping to be specified in the
        config file (``sample_barcodes`` in FastqParserConfig), which defines
        the relationship between sample names and their corresponding barcode
        sequences.
        
        This operation creates new Sample objects for each barcode, named using the
        original sample name concatenated with the barcode name. The original Data
        object is replaced with a new Data instance containing the demultiplexed
        samples.
        
        Args:
            barcode_loc (list): List of integers specifying the DNA region(s)
                containing the barcode sequence. Must reference positions in the
                library design.
                
            barcode_tol (int): Maximum Hamming distance allowed between observed and
                expected barcodes. Reads exceeding this threshold are not assigned
                to any barcode and are discarded.
        
        Returns:
            callable: Operation that accepts a Data object, performs sample
                demultiplexing based on barcodes, and returns a new Data object
                containing the demultiplexed samples.
        
        Note:
            Requires a sample barcode mapping to be specified in the config
            file (``sample_barcodes``) as a dictionary mapping bacrodes names to
            their sequences.
        
        Example:
            >>> # demultiplex using barcode in region 2, allowing up to 1 mismatch
            >>> demux_op = C.fastq_parser.demultiplex_sample_barcodes(
            ...     barcode_loc=[2], barcode_tol=1
            ... )
            >>> demultiplexed_data = demux_op(data)
        """
        self._infer_design("dna")
        self._loc_check(barcode_loc, self.D_design)
    
        if not isinstance(barcode_tol, int):
            msg = f"<demultiplex_sample_barcodes> expected to receive parameter barcode_tol as as int; received: {type(barcode_tol)}"
            self.logger.error(msg)
            raise ValueError(msg)

        if not hasattr(self, "sample_barcodes"):
            msg = "<demultiplex_sample_barcodes> lacks the necessary barcode information ('sample_barcodes' in FastqParserConfig)!"
            self.logger.error(msg)
            raise ValueError(msg)

        def split_by_barcode(data):
            all_barcodes = list(self.sample_barcodes.values())
            barcode_names = list(self.sample_barcodes.keys())

            demultiplexed_samples = []
            for sample in data:
                self._empty_array_check(sample.dna, inspect.stack()[0][3])                

                shape = (len(all_barcodes),) + sample._internal_state.shape
                bc_match = np.zeros(shape, dtype=bool)

                # iterate over all barcodes
                for b, bar in enumerate(all_barcodes):
                    min_match = len(bar) - barcode_tol
                
                    for i, template in enumerate(self.D_design):
                        # use internal state to figure out which
                        # entries are worth focusing on
                        row_mask = sample._internal_state[:, i]
                        row_idx = np.where(row_mask)[0]
                        col_idx = np.array(template(barcode_loc, return_mask=True))
                        local_bc = self._cast(sample.dna, "2d")[np.ix_(row_idx, col_idx)]
    
                        mask = (local_bc == bar).sum(axis=1) >= min_match
                        bc_match[b, row_idx[mask], i] = 1
                
                # after barcode iterations are done, we can reduce along the last
                # axis: if any of the designs match a barcode, then the whole read  
                # matches the barcode
                bc_match = bc_match.any(axis=-1)
                
                # if a read matches multiple barcodes, raise a warning and 
                # collapse the state to one read = one barcode (use the first match)
                
                # because of this possibility, we don't immediately create new samples
                # while still iterating
                if np.any(bc_match.sum(axis=0) > 1):
                    msg = f"<demultiplex_sample_barcodes>: at least some reads in sample {sample.name} match several barcodes simultaneously! The first matched barcode will be used for demultiplexing. . ."
                    self.logger.warning(msg)
                    bc_match = bc_match & (bc_match.cumsum(axis=1) == 1)                        
                
                for b in range(len(all_barcodes)):
                    bc_sub_sample_name = f'{sample.name}_{barcode_names[b]}'
                    bc_sub_sample = sample.ind_filter(bc_match[b], 
                                                      in_place=False,
                                                      new_name=bc_sub_sample_name
                                                     )
                    
                    demultiplexed_samples.append(bc_sub_sample)
                
            return Data(samples=demultiplexed_samples)
        return split_by_barcode


    def demultiplex_aa_barcodes(self, barcode_loc=None, barcode_tol=None):
        """
        Demultiplex degenerate amino acids using DNA barcodes.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Resolves ambiguous amino acid identities in non-proteinogenic saturation
        mutagenesis workflows by reading DNA barcodes. The procedure is first described
        in https://doi.org/10.1073/pnas.1809901115

        Briefly, during npAA saturation mutagenesis, several peptide libraries
        containing different genetic codes are prepared and analyzed together.
        For example, in sublibrary 1 an ``AUG`` codon may encode Met, in
        sublibrary 2 some other amino acid (e.g., N-Me-Gly), in sublibrary 3 -
        β-Ala, and so on. Therefore, the nature of the amino acid encoded by
        the AUG codon in this example is ambiguous. When the sublibraries are
        prepared, usually their mRNA are barcoded in the UTR region to encode
        the nature of such degenerate amino acids. This barcode can be read out by
        NGS to determine the nature of the degenerate amino acid in any given
        peptide.

        This operation requires a "barcode translation table" to be specified in the
        config file as part of the "constants" declaration. The barcode translation 
        table is a dictionary which denotes the mapping between the barcodes and the 
        degenerate amino acids. It looks like this::

            {"CACGAT":
              {"$": "f"}
            "ATGTCG":
              {"$": "g"}
            "AGGCTT":
              {"$": "c"}
            }

        Here, ``CACGAT`` etc are mRNA barcodes, ``$`` is the encoding used for the
        amino acid to be demultiplexed (translation_table in the config file
        will have to map the ``AUG`` codon to some placeholder token like ``$`` - since
        during translation, it's impossible to tell what actual amino acid should
        correspond to ``AUG``; the ``AUG`` codon is used here as an example, it can be any
        other codon), and ``f, g, c`` are the specific npAAs corresponding to these
        barcodes. More than 1 amino acid can be  demultiplexed at a time.

        In the example above, this operation will do a ``$`` -> ``f`` swapping for all
        reads containing barcode ``CACGAT`` in the DNA region specified by ``barcode_loc``.

        The read may mismatch the specified barcode (``CACGAT``) by no more than
        ``barcode_tol`` (maximum Hamming distance). Analogous procedure will be
        carried out for all other barcodes indicated in the barcode translation
        table.

        This operation collapses the sample's internal state.

        Args:
            barcode_loc (list): List of integers specifying the DNA region containing
                the barcode sequence.

            barcode_tol (int): Maximum Hamming distance allowed between observed and
                expected barcodes. Sequences exceeding this threshold are not
                demultiplexed.

        Returns:
            callable: Operation that accepts a Data object, performs barcode-based
                demultiplexing, and returns the modified Data object with resolved
                amino acid identities.

        Note:
            Requires a barcode translation table to be specified in the configuration
            file, mapping barcode sequences to amino acid substitutions.

        Example:
            >>> demux_op = C.fastq_parser.demultiplex_aa_barcodes(
            ...     barcode_loc=[3], barcode_tol=2
            ... )
            >>> data = demux_op(data)
        """
        self._infer_design("dna")
        self._loc_check(barcode_loc, self.D_design)
    
        if not isinstance(barcode_tol, int):
            msg = f"<demultiplex_aa_barcodes> expected to receive parameter barcode_tol as as int; received: {type(barcode_tol)}"
            self.logger.error(msg)
            raise ValueError(msg)

        if not hasattr(self.constants, "barcode_table"):
            msg = "<demultiplex_aa_barcodes>: a barcode translation table was not specified in the config. . ."
            self.logger.error(msg)
            raise ValueError(msg)

        def translate_aa_barcodes(data):
            T = self.constants.barcode_table
            allowed_barcodes = list(T.keys())

            # the 4-for loop is a bit out there, but whatchu gonna do
            for sample in data:
                self._empty_array_check(sample.dna, inspect.stack()[0][3])
                self._empty_array_check(sample.pep, inspect.stack()[0][3])

                if not sample._is_collapsed():
                    msg = f"<translate_aa_barcodes> op will collapse sample {sample.name}'s internal state"
                    self.logger.debug(msg)
                    sample._collapse_internal_state()

                for i, template in enumerate(self.D_design):
                    # use internal state to figure out which
                    # entries are worth focusing on
                    row_mask = sample._internal_state[:, i]
                    row_idx = np.where(row_mask)[0]
                    col_idx = np.array(template(barcode_loc, return_mask=True))
                    barcode = self._cast(sample.dna, "2d")[np.ix_(row_idx, col_idx)]

                    # iterate over allowed barcodes
                    for bc in allowed_barcodes:
                        min_match = len(bc) - barcode_tol
                        mask = (barcode == bc).sum(axis=1) >= min_match

                        if not np.any(mask):
                            continue

                        for degenerate_aa, new_aa in T[bc].items():
                            # code below is to modify in-place
                            row_inds = row_idx[mask]
                            aa_mask = (
                                self._cast(sample.pep, "2d")[row_inds] == degenerate_aa
                            )

                            row, col = np.where(aa_mask)

                            arr = self._cast(sample.pep, "2d")
                            arr[row_inds[row], col] = new_aa
                            sample.pep = arr

            return data

        return translate_aa_barcodes

    def save(self, where=None, fmt=None):
        """
        Save datasets to a file.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Writes the specified dataset from each sample to a file in the parser output
        directory. Files are named using the sample name and timestamp. Does not
        transform the Data object.

        Args:
            where (str): Dataset to save. Must be 'dna' or 'pep'.

            fmt (str): Output file format. Must be 'npy', 'csv', or 'fasta'.

        Returns:
            callable: Operation that accepts a Data object, saves the specified data, and returns the unmodified Data object.

        Example:
            >>> #save peptide datasets as .fasta files
            >>> save_op = C.fastq_parser.save(where='pep', fmt='fasta')
            >>> data = save_op(data)
        """
        if fmt not in ("npy", "csv", "fasta"):
            msg = f"<save_data> op received invalid fmt argument. Acceted any of ('npy', 'csv', 'fasta'); received: {fmt}"
            self.logger.error(msg)
            raise ValueError(msg)

        self._where_check(where)

        def _writer(arr, fmt, path):
            if fmt == "npy":
                np.save(path + ".npy", arr)
                return

            arr_1d = [x.decode("ascii") for x in arr]

            if fmt == "csv":
                with open(path + ".csv", "w") as f:
                    for seq in arr_1d:
                        f.write(f"{seq},\n")
                return

            if fmt == "fasta":
                with open(path + ".fasta", "w") as f:
                    for i, seq in enumerate(arr_1d):
                        f.write(f">sequence_{i}\n")
                        f.write(f"{seq}\n")
                return

        def save_data(data):
            self._prepare_destinations(root=self.dirs.parser_out, data=data)
            for sample in data:
                arr = self._cast(sample[where], "1d")
                self._empty_array_check(arr, inspect.stack()[0][3])

                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f"{sample.name}_{where}" + self.timestamp
                path = os.path.join(destination, fname)

                _writer(arr, fmt, path)

            return data

        return save_data

    def count_summary(self, where=None, top_n=None, fmt=None):
        """
        Generate sequence count summary for each sample.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Counts the occurrences of each unique sequence in the specified dataset and
        writes results to a file in the parser output directory. Does not transform
        the Data object.

        Args:
            where (str): Dataset to analyze. Must be 'dna' or 'pep'.

            top_n (int, optional): If specified, only the top N most abundant sequences
                are included in the output. If None, all unique sequences are reported.

            fmt (str): Output file format. Must be 'csv' or 'fasta'.

        Returns:
            callable: Operation that accepts a Data object, generates count summaries, and returns the unmodified Data object.

        Example:
            >>> #save top 100 peptides by count from each sample as .csv files
            >>> count_op = C.fastq_parser.count_summary(where='pep', top_n=100, fmt='csv')
            >>> data = count_op(data)
        """
        self._where_check(where)
        if fmt not in ("csv", "fasta"):
            msg = f"<fastq_count_summary> op received invalid fmt argument. Accepted any of ('csv', 'fasta'); received: {fmt}"
            self.logger.error(msg)
            raise ValueError(msg)

        if top_n is not None:
            if not isinstance(top_n, int):
                msg = f"<fastq_count_summary> op expected to receive parameter top_n as as int; received: {type(top_n)}"
                self.logger.error(msg)
                raise ValueError(msg)

        def _writer(sample, og_ind, counts, tot_count, fmt, path):
            if fmt == "csv":
                df = pd.DataFrame()
                df["Index"] = np.arange(counts.size) + 1

                if sample.pep.ndim:
                    df["Peptide"] = [
                        x.decode("ascii") for x in self._cast(sample.pep[og_ind], "1d")
                    ]
                if sample.dna.ndim:
                    df["DNA"] = [
                        x.decode("ascii") for x in self._cast(sample.dna[og_ind], "1d")
                    ]

                df[f"{where} count"] = counts
                df["Dataset %"] = 100 * df[f"{where} count"] / tot_count

                df.to_csv(path + ".csv", sep=",", index=False)

            if fmt == "fasta":
                arr_1d = [
                    x.decode("ascii") for x in self._cast(sample[where][og_ind], "1d")
                ]

                with open(path + ".fasta", "w") as f:
                    for i, seq in enumerate(arr_1d):
                        f.write(f">seq_{i + 1}_count_{counts[i]}\n")
                        f.write(f"{seq}\n")
                return

        def fastq_count_summary(data):
            self._prepare_destinations(root=self.dirs.parser_out, data=data)
            from clibas.misc import sorted_count

            for sample in data:
                self._empty_array_check(sample[where], inspect.stack()[0][3])
                _, og_ind, counts = sorted_count(
                    sample[where], top_n=top_n, return_index=True
                )
                tot_count = sample.size
                destination = os.path.join(self.dirs.parser_out, sample.name)
                fname = f"{sample.name}_{where}_count_summary" + self.timestamp
                path = os.path.join(destination, fname)

                _writer(sample, og_ind, counts, tot_count, fmt, path)

            return data

        return fastq_count_summary

    def library_design_match(self, where=None):
        """
        Summarize sequence matches to library design templates.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Analyzes how many sequences from each sample match each library design template.
        Results are written to a CSV file showing the distribution of sequences across
        different library designs. Does not transform the Data object.

        Args:
            where (str): Dataset to analyze. Must be 'dna' or 'pep'.

        Returns:
            callable: Operation that accepts a Data object, generates library design match summary, and returns the unmodified Data object.

        Example:
            >>> match_op = C.fastq_parser.library_design_match(where='pep')
            >>> data = match_op(data)
        """
        self._where_check(where)
        self._infer_design(where)

        def library_design_summary(data):
            design = self._infer_design(where)
            if not isinstance(data, Data):
                msg = f"<library_design_match_analysis> op expected data as Data type; received: {type(data)}"
                self.logger.error(msg)
                raise TypeError(msg)

            # summarize straight into a pandas dataframe
            sample_names = [sample.name for sample in data]
            templates = [template.lib_seq.decode("ascii") for template in design]

            df = pd.DataFrame(index=sample_names, columns=templates)

            # this op is really just an axis=0-wide sum of the internal states
            # the only difficulty is that internal states may not be collapsed
            # at this point, in which case results won't make much sense; raise
            # a warning                
            for sample in data:
                
                if not sample._is_collapsed():
                    msg = f"<library_design_summary> encountered a sample with {where} dataset entries matching multiple specified library templates; this will be reflected in the results. . ."
                    self.logger.warning(msg)
                
                df.loc[sample.name] = np.sum(sample._internal_state, axis=0)

            self._prepare_destinations(root=self.dirs.parser_out)
            fname = f"{self.exp_name}_by_template_breakdown" + self.timestamp
            path = os.path.join(self.dirs.parser_out, fname)
            df.to_csv(path + ".csv", sep=",")

            return data

        return library_design_summary

    def dataset_wide_count_summary(self, where=None, top_n=None):
        """
        Generate merged count summary across all samples.

        .. note::
           **Pipeline Operation** - Returns a callable for use in processing pipelines.

        Counts unique sequences across all samples and merges results into a single
        spreadsheet showing counts and percentages for each sample. Two CSV files are
        generated: one with absolute counts and one with percentages. Does not transform
        the Data object.

        Args:
            where (str): Dataset to analyze. Must be 'dna' or 'pep'.

            top_n (int, optional): If specified, only the top N most abundant sequences
                across all samples are included. If None, all unique sequences are reported.

        Returns:
            callable: Operation that accepts a Data object, generates dataset-wide count summary, and returns the unmodified Data object.

        Note:
            Requires at least 2 samples in the dataset. Single-sample datasets will
            log a warning and skip this operation.

        Example:
            >>> #summarize top 500 peptides (aggregate across all samples)
            >>> summary_op = C.fastq_parser.dataset_wide_count_summary(where='pep', top_n=500)
            >>> data = summary_op(data)
        """
        self._where_check(where)

        if top_n is not None:
            if not isinstance(top_n, int):
                msg = f"<dataset_wide_summary> op expected to receive parameter top_n as as int; received: {type(top_n)}"
                self.logger.error(msg)
                raise ValueError(msg)

        def joint_dataset_count_summary(data):
            if data.size < 2:
                msg = f"<joint_dataset_summary>: the dataset size is too small for a joint summary; data contains only {data.size} samples. . ."
                self.logger.warning(msg)
                return data

            from clibas.misc import sorted_count

            # count everything first
            dfs = []
            for sample in data:
                arr = self._cast(sample[where], "1d")
                self._empty_array_check(arr, inspect.stack()[0][3])
                
                X, C = sorted_count(arr)
                dfs.append(
                    {
                        f"{where} sequence": [s.decode("ascii") for s in X],
                        sample.name: C,
                    }
                )

            # merge and save
            M = pd.DataFrame.from_dict(dfs[0])
            for df in dfs[1:]:
                M = M.merge(
                    pd.DataFrame.from_dict(df), how="outer", on=f"{where} sequence"
                )
            M = M.fillna(0)
            
            # sort the dataframe: cc: counts columns
            cc = M.select_dtypes(include="number").columns
            M["Total counts"] = M[cc].sum(axis=1)
            M = M.sort_values(by="Total counts", ascending=False)
            
            M_p = M.copy() # percentages dataframe
            sizes = M[cc].sum(axis=0)
            M_p[cc] = 100 * M_p[cc].div(sizes)
            M_p["Fraction of total"] = M["Total counts"] / sizes.sum()

            if top_n is not None:
                M = M.head(top_n)
                M_p = M_p.head(top_n)

            # save the results
            self._prepare_destinations(root=self.dirs.parser_out)
            t = self.timestamp
            fname = f"{self.exp_name}_joint_dataset_count_summary" + t
            path = os.path.join(self.dirs.parser_out, fname)
            M.to_csv(path + ".csv", sep=",", index=False)

            # save an analogous percentage spreadsheet
            # M[columns] = 100 * M[columns].div(sizes)
            # M = M.rename(columns={"Total counts": "Fraction of total"})
            
            fname = f"{self.exp_name}_joint_dataset_percentage_summary" + t
            path = os.path.join(self.dirs.parser_out, fname)
            M_p.to_csv(path + ".csv", sep=",", index=False)

            return data

        return joint_dataset_count_summary
