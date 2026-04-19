"""
Library design specification and management.

Provides Template and LibraryDesign classes for defining DNA and peptide
library structures with variable and constant regions. Used internally by
parsers for sequence validation and region extraction.
"""


# helper functions that mimic string methods for byte-strings
def isdigit_b(b):
    return all(48 <= byte <= 57 for byte in b) and len(b) > 0


def b_iter(b):
    for x in b:
        yield bytes([x])


def b_enumerate(b):
    for idx, x in enumerate(b):
        yield (idx, bytes([x]))


import numpy as np


class Template:
    """
    Single library design template specification.

    Defines the structure of a single library template including variable
    regions (encoded as numerals 0-9) and constant regions (encoded as letters).
    Internal class used by LibraryDesign.

    Args:
        lib_seq (bytes): Template sequence (e.g., b'ACDEF11133211AWVFRTQ').

        monomers (dict): Mapping of numerals to allowed monomers.

        lib_type (str): 'dna' or 'pep'.

        val_monomer (tuple): Valid monomers for validation.

    Note:
        Not typically instantiated directly. Use LibraryDesign instead.
    """

    def __init__(self, lib_seq="", monomers={}, lib_type=None, val_monomer=None):
        self.lib_seq = lib_seq
        self.monomers = monomers
        self.lib_type = lib_type

        self._typecheck(val_monomer)
        self._build()
        return

    def _typecheck(self, val_monomer):
        if not isinstance(self.lib_seq, bytes):
            raise ValueError(
                "Library design must be specified as a byte-string (dtype=bytes). . ."
            )

        # make sure that monomer classes match what's specified by the library sequence
        lib_seq_monomer_types = set(
            int(x) for x in b_iter(self.lib_seq) if isdigit_b(x)
        )
        specified_monomer_types = set(self.monomers.keys())

        if not lib_seq_monomer_types.issubset(specified_monomer_types):
            raise ValueError(
                "Variable region monomer names must match what's specified in monomers. . ."
            )

        if self.lib_type != "dna" and self.lib_type != "pep":
            raise ValueError(
                'Library type can only accept either "dna" or "pep" as valid values'
            )

        if val_monomer is not None:
            for m in b_iter(self.lib_seq):
                if not isdigit_b(m):
                    if m not in val_monomer:
                        raise ValueError(
                            f"All library design monomers must be specified in the lookup tables! Monomer {m} was not understood. . ."
                        )

            for key in self.monomers:
                for m in self.monomers[key]:
                    if m not in val_monomer:
                        raise ValueError(
                            f"All library design monomers must be specified in the lookup tables! Monomer {m} was not understood. . ."
                        )
        return

    def _build(self):
        # L is expected length of the library sequence
        self.L = len(self.lib_seq)

        # build the key params: is_vr, region, mask
        is_vr = []
        region = []
        mask = []

        current_region = list()
        current_mask = list()

        for i, m in b_enumerate(self.lib_seq):
            if i == 0:
                is_vr.append(isdigit_b(m))

            if isdigit_b(m):
                if is_vr[-1] is True:
                    current_region.append(int(m))
                    current_mask.append(i)

                else:
                    region.append(current_region)
                    mask.append(current_mask)

                    current_region = [int(m)]
                    current_mask = [i]
                    is_vr.append(True)

            else:
                if is_vr[-1] is not True:
                    current_region.append(m)
                    current_mask.append(i)

                else:
                    region.append(current_region)
                    mask.append(current_mask)

                    current_region = [m]
                    current_mask = [i]
                    is_vr.append(False)

        region.append(current_region)
        mask.append(current_mask)

        self.is_vr = np.array(is_vr, dtype=bool)
        self.loc = np.arange(self.is_vr.size)
        self.mask = mask
        self.region = region
        return

    def __repr__(self):
        seq = b"".join(b_iter(self.lib_seq))
        return (
            f"<Template container for {seq.decode('ascii')} lib_type={self.lib_type}>"
        )

    def _fancy_index(self, arr, loc):
        out = []
        for x in loc:
            ind = np.where(self.loc == x)[0][0]
            out.extend(arr[ind])
        return out

    def __call__(self, loc, return_mask=False):
        """
        Extract regions or position masks at specified locations.

        Fancy indexing to retrieve sequence content or column indices for
        specified region locations.

        Args:
            loc (list): Region indices to extract.

            return_mask (bool): If True, returns column indices. If False,
                returns sequence content. Default is False.

        Returns:
            list: Region content or column mask depending on return_mask.
        """
        if not np.all(np.isin(loc, self.loc)):
            raise ValueError(
                "Library design: a call to a non-existent region was made. . ."
            )

        if return_mask:
            arr = self.mask
        else:
            arr = self.region

        return self._fancy_index(arr, loc)

    def truncate_and_reindex(self, loc):
        """
        Keep only specified regions and reindex template.

        Modifies template in place to retain only regions at specified
        locations. Used during 'fetch_at' operations.

        Args:
            loc (list): Region indices to keep.
        """

        def remask():
            mask = list()
            ind = 0
            for reg in self.region:
                current = list()
                for x in reg:
                    current.append(ind)
                    ind += 1
                mask.append(current)

            self.mask = mask
            self.L = ind + 1
            return

        self.is_vr = self.is_vr[loc]
        self.loc = np.array(loc)

        reg = list()
        for ind in self.loc:
            reg.append(self.region[ind])

        self.region = reg
        remask()
        return

    def rebuild(self):
        """
        Rebuild template from original sequence specification.

        Reverts template to original state after truncation.
        """
        return self._build()


class LibraryDesign:
    """
    Library design specification with multiple templates.

    Defines library structure including variable and constant regions for
    DNA or peptide libraries. Supports multiple templates with different
    variable region lengths but identical topology (same number and type
    of regions).

    Variable regions use numerals (0-9) to indicate randomized positions,
    with each numeral corresponding to a specific monomer set. Constant
    regions use standard single-letter codes.

    Args:
        templates (list): List of template sequences as bytes.

        monomers (dict): Mapping of numerals to allowed monomer lists.

        lib_type (str): 'dna' or 'pep'.

        val_monomer (tuple): Acceptable monomers to validate the constructs.
            Usually, amino acids specified in the config translation table or
            a list of acceptable bases from the config file.

    Example:
        Template structure::

                    seq:  ACDEF11133211AWVFRTQ12345YTPPK
                 region:  [-0-][---1--][--2--][-3-][-4-]
            is_variable:  False  True   False True False

        >>> lib = LibraryDesign(
        ...     templates=[b'211113GSGSGS', b'2111113GSGSGS'],
        ...     monomers={
        ...               1: [b'A', b'C', b'D'],
        ...               2: [b'M'],
        ...               3: [b'C']
        ...     },
        ...     lib_type='pep',
        ...     val_monomer=None
        ... )

    Note:
        - All templates must have identical topology
        - Variable region numerals must be defined in monomers dict
    """

    def __init__(self, templates=[], monomers={}, lib_type=None, val_monomer=None):
        self.monomers = monomers
        self.lib_type = lib_type
        self.templates = tuple(
            Template(
                lib_seq=x,
                monomers=self.monomers,
                lib_type=self.lib_type,
                val_monomer=val_monomer,
            )
            for x in templates
        )

        self.L = list(set([x.L for x in self.templates]))
        self._topology_check()

        self.loc = self.templates[0].loc
        self.is_vr = self.templates[0].is_vr
        return

    def __repr__(self):
        return f"<Library design container for {len(self.templates)} templates. lib_type={self.lib_type}>"

    def __iter__(self):
        for template in self.templates:
            yield template

    def __len__(self):
        return len(self.templates)

    def __getitem__(self, item):
        return self.templates[item]

    def _topology_check(self):
        """
        Validate that all templates share identical region topology.

        Ensures all templates have the same arrangement of variable and constant
        regions (e.g., constant-variable-constant). Checks that is_vr patterns
        match across all templates.

        Raises:
            ValueError: If templates have different topologies.
        """
        topologies = []
        for t in self.templates:
            topologies.append(tuple(t.is_vr))

        if len(set(topologies)) != 1:
            raise ValueError("All library templates should have the same topology. . .")

        return

    def truncate_and_reindex(self, loc):
        """
        Truncate all templates to specified regions.

        Args:
            loc (list): Region indices to keep.
        """
        for template in self.templates:
            template.truncate_and_reindex(loc)

        return

    def rebuild(self):
        """Rebuild all templates from original specifications."""
        for template in self.templates:
            template.rebuild()

        return
