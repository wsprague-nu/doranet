"""Essential Molecule object type."""

import abc

from ._dataunit import DataUnit


class Molecule(DataUnit):
    """Interface for objects representing a molecule."""

    @abc.abstractmethod
    def get_inchikey(self) -> str:
        """
        Return nearly-unique InChIKey hash of molecule.

        Returns
        -------
        str
            InChIKey hash.
        """

    @abc.abstractmethod
    def get_inchi(self) -> str:
        """
        Return InChI string of molecule.

        Returns
        -------
        str
            InChI hash.
        """

    @abc.abstractmethod
    def get_smiles(self) -> str:
        """
        Return Canonical SMILES string of molecule.

        Returns
        -------
        str
            SMILES string.
        """

    @abc.abstractmethod
    def get_cxsmiles(self) -> str:
        """
        Return CXSMILES string of molecule.

        Returns
        -------
        str
            CXSMILES string.
        """
