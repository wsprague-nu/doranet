from ._dataunit import DataUnit


class Molecule(DataUnit):
    """
    Empty interface signalling molecule status of DataUnit class.

    Classes implementing this interface manage information about a single
    molecule, allowing for memory management and lumped molecule frameworks.

    Attributes
    ----------
    blob : bytes
        Binary representation of object.
    uid : doranet.interfaces.Identifier
        Unique identifier of object.
    """

    __slots__ = ()


MolDatBase = Molecule
