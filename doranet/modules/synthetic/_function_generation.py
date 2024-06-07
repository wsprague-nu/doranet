import pathlib
from collections.abc import Iterable

import doranet as dn


def generate_network(
    starters: str
    | pathlib.Path
    | Iterable[
        dn.interfaces.MolDatBase
        | dn.interfaces.MetaStruct[dn.interfaces.MolDatBase]
    ],
    operators: str
    | pathlib.Path
    | Iterable[
        dn.interfaces.OpDatBase
        | dn.interfaces.MetaStruct[dn.interfaces.OpDatBase]
    ],
    targets: str
    | pathlib.Path
    | Iterable[
        dn.interfaces.MolDatBase
        | dn.interfaces.MetaStruct[dn.interfaces.MolDatBase]
    ]
    | None = None,
    helpers: str
    | pathlib.Path
    | Iterable[
        dn.interfaces.MolDatBase
        | dn.interfaces.MetaStruct[dn.interfaces.MolDatBase]
    ]
    | None = None,
    max_depth: int | None = 1,
) -> dn.interfaces.ChemNetwork:
    """
    Generate a synthetic chemical network.

    Specifications for molecules and operators may be provided in the form of a
    file path or an Iterable of the relevant datastructure.

    Parameters
    ----------
    starters
        Starting set of molecules to seed the network.
    operators
        Operators to be recursively used on the set of molecules.
    targets
        Target molecules. These molecules will not be used as an argument to
        any operator, and (depending on strategy) will shape the network
        expansion. Should not overlap with the `starters` or `helpers` provided.
    helpers
        Helper molecules. These molecules can fulfill the role of reactant, but
        there must be at least one non-helper molecule used as a reactant in
        each reaction.
    max_depth
        Maximum number of steps from starter molecules (analogous to number of
        recursive iterations of Cartesian expansion).
    """
    raise NotImplementedError("This function not yet implemented")
