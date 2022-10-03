# Tutorial 4: Creating a Basic Network

In this tutorial, you will learn how to set up a network with some pre-existing data, how to store it in a file, and how to assign metadata.

Starting off, let's import pickaxe_generic and set up some molecules and operators.

```python
import pickaxe_generic as pg

engine = pg.create_engine()

water = engine.mol.rdkit("O")
ethanol = engine.mol.rdkit("CCO")
acetone = engine.mol.rdkit("CC(C)=O")
butanone = engine.mol.rdkit("CCC(C)=O")
methyl_butanoate = engine.mol.rdkit("CCCC(=O)OC")
delta_valerolactone = engine.mol.rdkit("O=C1CCCCO1")
hydroxyvaleric_acid = engine.mol.rdkit("O=C(O)CCCCO")

aldol_condensation = engine.op.rdkit("[O&+0:1]=[C&+0:2]-[C&+0;H2,H3:3].[C&+0:4]=[O&+0:5]>>[*:1]=[*:2]-[*:3]=[*:4].[*:5]")
ester_hydrolysis_nonring = engine.op.rdkit("[O&+0:1]=[C&+0:2]-&!@[O&+0&H0:3].[O&+0&H2:4]>>[*:1]=[*:2]-[*:4].[*:3]")
ester_hydrolysis_ring = engine.op.rdkit("[O&+0:1]=[C&+0:2]-&@[O&+0&H0:3].[O&+0&H2:4]>>([*:1]=[*:2]-[*:4].[*:3])")
esterification = engine.op.rdkit("[O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4]>>[*:1]=[*:2]-[*:4].[*:3]")
esterification_intra = engine.op.rdkit("([O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4])>>[*:1]=[*:2]-[*:4].[*:3]")
```

## Initializing Network

The first step in setting up a new network is to get a new one from the engine.

```python
network = engine.new_network()
```

The network object provides access to the molecules, operators, and reactions of which it is constituted.

At the moment, there are no molecules inside the network.  We must add our molecules to it ourselves.  Let's try adding water.

```python
water_i = network.add_mol(water)
```
```sh
>>> water_i
0
```

We've successfully added water to the network.  When you add an object to the network, the network will return an integer representing its index inside the network.  We can retrieve the molecule from the network by its index.

```sh
>>> network.mols[water_i]
MolDatBasic('O')
```

### Notes for Programmers

The index always starts at zero and increases for every subsequent addition.  Objects cannot be removed from the network once added.  Pruning is accomplished by creating a new network with a more limited subset of items.

## Retrieving Objects From Network

This is all well and good, but what if someone else put water into the network before us, and we want to see if it's in there?  Fortunately, there are a couple simple solutions.

One way to check for water is by making your own water molecule and then checking to see if it is present.  Using the UID of the molecule also works.

```sh
>>> water in network.mols
True
>>> water.uid in network.mols
True
>>> butanone in network.mols
False
```

If you want to retrieve the molecule's index from the network, you can do so using the `.i()` method.  At the moment, only lookup by UID is supported.

```sh
>>> network.mols.i(water.uid)
0
```

To retrieve the molecule, simply use its index.

```sh
>>> network.mols[water_i]
MolDatBasic('O')
```

What happens if we try to add another water molecule to the network?

```sh
>>> network.add_mol(water)
0
```

It turns out there can only be one molecule with a particular UID within the network at a time.  If it is already present, the existing index is then returned.

You can also check the number of molecules in the network using `len()`.

```sh
>>> len(network.mols)
1
```

Interacting with operators on the network is largely the same as for molecules.

```python
ester_hydrolysis_ring_i = network.add_op(ester_hydrolysis_ring)
```
```sh
>>> print(ester_hydrolysis_ring_i)
0
```

### Reactions on the Network

Reactions work a little bit differently on the network.  Since they are a purely associative construct, they have a slightly different interface for adding them.  Before we do, let's add an additional molecule so that we can generate a reaction.

```python
network.add_mol(delta_valerolactone)
```

When we react these molecules using the ester hydrolysis operator above, we get one product, 5-hydroxyvaleric acid.

```python
hydroxyvaleric_acid = ester_hydrolysis_ring(delta_valerolactone,water)[0][0]
```
```sh
>>> hydroxyvaleric_acid
MolDatBasic('O=C(O)CCCCO')
```

Once we add the molecule to the network, we can use the `.add_rxn()` method to add the reaction.  We know implicitly the indices of the molecules added because of the order in which they were added, but we can make sure by constructing our reaction as below.

```python
water_i = network.mols.i(water.uid)
delta_valerolactone_i = network.mols.i(delta_valerolactone.uid)
ester_hydrolysis_ring_i = network.ops.i(ester_hydrolysis_ring.uid)

hydroxyvaleric_acid = ester_hydrolysis_ring(delta_valerolactone,water)[0][0]
hydroxyvaleric_acid_i = network.add_mol(hydroxyvaleric_acid)

reaction_i = network.add_rxn(
    operator=ester_hydrolysis_ring_i,
    reactants=(delta_valerolactone_i, water_i),
    products=(hydroxyvaleric_acid_i,),
)
```
```sh
>>> reaction_i
0
>>> network.rxns[0]
Reaction(operator=0, reactants=(1, 0), products=(2,))
```

The Reaction object is basically just a [named tuple](https://docs.python.org/3/library/typing.html#other-special-directives) which associates the indices of the components involved, so it takes up only a minimum of storage space.


## Connectivity of the Network

The network provides a couple of methods for exploration.  The first method is simply that `network.mols`, `network.ops`, and `network.rxns` are all iterable and sliceable, so all of their information can be displayed at once if so desired.

Another option, for example if you have an algorithm to find connections within an existing chemical network, is to list the "consumer" and "producer" reactions for a particular molecule.  This is, fittingly, done by calling the `.consumer()` and `.producer()` methods.

```sh
>>> network.consumers(0)
[0]
>>> network.producers(0)
[]
>>> network.producers(2)
[0]
```

The results shown are the indices of the reactions which produce or consume the molecule queried.  Molecule 0 (water) is consumed by the only reaction we have added so far (reaction 0).  However, it is not produced by anything.  However, molecule 2 (5-hydroxyvaleric acid) is produced by the same reaction which consumed the water.  Though this information can be easily extracted and processed from `network.rxns` when necessary, we have made it available as native functionality.

## Saving Network to a File

It is very straightforward to save a chemical network to a file using the `.save_to_file()` method.  Simply provide a filename (without extension) and an optional path, and the function will save all network information.

```python
network.save_to_file("saved_network")
```
```sh
> ls saved_network*
saved_network.pgnet
```

## Loading Network from File

Network files are version-controlled, so networks saved using an older version of Pickaxe-Generic should be loaded properly.  Loading a network file from a newer version of Pickaxe-Generic should either work properly or raise an error.  Please report bugs to the GitHub.

Networks may be loaded from files via the engine.

```python
network_loaded = engine.network_from_file("saved_network")
```
```sh
>>> print(list(network_loaded.mols))
[MolDatBasic('O'), MolDatBasic('O=C1CCCCO1'), MolDatBasic('O=C(O)CCCCO')]
```

## Metadata

The network also stores "metadata," which consists of [Mappings](https://docs.python.org/3/library/stdtypes.html#mapping-types-dict) associated with every object in the network, including reactions.  Using metadata, you can store a key-value pair associated with a particular molecule, like enthalpy, or a cumulative cost calculation based on a particular network expansion method.

To access metadata assigned to an object, we must first add metadata to the object.  The easiest way to do this is when the object is added to the network, via the `meta` argument.  As seen below, this can also be done for objects which are already in the network, though key collisions will be resolved by overwriting with the newest value.

```python
ethanol_i = network.add_mol(ethanol, meta={"is_alcohol": True})
water_i = network.add_mol(water, meta={"is_alcohol": False, "solvent": True})
```

Metadata can be then accessed via the `.mols.meta()` method.  This can either be
targeted at a specific molecule, or by having all the molecules print their information.

```sh
>>> print(network.mols.meta(water_i))
{'is_alcohol': False, 'solvent': True}
>>> print(network.mols.meta(water_i, ["is_alcohol"]))
{'is_alcohol': False}
>>> print(network.mols.meta([ethanol_i], ["is_alcohol"]))
({'is_alcohol': True},)
>>> print(network.mols.meta([water_i, ethanol_i], keys=["is_alcohol"]))
({'is_alcohol': False}, {'is_alcohol': True})
>>> print(network.mols.meta([water_i, ethanol_i], keys=["solvent"]))
({'solvent': True}, {})
>>> print(network.mols.meta(keys=["is_alcohol"]))
({'is_alcohol': False}, {}, {}, {'is_alcohol': True})
>>> print(network.mols.meta())
({'is_alcohol': False, 'solvent': True}, {}, {}, {'is_alcohol': True})
```

