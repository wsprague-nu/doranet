# Tutorial 5: Cartesian Expansion

In this tutorial, you will learn how to perform an automated network expansion using a Cartesian strategy.

In Pickaxe-Generic, a Recipe is defined by the unique combination of an operator and an ordered set of reactants.  Given a set of operators and a set of molecules, the Cartesian product refers to all of the Recipes which may be obtained via different combinations of operators and reactants.

A "Cartesian expansion" is the iterative method of combining all available reactants and operators to create new products, adding those products to the initial set, and performing the operation again.  Depending on the iteration which produced them, the molecules in the network may be assigned different "generations."

### Math Note

Technically, the Cartesian product refers to the Cartesian product of all operators with all possible subsets of reactant molecules available.  Once the incompatible subsets are filtered out, the Cartesian product as Pickaxe-Generic understands it is obtained.

## Setting Up Network Expansion

The first step in an automated expansion is to create some initial reactants and operators.  Here we choose hydrogenation as our initial operator.

```python
import pickaxe_generic as pg

engine = pg.create_engine()

initial_reactant_smiles = [
    "[H][H]", # hydrogen
    "CC=O", # acetaldehyde
    "CC(C)=O", # acetone
    "CCCO", # propanol
    "C=CC=C", # butadiene
]

operator_smarts = {
    "hydrogenation of alkene/carbonyl": "[C,O;+0:1]=[C&+0:2].[#1][#1]>>[*:1]-[*:2]"
}
```

Now we create a network and add these components to it.  Giving your operators unique names, and storing those names in metadata can make filtering and post-processing easier, so we show it here.

```python
network = engine.new_network()

for smiles in initial_reactant_smiles:
    network.add_mol(engine.mol.rdkit(smiles))

for name, smarts in operator_smarts.items():
    network.add_op(engine.op.rdkit(smarts), meta={"name": name})
```

We will save this initial network to a file, so that we may restore it for testing out new strategies and configuration options.

```python
network.save_to_file("5-cartesian-expansion-initial")
```

## Running a Basic Expansion

A strategy is initialized with some configurational elements, including the target network.

```python
network = engine.network_from_file("5-cartesian-expansion-initial")
strat = engine.strat.cartesian(network)
```

It does not do any work initially, but must have its `.expand()` method called in order to perform an expansion.  We can also display the molecules before and after expansion.

```sh
>>> from pprint import pprint
>>> pprint(list(enumerate(network.mols)))
[(0, MolDatBasic('[H][H]')),
 (1, MolDatBasic('CC=O')),
 (2, MolDatBasic('CC(C)=O')),
 (3, MolDatBasic('CCCO')),
 (4, MolDatBasic('C=CC=C'))]
```
```python
strat.expand()
```
```sh
>>> pprint(list(enumerate(network.mols)))
[(0, MolDatBasic('[H][H]')),
 (1, MolDatBasic('CC=O')),
 (2, MolDatBasic('CC(C)=O')),
 (3, MolDatBasic('CCCO')),
 (4, MolDatBasic('C=CC=C')),
 (5, MolDatBasic('CCO')),
 (6, MolDatBasic('CC(C)O')),
 (7, MolDatBasic('C=CCC')),
 (8, MolDatBasic('CCCC'))]
>>> pprint(list(enumerate(network.rxns)))
[Reaction(operator=0, reactants=(1, 0), products=(5,)),
 Reaction(operator=0, reactants=(2, 0), products=(6,)),
 Reaction(operator=0, reactants=(4, 0), products=(7,)),
 Reaction(operator=0, reactants=(7, 0), products=(8,))]
```

As you can see, all possible molecules were hydrogenated, and butadiene was hydrogenated twice to form butane.

## Limiting Network Size

The example above was fairly simple.  Hydrogenation is, by definition, limited by the saturation of the targeted molecules.  However, with operators which create larger molecules or operate on very generic reaction sites, the number of molecules generated can quickly explode.

The cartesian expansion, and indeed most expansion strategies, comes with several ways to limit the number of reactions.  The first is by limiting the number of recipes which are tested.