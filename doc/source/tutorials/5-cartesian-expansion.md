# Cartesian Expansion

In this tutorial, you will learn how to perform an automated network expansion using a Cartesian strategy.

In DORAnet, a Recipe is defined by the unique combination of an operator and an ordered set of reactants. Given a set of operators and a set of molecules, the Cartesian product refers to all of the Recipes which may be obtained via different combinations of operators and reactants.

A "Cartesian expansion" is the iterative method of combining all available reactants and operators to create new products, adding those products to the initial set, and performing the operation again. Depending on the iteration which produced them, the molecules in the network may be assigned different "generations."

### Math Note

Technically, the Cartesian product refers to the Cartesian product of all operators with all possible subsets of reactant molecules available. Once the incompatible subsets are filtered out, the Cartesian product as DORAnet understands it is obtained.

## Setting Up Network Expansion

The first step in an automated expansion is to create some initial reactants and operators. Here we choose hydrogenation as our initial operator.

```python
import doranet as dn

engine = dn.create_engine()

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

Now we create a network and add these components to it. Giving your operators unique names, and storing those names in metadata can make filtering and post-processing easier, so we show it here.

```python
network = engine.new_network()

for smiles in initial_reactant_smiles:
    network.add_mol(engine.mol.rdkit(smiles))

for name, smarts in operator_smarts.items():
    network.add_op(engine.op.rdkit(smarts), meta={"name": name})
```

As an example of using the name of the operator, we can list the operators by name.

```sh
>>> [x["name"] for x in network.ops.meta(keys=["name"])]
['hydrogenation of alkene/carbonyl']
```

This becomes more useful when viewing reactions (which only indicate the index of the operator).

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

It does not do any work initially, but must have its `.expand()` method called in order to perform an expansion. We can also display the molecules before and after expansion.

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
>>> pprint(list((rxn, network.ops.meta(rxn.operator,"name")) for rxn in network.rxns))
[Reaction(operator=0, reactants=(1, 0), products=(5,)),
 Reaction(operator=0, reactants=(2, 0), products=(6,)),
 Reaction(operator=0, reactants=(4, 0), products=(7,)),
 Reaction(operator=0, reactants=(7, 0), products=(8,))]
```

As you can see, all possible molecules were hydrogenated, and butadiene was hydrogenated twice to form butane.

### Developer Note

Under the hood, the Cartesian strategy is using the much more complex PriorityQueue strategy with no ranker function and a custom global hook function. Otherwise, much of the functionality is the same.

## Limiting Network Size

The example above was fairly simple. Hydrogenation is, by definition, limited by the saturation of the targeted molecules. However, with operators which create larger molecules or operate on very generic reaction sites, the number of molecules generated can quickly explode.

The cartesian expansion, like other strategies, comes with several ways to limit the size of the generated network. The first is by limiting the number of recipes which are tested.

### Limiting Number of Recipes

Using the `max_recipes` argument, we can see that limiting the number of recipes to `max_recipes=2` results in only two reactions appearing in the system.

```python
network = engine.network_from_file("5-cartesian-expansion-initial")
strat = engine.strat.cartesian(network)
strat.expand(max_recipes=2)
```

```sh
>>> pprint(list(network.rxns))
[Reaction(operator=0, reactants=(1, 0), products=(5,)),
 Reaction(operator=0, reactants=(2, 0), products=(6,))]
```

This method is effective for limiting the overall size of the network in an absolute sense. However, of the reactions which are available at the start of the program, some will inevitably be prioritized over others.

### Limiting Number of Cartesian Products

The Cartesian strategy will naturally perform multiple iterations over the sets of molecules and operators. Limiting these iterations provides a way to limit network size, while not favoring any particular combination of reactants and operators during an iteration. However, the cost is that if there are very many Recipes which are possible in an iteration, they will all be evaluated before the stopping condition is met.

Setting a limit on the number of iterations is done via the `num_iter` argument.

```python
network = engine.network_from_file("5-cartesian-expansion-initial")
strat = engine.strat.cartesian(network)
strat.expand(num_iter=1)
```

```sh
>>> pprint(list(network.rxns))
[Reaction(operator=0, reactants=(1, 0), products=(5,)),
 Reaction(operator=0, reactants=(2, 0), products=(6,)),
 Reaction(operator=0, reactants=(4, 0), products=(7,))]
```

## Takeaways

1. Each iteration, the Cartesian strategy combines all operators with their compatible molecules.
2. The network size can be limited either by limiting the number of iterations or by limiting the total number of recipes which may be tested.

Congratulations! You have finished the fifth part of the DORAnet tutorial. Proceed to the [next part](./6-filters.md) to learn how to use filters in order to restrict your network expansion to only the most relevant reactions.
