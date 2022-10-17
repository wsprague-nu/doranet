# Tutorial 7: Metadata

In this tutorial, you will learn how to generate metadata during network expansion.

First, create an engine and network with some reactants and initial reagents, saving the network to a file in order to run multiple experiments from the same initial state.

```python
import pickaxe_generic as pg

engine = pg.create_engine()

network = engine.new_network()

reagents = [
    "[H][H]",  # hydrogen
    "O",  # water
    "CO",  # methanol
    "CCO",  # ethanol
    "CC(O)=O",  # acetic acid
]

operator_smarts = {
    "ester_hydrolysis_nonring": "[O&+0:1]=[C&+0:2]-&!@[O&+0&H0:3].[O&+0&H2:4]>>[*:1]=[*:2]-[*:4].[*:3]",
    "ester_hydrolysis_ring": "[O&+0:1]=[C&+0:2]-&@[O&+0&H0:3].[O&+0&H2:4]>>([*:1]=[*:2]-[*:4].[*:3])",
    "esterification": "[O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4]>>[*:1]=[*:2]-[*:4].[*:3]",
    "esterification_intra": "([O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4])>>[*:1]=[*:2]-[*:4].[*:3]",
    "hydrogenation of carbonyl": "[C+0:1]=[O+0:2].[H][H]>>[*:1][*:2]",
}

for smiles in reagents:
    network.add_mol(engine.mol.rdkit(smiles))

for name, smarts in operator_smarts.items():
    network.add_op(engine.op.rdkit(smarts), meta={"name": name})

network.save_to_file("6-filters")
```

## Types of Metadata

Metadata in Pickaxe-Generic is information associated with a molecule, operator, or reaction (here referenced as "data objects").  Despite such a general definition, we can still outline some classifications of metadata which may be calculated "on-the-fly".

The first way to classify metadata is based on its "purity" by separating metadata into *inherent* and *contextual* categories.  These categories are related to the concept of a [pure function](https://en.wikipedia.org/wiki/Pure_function).  For a given data object, metadata which is always the same for that object, regardless of context, is termed *inherent* metadata.  By contrast, *contextual* metadata relies on the data object's context within the network, and therefore might change during a given network expansion.  Some examples of each of these types are given below.

* Inherent Metadata
  * Gas-phase enthalpy of formation at 298 K
  * Gas-phase entropy of formation at 298 K
  * Molecular formula
  * Molecular fingerprint
  * Tanimoto similarity of molecule to target
  * Toxicity of molecule
  * Name of operator
  * Gas-phase enthalpy of reaction at 298 K
  * Gas-phase Gibbs free energy of reaction at 298 K
* Contextual Metadata
  * Maximum atom economy to produce molecule from particular reagents, assuming side products are discarded
  * Minimum cost to produce molecule from particular reagents, assuming side products are discarded
  * Net flux of molecule (for microkinetic models)

The inherent vs. contextual metadata classification is most useful when deciding how it ought to be calculated.  However, there is another classification which is more useful when determining if metadata calculation may be achieved in parallel during network expansion (the most efficient place for it).

The second way to classify metadata is by separating it based on the information necessary to calculate it.  Metadata which requires only the information contained in a single reaction to calculate (plus some way to resolve between different values) is termed *local* metadata.  Metadata which requires more information, typically about the network as a whole, is termed *global* metadata.

All *inherent* metadata is *local*, since it can be calculated easily requiring not just a single reaction, but a single data object.  However, some of the previously defined *contextual* metadata may be calculated locally.  One example is maximum atom economy (to produce a particular molecule from particular reagents which are not necessarily its immediate precursors).  If the maximum atom economy to produce each reactant is known, then the atom economy of the products of a reaction may be calculated using only that information and the reaction's stoichiometry.  If two reactions both produce the same molecule, but different atom economy values for that molecule, then clearly the maximum between those two is the correct one.  Because a value may be calculated from an individual reaction, and there exists a function to consistently resolve conflicts between that value and others which are calculated, maximum atom economy is *local* metadata.  However, net flux for batch reactors is based on the state of the entire microkinetic network, and so cannot be calculated locally.  In fact, depending on the model used the net flux may change regularly for all molecules.

This model of *local* vs. *global* metadata may seem complex, but is required in order to determine which values may be calculated during network generation, instead of being globally updated between network generation steps.  Only *local* metadata may be calculated during network generation.

### Developer's Note

You may realize that inherent metadata does not necessarily have to be stored as metadata, since it can be recalculated given only the data object.  While this is technically true, having to recalculate complex inherent metadata may take more time than simply using the metadata framework as a cache for this information.  In addition, when running in parallel Pickaxe-Generic will not pass the data object to Recipe enumeration processes if only the metadata will suffice for the recipe filters and recipe ranking functions.  This saves on I/O, which is often the bottleneck of parallel processing.

## The Reaction Analysis Plan

## Inherent Metadata

## Contextual, Local Metadata

## Global Metadata

## Filtering on Metadata

## The Reaction Analysis Plan, Revisited

## Takeaways

