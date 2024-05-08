# Global Hooks

In this tutorial, you will learn how to write, apply, and compose global hook functions in order to implement stopping criteria and global metadata updates.

First, create an engine and network with some reactants and initial reagents, saving the network to a file in order to run multiple experiments from the same initial state.

```python
import doranet as dn

engine = dn.create_engine()

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
    network.add_mol(engine.mol.rdkit(smiles), meta={"gen": 0})

for name, smarts in operator_smarts.items():
    network.add_op(engine.op.rdkit(smarts), meta={"name": name})

network.save_to_file("8-global-hooks")
```

## Stopping Criteria

When expanding the synthetic network, there may be some criterion which, when met, causes the expansion to cease. One of these may be that a certain maximum number of molecules is reached (possibly for memory reasons). Another may be that a target molecule has been produced. In each of these scenarios, a global hook function can provide the answer.

In the example below, a hook function is added to one strategy which sets a total molecule threshold at 10. When this threshold is exceeded, the network will not be expanded by further iterations. This strategy is compared to one without the global hook function.

Note also that the hook function is put in a list; it must be passed in as a sequence so that multiple hook functions can be used in the same context.

```python
from pprint import pprint

network_no_hook = engine.network_from_file("8-global-hooks")
network_with_hook = engine.network_from_file("8-global-hooks")

strat_no_hook = engine.strat.cartesian(network_no_hook)
strat_with_hook = engine.strat.cartesian(network_with_hook)

mol_limit_hook = engine.hook.max_mols(10)
gen_calc = engine.meta.generation("gen")

strat_no_hook.expand(num_iter=3, reaction_plan=gen_calc)
strat_with_hook.expand(
    num_iter=3, global_hooks=[mol_limit_hook], reaction_plan=gen_calc
)
```

```sh
>>> pprint(
...   [
...     (i, v[0].smiles, v[1])
...     for i, v in enumerate(
...       zip(network_no_hook.mols, network_no_hook.mols.meta(keys=["gen"]))
...     )
...   ]
... )
[(0, '[H][H]', {'gen': 0}),
 (1, 'O', {'gen': 0}),
 (2, 'CO', {'gen': 0}),
 (3, 'CCO', {'gen': 0}),
 (4, 'CC(=O)O', {'gen': 0}),
 (5, 'CC(O)O', {'gen': 1}),
 (6, 'COC(C)=O', {'gen': 1}),
 (7, 'CCOC(C)=O', {'gen': 1}),
 (8, 'CC(=O)OC(C)=O', {'gen': 1}),
 (9, 'CC(=O)OC(C)O', {'gen': 2}),
 (10, 'COC(C)O', {'gen': 2}),
 (11, 'CCOC(C)O', {'gen': 2}),
 (12, 'CC(O)OC(C)O', {'gen': 3}),
 (13, 'CC(=O)OC(C)OC(C)=O', {'gen': 3}),
 (14, 'COC(C)OC(C)=O', {'gen': 3}),
 (15, 'CCOC(C)OC(C)=O', {'gen': 3})]
>>> pprint(
...   [
...     (i, v[0].smiles, v[1])
...     for i, v in enumerate(
...       zip(
...         network_with_hook.mols,
...         network_with_hook.mols.meta(keys=["gen"]))
...     )
...   ]
... )
[(0, '[H][H]', {'gen': 0}),
 (1, 'O', {'gen': 0}),
 (2, 'CO', {'gen': 0}),
 (3, 'CCO', {'gen': 0}),
 (4, 'CC(=O)O', {'gen': 0}),
 (5, 'CC(O)O', {'gen': 1}),
 (6, 'COC(C)=O', {'gen': 1}),
 (7, 'CCOC(C)=O', {'gen': 1}),
 (8, 'CC(=O)OC(C)=O', {'gen': 1}),
 (9, 'CC(=O)OC(C)O', {'gen': 2}),
 (10, 'COC(C)O', {'gen': 2}),
 (11, 'CCOC(C)O', {'gen': 2})]
```

Notice that the second run, with the maximum molecules hook, stopped at an earlier generation than the first. It has more than 10 molecules, but the hook is not called between iterations (here generations), so it has no way of stopping the additional molecules from being generated.

If this behavior is confusing, check the flow diagram from the [filters tutorial](./6-filters.md#using-filters-to-mitigate-network-growth).

Another hook which may be useful is one which stops expansion when a target molecule has been generated. An example is shown below.

```python
network = engine.network_from_file("8-global-hooks")

strat = engine.strat.cartesian(network)

target_hook = engine.hook.target(engine.mol.rdkit("CC(O)O"))
gen_calc = engine.meta.generation("gen")

strat.expand(num_iter=3, global_hooks=[target_hook], reaction_plan=gen_calc)
```

```sh
>>> pprint(
...   [
...     (i, v[0].smiles, v[1])
...     for i, v in
...     enumerate(zip(network.mols, network.mols.meta(keys=["gen"])))
...   ]
... )
[(0, '[H][H]', {'gen': 0}),
 (1, 'O', {'gen': 0}),
 (2, 'CO', {'gen': 0}),
 (3, 'CCO', {'gen': 0}),
 (4, 'CC(=O)O', {'gen': 0}),
 (5, 'CC(O)O', {'gen': 1}),
 (6, 'COC(C)=O', {'gen': 1}),
 (7, 'CCOC(C)=O', {'gen': 1}),
 (8, 'CC(=O)OC(C)=O', {'gen': 1})]
```

Even with all of the same arguments besides the hook, once the target (molecule 5) has been generated, the expansion halts.

## Global Metadata

Global hooks can also be used in order to refresh/update metadata in order to reconcile conflicts that may have accumulated. For example, in the mass waste example above there may be a lower mass waste calculated for a molecule which was previously generated and assigned a mass waste. Therefore, molecule produced from that one must have their metadata updated. A global hook function is useful for such an update.

🚧UNDER CONSTRUCTION🚧

There are not currently any global hook functions providing such a functionality. However, the implementation is relatively simple and looking over the examples in [hooks.py](../pickaxe_generic/hooks.py) should provide some suggestions as to how these may be implemented. In addition, the interface that global hooks subclass from (located in [interfaces.py](../pickaxe_generic/interfaces.py)) contains a docstring for the `__call__` method which describes the meaning of the different `Enum` return types.

## Addendum: Global Hook Functions and the Cartesian Strategy

It may be of some interest that the Cartesian strategy is actually the same as the [priority queue](./9-priority-queue.md) strategy. It accomplishes this by adding a default global hook which limits the number of iterations, as well as defining both the heap size and beam size as uncapped.

This also means that the Cartesian strategy can run out of RAM if the number of possible recipes in a particular generation becomes too large. If a higher-memory solution is desired, a method involving recipe filtering and generation metadata may be used with limited beam/heap sizes in order to guarantee a replication of the Cartesian strategy with a lower memory requirement. However this has not been implemented as of the present time.

## Takeaways

1. Between iterations, a sequence of global hook functions will execute which provide the user a way to intervene in the network expansion in order to reconcile metadata, pause the expansion, limit the expansion based on certain global conditions, or any other functionality which is required.
1. Global hook functions are fairly straightforward to implement, requiring only a single function implementation.
1. The Cartesian strategy is a thin layer over the Priority Queue strategy, and uses its own global function to limit the number of iterations.

Congratulations! You have finished the eighth part of the DORAnet tutorial. Proceed to the [next part](./9-priority-queue.md) to learn how the advanced Priority Queue strategy works and use it to implement a number of network search algorithms.
