# DORAnet Quick Start: Hybrid Pathways

In this quick start, you’ll install DORAnet, generate an enzymatic (forward) and a synthetic (retro) network, and run post-processing to find and visualize paths from your starters to a target.

## Install

Install DORAnet from PyPI:

```bash
pip install doranet
```

> For the best pathway PDF layout later, install Graphviz + PyGraphviz (optional but recommended):
>
> ```bash
> conda install conda-forge::pygraphviz
> ```

## Imports

```python
import doranet.modules.enzymatic as enzymatic
import doranet.modules.synthetic as synthetic
import doranet.modules.post_processing as post_processing
```

* The `enzymatic` module wraps assumptions for enzymatic network expansion.
* The `synthetic` module is for chemical expansion.
* The `post_processing` module collects functions that find, rank, and visualize pathways.

## Define Inputs (SMILES)

We’ll choose the **starters**, **helpers**, and **target** as SMILES strings. SMILES for common molecules can be found on [PubChem](https://pubchem.ncbi.nlm.nih.gov/) and [Wikipedia](https://www.wikipedia.org/).

You can also provide filenames (e.g., `.txt`, `.csv`) where each line is a SMILES. For example: `starters = "starters.txt"`.

> **Helpers (used in chemical expansions only).** Helpers are molecules (e.g., H₂, H₂O, O₂) that can react with starters but cannot be the only reactants in a reaction. They are optional, but many chemical rules require common helpers such as oxygen and water. In the final pathway PDF, helpers are listed next to reactions (not rendered as molecule images).

```python
user_starters = {"CCO"}                    # "CCO" for ethanol
user_helpers = {"O", "O=O", "[H][H]"}      # "O", "O=O", "[H][H]" for water (O), oxygen (O=O), and hydrogen([H][H])
user_target  = {"CC(O)=O"}                 # "CC(O)=O" for acetic acid
job_name = "Acetic_acid_hybrid"
```

## Generate Networks

### 1) Enzymatic (Forward) Expansion

Generate a one-generation forward network using built-in enzymatic operators (reaction rules). The operators are stored in a TSV under `doranet/modules/enzymatic`, and entries include UniProt IDs for the source reactions.

```python
forward_network = enzymatic.generate_network(
    job_name=job_name,
    starters=user_starters,
    gen=1,
    direction="forward",
)
```

### 2) Synthetic (Retro) Expansion

Generate a one-generation retrosynthetic network using chemical rules. Here we set the “starter” set to the target molecule.

```python
retro_network = synthetic.generate_network(
    job_name=job_name,
    starters=user_target,
    helpers=user_helpers,
    gen=1,
    direction="retro",
)
```

## Post-Processing (One Step)

`post_processing.one_step` searches for connections between your starters and target across the two networks, then ranks and visualizes pathways. It will write several outputs to disk, including a PDF.

If Graphviz/PyGraphviz is **not** installed, DORAnet will use a custom layout, which may be less readable for complex pathways.

You can also use pre-existing network files on disk by passing their names (e.g., `networks = {"Acetic_acid_hybrid__forward_saved_network"}`).

```python
if __name__ == "__main__":
    post_processing.one_step(
        networks={
            forward_network,
            retro_network,
        },
        total_generations=2,
        starters=user_starters,
        helpers=user_helpers,
        target=user_target,
        job_name=job_name,
    )
```

---

# Extra 1 — Post-Processing in Steps

`one_step` internally runs four stages:

1. **pretreat\_networks** – merges networks, sanitizes reactions, and saves a JSON.
2. **pathway\_finder** – searches for pathways; saves a TXT of all pathways and files for Reaxys queries.
3. **pathway\_ranking** – ranks pathways; saves a ranked TXT.
4. **pathway\_visualization** – saves a PDF of all pathways.

If you want to tweak ranking weights (or insert a Reaxys step) without redoing pretreatment/search, run the stages manually:

```python
# 1) Pretreatment
post_processing.pretreat_networks(
    networks={
        forward_network,
        retro_network,
    },
    total_generations=2,
    starters=user_starters,
    helpers=user_helpers,
    job_name=job_name,
)

# 2) Search
post_processing.pathway_finder(
    starters=user_starters,
    helpers=user_helpers,
    target=user_target,
    search_depth=2,
    max_num_rxns=2,
    min_rxn_atom_economy=0.5,
    job_name=job_name,
)

# 3) Ranking
post_processing.pathway_ranking(
    starters=user_starters,
    helpers=user_helpers,
    target=user_target,
    num_process=2,
    job_name=job_name,
)

# 4) Visualization
post_processing.pathway_visualization(
    starters=user_starters,
    helpers=user_helpers,
    num_process=2,
    job_name=job_name,
)
```

---

# Extra 2 — Thermodynamic Filters

You can filter reactions during expansion by providing thermodynamic calculators and a maximum allowed change. DORAnet doesn’t ship calculators, but you can pass your own.

* **Chemical expansion**: Provide a **molecule** property calculator (e.g., enthalpy of formation for a SMILES).
* **Enzymatic/bio expansion**: Provide a **reaction** property calculator that consumes a dict with `"reactants"` and `"products"` and returns the thermodynamic change (e.g., ΔG).

```python
# Chemical expansion
def mol_dH(smiles: str) -> float:
    # Compute and return a thermodynamic value for this molecule.
    return 0.0

retro_network = synthetic.generate_network(
    job_name=job_name,
    starters=user_target,
    helpers=user_helpers,
    gen=1,
    direction="retro",
    molecule_thermo_calculator=mol_dH,
    max_rxn_thermo_change=15,
)

# Bio/enzymatic expansion
def rxn_dG(rxn_dict: dict) -> float:
    reactants = rxn_dict["reactants"]
    products  = rxn_dict["products"]
    # Compute and return the thermodynamic change for this reaction.
    return 0.0

retro_network = enzymatic.generate_network(
    job_name=job_name,
    starters={"OC(=O)C(=O)CCCO"},
    gen=1,
    direction="retro",
    rxn_thermo_calculator=rxn_dG,
    max_rxn_thermo_change=0,
)
```

---

# Extra 3 — Optional Arguments and Their Default Values (Reference)

**`synthetic.generate_network`**

```python
generate_network(
    job_name="default_job",
    starters=False,
    helpers=False,
    gen=1,
    direction="forward",
    molecule_thermo_calculator=None,
    max_rxn_thermo_change=15,
    max_atoms=None,                      # e.g., max_atoms={"C": 20, "O": 10} caps per-element atom counts in any product molecule.
    allow_multiple_reactants="default",  # Default True in forward, False in retro. If False, a reactant may react with itself or helpers, but not with other reactants.
    targets=None,                        # String/list/set... Check presence after expansion.
)
```

**`enzymatic.generate_network`**
*No user helpers in enzymatic expansion. DORAnet ships a cofactor list used by bio rules (similar role to helpers).*

```python
generate_network(
    job_name="default_job",
    starters=False,
    gen=1,
    direction="forward",
    rxn_thermo_calculator=None,
    max_rxn_thermo_change=15,
    max_atoms=None,
    allow_multiple_reactants=False,
    targets=None,
)
```

**`post_processing`**

```python
pretreat_networks(
    networks=None,                     # Can be more than 2 networks. e.g., {forward_1, forward_2, forward_3, retro}.
    total_generations=1,               # e.g., 2-gen forward + 3-gen retro expansion, set total_generations=5.
    starters=None,
    helpers=None,
    job_name="default_job_name",
    remove_pure_helpers_rxns=False,    # If True, remove reactions whose reactants are only helpers.
    sanitize=True,                     # If True, remove molecules unreachable from starters/helpers within total_generations.
    transform_enols_flag=False,        # If True, convert enol products to keto forms.
    molecule_thermo_calculator=None,   # Used to assess thermo change during enol transformation.
)

pathway_finder(
    starters=None,
    helpers=None,
    target=None,
    search_depth=1,                 # Should not exceed total_generations from pretreat.
    max_num_rxns=1,                 # Max reactions per pathway.
    min_rxn_atom_economy=0.3,       # 0–1.
    job_name="default_job_name",
    consider_name_difference=True,  # If True, same reactions with different names are treated as distinct.
)

pathway_ranking(
    starters=None,
    helpers=None,
    target=None,
    weights=None,                     # Defaults:
                                      # {"reaction_thermo": 2,
                                      #  "number_of_steps": 4,
                                      #  "by_product_number": 2,
                                      #  "atom_economy": 1,
                                      #  "salt_score": 0,
                                      #  "in_reaxys": 0,
                                      #  "coolness": 0}
    num_process=1,                    # Number of CPU processes
    reaxys_result_name=None,          # CSV filename.
    job_name="default_job_name",
    cool_reactions=None,
    molecule_thermo_calculator=None,  # For by-product scoring.
    max_rxn_thermo_change=15,
)

pathway_visualization(
    starters=None,
    helpers=None,
    num_process=1,
    reaxys_result_name="default",
    job_name="default_job_name",
    exclude_smiles=None,         # Set/list... Exclude pathways containing these molecules.
    reaxys_rxn_color="blue",
    normal_rxn_color="black",
)
```

---

# Extra 4 — Reaxys Batch Query

Reference (batch query submission): [https://service.elsevier.com/app/answers/detail/a\_id/26151/supporthub/reaxys/p/10958/](https://service.elsevier.com/app/answers/detail/a_id/26151/supporthub/reaxys/p/10958/)

`pathway_finder` writes:

* a **TXT** file listing all reactions across all pathways (suitable as a Reaxys batch query), and
* a **CSV** with placeholder zeros for Reaxys results.

If you have Reaxys access, upload the batch query and copy the result log into the CSV. Those results can then inform ranking and visualization.

---

# Extra 5 — Install Environment

It’s recommended to use a virtual environment (Conda is popular).

1. **Install Miniconda or Anaconda**

   * Miniconda: [https://docs.anaconda.com/free/miniconda/](https://docs.anaconda.com/free/miniconda/)
   * Anaconda: [https://docs.anaconda.com/free/anaconda/install/](https://docs.anaconda.com/free/anaconda/install/)

2. **Run the installer**
   Follow the on-screen instructions.

3. **Open a terminal**

   * Linux & macOS: built-in Terminal
   * Windows: **Anaconda PowerShell Prompt**

4. **Create and activate an environment**

   * Use an existing env:

     ```bash
     conda activate your_env_name
     ```
   * Or create a new env:

     ```bash
     conda create -n your_env_name python=3.10
     conda activate your_env_name
     ```

5. **Install DORAnet**

   ```bash
   pip install doranet
   ```

---

# Extra 6 — Common Issues

* **Windows + multiprocessing:** Wrap post-processing calls in a main guard to avoid multiprocessing issues.

```python
if __name__ == "__main__":
    post_processing.one_step(...)
```

---

## What’s Next?

* The functions introduced here provide a simplified interface to DORAnet’s core functions. See the other tutorial pages for in-depth operations.

