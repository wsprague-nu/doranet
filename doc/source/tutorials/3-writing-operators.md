# Writing Operators

In this tutorial, you will learn how to write RDKit SMARTS-based operators. Custom operators are also possible but require more advanced programming knowledge, and come later in the documentation.

Let's create a couple of sample molecules to use in this tutorial.

```python
import doranet as dn

engine = dn.create_engine()

water = engine.mol.rdkit("O")
ethanol = engine.mol.rdkit("CCO")
acetone = engine.mol.rdkit("CC(C)=O")
butanone = engine.mol.rdkit("CCC(C)=O")
methyl_butanoate = engine.mol.rdkit("CCCC(=O)OC")
delta_valerolactone = engine.mol.rdkit("O=C1CCCCO1")
hydroxyvaleric_acid = engine.mol.rdkit("O=C(O)CCCCO")
```

## Reaction Site Matching

Following up on the question from [last time](./2-molecules-and-operators.md), why were there two sets of products?

We can see why using the aldol condensation operator from [before](./2-molecules-and-operators.md#creating-operators).

```python
aldol_condensation = engine.op.rdkit("[O&+0:1]=[C&+0:2]-[C&+0;H2,H3:3].[C&+0:4]=[O&+0:5]>>[*:1]=[*:2]-[*:3]=[*:4].[*:5]")
```

If we perform the aldol condensation on two acetone molecules, we get two product sets which are the same, mesityl oxide and water.

```sh
>>> aldol_condensation(acetone,acetone)
((MolDatBasic('CC(=O)C=C(C)C'), MolDatBasic('O')), (MolDatBasic('CC(=O)C=C(C)C'), MolDatBasic('O')))
```

However, if we choose to use butanone as the first reactant, something changes.

```sh
>>> aldol_condensation(butanone,acetone)
((MolDatBasic('CC(=O)C(C)=C(C)C'), MolDatBasic('O')), (MolDatBasic('CCC(=O)C=C(C)C'), MolDatBasic('O')))
```

Now we have two different products. What's going on?

The answer is because, when matching the first template to butanone, RDKit has a choice of two different carbons for atom 3. In acetone, these carbon atoms are symmetric in the molecule, so the choice doesn't make a difference. In butanone, however, one of them is tertiary and the other is a secondary carbon. Therefore, the product sets are different based on which carbon was picked to match the template.

RDKit expands all permutations of reaction sites. Note also that in the second experiment, acetone was not matched to the first template, because it was the second argument.

Ultimately, if two product sets are the same, they will be considered the same and will not be stored differently. However, if degenerate template matches are of interest to you, there are opportunities to consider this information during a network expansion. This will be discussed later in the metadata section.

## Testing Compatibility

Operators in DORAnet have an additional function, which is used to test the compatibility of a molecule with a particular argument of an operator. This allows the program to cache which molecules can be used with which operators, and with which arguments. Doing this often speeds up network expansion tremendously, so if you are trying to implement a new type of operator, be sure to keep this in mind. Side note: this functionality is also useful for testing out new SMARTS reaction strings.

The function to test argument compatibility is called `.compat()`, and here it is in action.

```sh
>>> aldol_condensation.compat(acetone,0)
True
```

This shows that `acetone` is compatible with the first argument of the `aldol_condensation` operator (numbering starts at 0).

We can also see an example of a molecule which is **not** compatible with `aldol_condensation`.

```sh
>>> aldol_condensation.compat(water,1)
False
```

Water clearly does not possess a carbonyl group, so it fails the check. If an operator is called using an incompatible molecule as an argument, undefined behavior results, and it may even throw an error. In this case water is not matched to a template and no error results, but there are no products.

```sh
>>> aldol_condensation(acetone,water)
()
```

The only guarantee is that if a molecule is compatible via `.compat`, it will produce a correct result unless the operator itself has a bug.

## # of Reactants

If we had an operator, but didn't know the SMARTS used to make it, how would we know how many arguments it has? One thing we could try would be to query the `.smarts` property (only available on RDKit SMARTS objects).

```sh
>>> print(aldol_condensation.smarts)
[O&+0:1]=[C&+0:2]-[C&+0;H2,H3:3].[C&+0:4]=[O&+0:5]>>[*:1]=[*:2]-[*:3]=[*:4].[*:5]
```

From this, we can easily see that there are two reactant templates, separated by `.`. However, a SMARTS reaction string is not always available. The generalized way to test the length of an operator is to call the `len()` function on it.

```sh
>>> print(len(aldol_condensation))
2
```

This shows that the operator has exactly two arguments. If we try to call this operator to react acetone, but expect it to infer that both are acetone, it will raise an exception.

```sh
>>> aldol_condensation(acetone)

Traceback (most recent call last):
  File "C:\Users\sprag\GitHub-Repos\doranet\pickaxe_generic\datatypes.py", line 301, in __call__
    for products in self._rdkitrxn.RunReactants(
ValueError: ChemicalParserException: Number of reactants provided does not match number of reactant templates.

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "<string>", line 1, in <module>
  File "C:\Users\sprag\GitHub-Repos\doranet\pickaxe_generic\datatypes.py", line 306, in __call__
    raise RuntimeError(
RuntimeError: Error occurred when using operator OpDatBasic('[O&+0:1]=[C&+0:2]-[C&+0;H2,H3:3].[C&+0:4]=[O&+0:5]>>[*:1]=[*:2]-[*:3]=[*:4].[*:5]') on (MolDatBasic('CC(C)=O'),)
```

However, you might be realizing a potential problem with this expectation of exact numbers of arguments. What about ring-closing reactions?

## Ring-Forming Reactions in RDKit

Testing one example of esterification, ethanol with 5-hydroxyvaleric acid, shows that the expected products are generated.

```python
esterification = engine.op.rdkit("[O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4]>>[*:1]=[*:2]-[*:4].[*:3]")
```

```sh
>>> esterification(hydroxyvaleric_acid,ethanol)
((MolDatBasic('CCOC(=O)CCCCO'), MolDatBasic('O')),)
```

However, there is an additional reaction which is possible given these reactants. This would be the intramolecular esterification to produce δ-valerolactone. How would such an operation be represented within this system? While you could produce a custom operator type that operates on all subsets of reactants, the easier method is to simply create a new operator which performs the intramolecular variant.

```python
esterification_intra = engine.op.rdkit("([O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4])>>[*:1]=[*:2]-[*:4].[*:3]")
```

Note the parentheses surrounding the reactants. Their presence indicates that even though the `.` character is present, it indicates only a disconnect within the template, but there is ultimately only one reactant template. We can show that the result of this formulation is the intramolecular condensation.

```sh
>>> esterification_intra(hydroxyvaleric_acid)
((MolDatBasic('O=C1CCCCO1'), MolDatBasic('O')),)
```

Trying to perform the intermolecular reaction with the intramolecular operator returns an error as expected.

```sh
>>> esterification_intra(hydroxyvaleric_acid,ethanol)
Traceback (most recent call last):
  File "C:\Users\sprag\GitHub-Repos\doranet\pickaxe_generic\datatypes.py", line 301, in __call__
    for products in self._rdkitrxn.RunReactants(
ValueError: ChemicalParserException: Number of reactants provided does not match number of reactant templates.

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "C:\Users\sprag\GitHub-Repos\doranet\pickaxe_generic\datatypes.py", line 306, in __call__
    raise RuntimeError(
RuntimeError: Error occurred when using operator OpDatBasic('([O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4])>>[*:1]=[*:2]-[*:4].[*:3]') on (MolDatBasic('O=C(O)CCCCO'), MolDatBasic('CCO'))
```

## Ring-Breaking Reactions in RDKit

There is a potentially more dangerous issue when writing Reaction SMARTS, and that is when a bond is broken permanently. The naive implementation below works for breaking non-ring bonds, but struggles when presented with a reactant which undergoes decyclization.

```python
ester_hydrolysis_incorrect = engine.op.rdkit("[O&+0:1]=[C&+0:2]-[O&+0&H0:3].[O&+0&H2:4]>>[*:1]=[*:2]-[*:4].[*:3]")
```

```sh
>>> ester_hydrolysis_incorrect(methyl_butanoate,water)
((MolDatBasic('CCCC(=O)O'), MolDatBasic('CO')),)
>>> ester_hydrolysis_incorrect(delta_valerolactone,water)
((MolDatBasic('CCCCC(=O)O'), MolDatBasic('CCCCO')),)
```

As you can see, the methyl butanoate is properly broken into pieces, but δ-valerolactone is somehow split into two molecules, the sum of which actually have a greater molecular weight than the original!

This is because the products are marked as being separate by the `.` between them. This issue is resolved with some judicious parentheses and use of the `@` ring-bond marker, to split the operator into a ring version and a non-ring version, [as for the bond-forming operator above](./3-writing-operators.md#ring-forming-reactions-in-rdkit).

```python
ester_hydrolysis_nonring = engine.op.rdkit("[O&+0:1]=[C&+0:2]-&!@[O&+0&H0:3].[O&+0&H2:4]>>[*:1]=[*:2]-[*:4].[*:3]")
ester_hydrolysis_ring = engine.op.rdkit("[O&+0:1]=[C&+0:2]-&@[O&+0&H0:3].[O&+0&H2:4]>>([*:1]=[*:2]-[*:4].[*:3])")
```

The results of using these two are shown below.

```sh
>>> ester_hydrolysis_nonring(methyl_butanoate,water)
((MolDatBasic('CCCC(=O)O'), MolDatBasic('CO')),)
>>> ester_hydrolysis_nonring(delta_valerolactone,water)
()
>>> ester_hydrolysis_ring(methyl_butanoate,water)
()
>>> ester_hydrolysis_ring(delta_valerolactone,water)
((MolDatBasic('O=C(O)CCCCO'),),)
```

This demonstrates that, with some additional syntax, the operators will correctly recognize the presence/non-presence of ring bonds and perform their role accordingly.

## Product Specificity and Explicit Hydrogens

One of the big problems you may encounter is if you accidentally overspecify your products. Specifying the number of hydrogens on your product (if you are using an implicit model in RDKit) is a good way to cause it to fail reaction site matching. Try to limit your product specifications to only the connectivity and, if applicable, changes in charge.

At the moment, there are no good examples of the explicit hydrogen problem, but if you are having consistent, inexplicable valence errors deep into your network, this may be why.

## Takeaways

1. RDKit SMARTS-based operators match as many sites as possible, even redundant ones.
2. Extra care must be taken for addition and elimination reactions; always consider the ring case.
3. RDKit SMARTS-based operators should specify as little in the product template as possible.

Congratulations! You have finished the third part of the DORAnet tutorial. Proceed to the [next part](./4-creating-a-basic-network.md) to learn how to store your molecules, operators, and reactions in a network.
