# Molecules and Operators

In this tutorial, you will learn how to create a molecule, create a chemical operator, and use the operator to perform a reaction. Let's start by importing doranet and giving it an alias.

## Importing/configuring

```python
import doranet as dn
```

Next, create an Engine; this object will initialize most of your objects and
handle things like database connections, parallel processing, etc. for you.

```python
engine = dn.create_engine()
```

## Creating molecules

Now we want to create some molecules. The preferred way for the end user to create molecules is from [Daylight SMILES](https://daylight.com/dayhtml/doc/theory/theory.smiles.html). Let's try making acetone as an RDKit-style molecule (other styles are not yet supported).

```python
acetone = engine.mol.rdkit("CC(=O)C")
```

Investigating what this molecule consists of shows that the SMILES string appears to be different. This is because the SMILES string has been "canonicalized" by RDKit. Because there are many SMILES strings which could represent the same molecule, it becomes difficult to quickly compare two molecules to see if they are the same one. Therefore, RDKit has a procedure such that for a particular molecule, only one of those SMILES strings is considered valid. This is the "canonical" SMILES string.

```sh
>>> acetone
MolDatBasic('CC(=O)C')
```

Let's investigate some further properties of a molecule in the DORAnet framework. For example, all molecules have the `uid` property. The `uid` property refers to the unique identifying code representing that molecule. In this case, it is equivalent to the canonical SMILES string. However, depending on your application this may be something more complex.

```sh
>>> acetone.uid
'CC(C)=O'
```

If you want to obtain a guaranteed SMILES string, and you are using molecules which support SMILES strings, it can be obtained via the `smiles` property. Since RDKit-based molecules inherit from `MolDatRDKit`, they support this property.

```sh
>>> acetone.smiles
'CC(C)=O'
```

## Creating operators

Now let's try creating an operator and performing a reaction.

The example reaction will be that of an aldol condensation. When applied to two acetone molecules, it will produce mesityl oxide. We will use an RDKit-style "templated operator" which is defined by [Reaction SMARTS](https://www.daylight.com/dayhtml_tutorials/languages/smarts/#RXN) with some [RDKit-specific extensions](https://www.rdkit.org/docs/RDKit_Book.html#smarts-support-and-extensions).

Just like molecule objects, operator objects are obtained via the engine.

```python
aldol_condensation = engine.op.rdkit("[O&+0:1]=[C&+0:2]-[C&+0;H2,H3:3].[C&+0:4]=[O&+0:5]>>[*:1]=[*:2]-[*:3]=[*:4].[*:5]")
```

We can run this operator on two acetone molecules to condense them, creating new molecules.

```sh
>>> aldol_condensation(acetone,acetone)
((MolDatBasic('CC(=O)C=C(C)C'), MolDatBasic('O')),
 (MolDatBasic('CC(=O)C=C(C)C'), MolDatBasic('O')))
```

You may notice that there are two sets of products, and that these are equivalent. This is because there are technically two choices for reactant atom 3 on acetone; this will be further explored in the next part.

## Takeaways

1. The Engine contains configuration information and can initialize most objects in DORAnet.
2. Molecules obtained via `engine.mol.rdkit` wrap an RDKit-style molecule using Daylight SMILES syntax.
3. Operators obtained via `engine.op.rdkit` wrap an RDKit-style template-based operator using Daylight Reaction SMARTS syntax.

Congratulations! You have finished the second part of the DORAnet tutorial. Proceed to the [next part](./3-writing-operators.md) to learn some nuances of writing operators in the RDKit style.
