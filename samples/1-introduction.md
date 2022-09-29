# Example 1: Introduction

Example of how the fundamental building blocks of Pickaxe-Generic work.

In this example, we will create some molecules, create a chemical operator, and use the operator to generate some product molecules.  Let's start by importing pickaxe_generic and giving it an alias.

```python
import pickaxe_generic as pg
```

Next, create an Engine; this object will initialize most of your objects and
handle things like database connections, parallel processing, etc. for you.

```python
engine = pg.create_engine()
```

Now we want to create some molecules.  The preferred way for the end user to create molecules is from Daylight SMILES.  Let's try making acetone as an RDKit-style molecule (other styles are not yet supported).

```python
acetone = engine.mol.rdkit("CC(O)C")
print(acetone)
```

Investigating what this molecule consists of shows that the SMILES string appears to be different.  This is because the SMILES string has been "canonicalized" by RDKit.  Because there are many SMILES strings which could represent the same molecule, it becomes difficult to quickly compare two molecules to see if they are the same one.  Therefore, RDKit has a procedure such that for a particular molecule, only one of those SMILES strings is considered valid.  This is the "canonical" SMILES string.

```python
>>> print(acetone)
MolDatBasic('CC(C)O')
```

Let's investigate some further properties of a molecule in the Pickaxe-Generic framework.  For example, all molecules have the `uid` property.  The `uid` property refers to the unique identifying code representing that molecule.  In this case, it is equivalent to the canonical SMILES string.  However, depending on your application this may be something more complex.

```python
>>> print(acetone.uid)
'CC(C)O'
```

If you want to obtain a guaranteed SMILES string, and you are using molecules which support SMILES strings, it can be obtained via the `smiles` property.  Since RDKit-based molecules inherit from `MolDatRDKit`, they support this property.

```python
>>> print(acetone.smiles)
'CC(C)O'
```


