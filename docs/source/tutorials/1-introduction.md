# Introduction

In this tutorial, you will learn what background information is necessary to use DORAnet, what purpose it is meant to serve, and the basic ideas behind it.

## Goals and Background

DORAnet is a generic reaction network generation framework, capable of putting together a number of heterogeneous molecule and operator types to explore chemical space. It also has the flexibility to store its data however needed and to expand the network according to numerous different strategies. While the initial codebase contains only a few implementations, these may be easily extended using the interfaces provided.

The creation of this package was motivated by a need to enable groups with different expertise to work on the same synthetic network generation program, with minimal computational experience or end-user overhead. It is also intended to permit rapid prototyping of different network expansion schemes and operator types without requiring extensive overhauls to core routines.

Despite these goals, there is still some background which is absolutely necessary to use this software.

1. Must be able to use [Python](https://www.python.org/) at a basic level.
2. Must understand how to interpret and write [Daylight SMILES](https://daylight.com/dayhtml/doc/theory/theory.smiles.html) and [Daylight SMARTS](https://www.daylight.com/dayhtml_tutorials/languages/smarts).
3. Must have sufficient chemistry knowledge to interpret results.

## Core Concepts

There are really only two primary building blocks in DORAnet: **molecules** and **operators**.

A **molecule** represents little other than a container of basic information. One example of such information is a SMILES string, which represents a molecular structure. However, a molecule could also be an association of tautomers, resonance structures, entire lumped molecular classes, or a statistical polymer. The set of all molecules is termed **chemical space**.

An **operator** represents a function which takes molecules as arguments and returns sets of product molecules. The reason there can be multiple sets of product molecules is due to the possibility for different reaction sites being activated on the same reagents, thereby producing different products. These are separated into groups so that quantities like enthalpy of reaction for each product set may be calculated.

A **reaction** represents the association of an ordered set of reactants, an ordered set of products, and an operator. It does not have any inherent/defining features beyond this, though calculated quantities like enthalpy of reaction may be associated with a reaction.

A graph which contains relationships between **molecules** and **operators** via associative **reaction** nodes is a **chemical network**.

Once you have defined your baseline of what type of initial molecules and operators to use, and specify the characteristics of each one, you can start generating a chemical network right away. This will be covered in the next tutorial.

## Takeaways

1. DORAnet is based on **molecules** and **operators**.
2. **Molecules** are elements of **chemical space**.
3. **Operators** are functions which map an ordered set of molecules (as reactants) to a set of ordered sets of molecules (as products).
4. **Reactions** are a purely associative construct which represent a particular ordered set of reactant molecules, an operator, and an ordered set of product molecules, where the set of products is one of the sets mapped to by the operator.

Congratulations! You have finished the first part of the DORAnet tutorial. If you have installed DORAnet already, proceed to the [next part](./2-molecules-and-operators.md) to start doing some code. If not, reference README.md in the top level of this project.
