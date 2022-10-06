# Example Operators

This is a list of sample SMARTS operators, which cover some reactions in Organic Chemistry Eighth Edition by Brown, Iverson, Anslyn, and Foote (2017).  Stereoselectivity and chirality are not included, so there is some overlap between operators.

WARNING: these have not necessarily been fully tested, but have been included to provide a reference.

## Reactions of Alkenes

1. Addition of HX (6.3A)
    * [C&+0:1]=[C&+0:2].[F,Cl,Br,I;H1&+0:3]>>[*:1]-[*:2]-[*:3]
    * Overlap with Haloalkanes, Halogenation, and Radical Reactions: HBr Addition to Alkenes Under Radical Conditions
1. Acid-Catalyzed Hydration (6.3B)
    * [C&+0:1]=[C&+0:2].[O&H2&+0:3]>>[*:1]-[*:2]-[*:3]
    * Overlap with Reactions of Alkenes: Oxymercuration-Reduction
    * Overlap with Reactions of Alkenes: Hydroboration-Oxidation
1. Addition of Bromine and Chlorine (6.3D)
    * [C&+0:1]=[C&+0:2].[Cl&+0:3]-[Cl&+0:4]>>[*:3]-[*:1]-[*:2]-[*:4]
    * [C&+0:1]=[C&+0:2].[Br&+0:3]-[Br&+0:4]>>[*:3]-[*:1]-[*:2]-[*:4]
1. Addition of HOCl and HOBr (6.3E)
    * [C&+0:1]=[C&+0:2].[O&H2&+0:3].[Cl&+0:4]-[Cl&+0:5]>>[*:3]-[*:1]-[*:2]-[*:4].[*:5]
    * [C&+0:1]=[C&+0:2].[O&H2&+0:3].[Br&+0:4]-[Br&+0:5]>>[*:3]-[*:1]-[*:2]-[*:4].[*:5]
1. Oxymercuration-Reduction (6.3F)
    * [C&+0:1]=[C&+0:2].[O&H2&+0:3]>>[*:1]-[*:2]-[*:3]
    * Overlap with Reactions of Alkenes: Acid-Catalyzed Hydration
    * Overlap with Reactions of Alkenes: Hydroboration-Oxidation
1. Hydroboration-Oxidation (6.4)
    * Overlap with Reactions of Alkenes: Acid-Catalyzed Hydration
    * Overlap with Reactions of Alkenes: Oxymercuration-Reduction
1. Oxidation to a Vicinal Diol by OsO<sub>4</sub> (6.5A)
    * [C&+0:1]=[C&+0:2].[O&H1&+0:3]-[O&H1&+0:4]>>[*:3]-[*:1]-[*:2]-[*:4]
    * Using method of hydrogen peroxide re-oxidation (Milas hydroxylation)
1. Oxidation by Ozone (6.5B)
    * [C&+0:1]=&!@[C&+0:2].[O&+0:3]=[O&+1:4]-[O&-1:5].[O&+0:6]=[O&+1:7]-[O&-1:8]>>[*:1]=[*:3].[*:2]=[*&+0:4].[*&+0:5]=[*:6].[*&+0:7]=[*&+0:8]
    * [C&+0:1]=&@[C&+0:2].[O&+0:3]=[O&+1:4]-[O&-1:5].[O&+0:6]=[O&+1:7]-[O&-1:8]>>([*:1]=[*:3].[*:2]=[*&+0:4]).[*&+0:5]=[*:6].[*&+0:7]=[*&+0:8]
1. Addition of H<sub>2</sub>; Catalytic Reduction (6.6)
    * [C&+0:1]=[C&+0:2].[#1]-[#1]>>[*:1]-[*:2]
    * Overlap with Reactions of Alkenes: Enantioselective Reduction
1. Enantioselective Reduction (6.7)
    * [C&+0:1]=[C&+0:2].[#1]-[#1]>>[*:1]-[*:2]
    * Overlap with Reactions of Alkenes: Addition of H<sub>2</sub>; Catalytic Reduction

## Alkynes
1. Acidity of Terminal Alkynes (7.4) -> Alkylation of Acetylide Anions (7.5A)
    * [C&+0:1]#[C&H1&+0:2].[C&+0:3]-[F,Cl,Br,I;+0:4]>>[*:1]#[*:2]-[*:3].[*:4]
    * ([C&+0:1]#[C&H1&+0:2].[C&+0:3]-[F,Cl,Br,I;+0:4])>>[*:1]#[*:2]-[*:3].[*:4]
1. Sythesis of an Alkyne from an Alkene (7.5B)
    * [Cl&+0:1]-[C&H1&+0:2]-[C&H1&+0:3]-[Cl&+0:4]>>[*:2]#[*:3].[*:1].[*:4]
    * [Br&+0:1]-[C&H1&+0:2]-[C&H1&+0:3]-[Br&+0:4]>>[*:2]#[*:3].[*:1].[*:4]
    * Technically the second phase, first involves brominating or chlorinating the C=C bond.
1. Addition of Br<sub>2</sub> and Cl<sub>2</sub> (7.6A)
    * [C&+0:1]#[C&+0:2].[Cl&+0:3]-[Cl&+0:4]>>[*:3]-[*:1]=[*:2]-[*:4]
    * [C&+0:1]#[C&+0:2].[Br&+0:3]-[Br&+0:4]>>[*:3]-[*:1]=[*:2]-[*:4]
    * [C&+0:1]#[C&+0:2].[Cl&+0:3]-[Cl&+0:4].[Cl&+0:5]-[Cl&+0:6]>>[*:3]-[*:1]\(-[*:5])-[*:2]\(-[*:6])-[*:4]
    * [C&+0:1]#[C&+0:2].[Br&+0:3]-[Br&+0:4].[Br&+0:5]-[Br&+0:6]>>[*:3]-[*:1]\(-[*:5])-[*:2]\(-[*:6])-[*:4]
1. Addition of HX (7.6B)
    * [C&+0:1]#[C&+0:2].[F,Cl,Br,I;H1&+0:3]>>[*:1]=[*:2]-[*:3]
    * [C&+0:1]#[C&+0:2].[F&H1&+0:3].[F&H1&+0:4]>>[*:1]-[*:2]\(-[*:3])-[*:4]
    * [C&+0:1]#[C&+0:2].[Cl&H1&+0:3].[Cl&H1&+0:4]>>[*:1]-[*:2]\(-[*:3])-[*:4]
    * [C&+0:1]#[C&+0:2].[Br&H1&+0:3].[Br&H1&+0:4]>>[*:1]-[*:2]\(-[*:3])-[*:4]
    * [C&+0:1]#[C&+0:2].[I&H1&+0:3].[I&H1&+0:4]>>[*:1]-[*:2]\(-[*:3])-[*:4]
1. Keto-Enol Tautomerism (7.7A)
    * [C&+0:1]=[C&+0:2]-[O&H1&+0:3]>>[*:1]-[*:2]=[*:3]
    * [C&H1&+0:1]-[C&+0:2]=[O&+0:3]>>[*:1]=[*:2]-[*:3]
1. Hydroboration-Oxidation (7.7A)
    * [C&+0:1]#[C&+0:2].[O&H1&+0:3]-[O&H1&+0:4]>>[*:1]=[*:2]-[*:3].[*:4]
    * [C&+0:1]#[C&+0:2].[O&H1&+0:3]-[O&H1&+0:4]>>[*:1]-[*:2]=[*:3].[*:4]
1. Acid-Catalyzed Hydration (7.7B)
    * [C&+0:1]#[C&+0:2].[O&H2&+0:3]>>[*:1]=[*:2]-[*:3]
    * [C&+0:1]#[C&+0:2].[O&H2&+0:3]>>[*:1]-[*:2]=[*:3]
1. Catalytic Reduction (7.8A)
    * [C&+0:1]#[C&+0:2].[#1][#1].[#1][#1]>>[*:1]-[*:2]
    * [C&+0:1]#[C&+0:2].[#1][#1]>>[*:1]=[*:2]
    * Overlap with Reactions of Alkynes: Hydroboration-Protonolysis
    * Overlap with Reactions of Alkynes: Reduction Using Na or Li Metal in NH<sub>3</sub>(*l*)
1. Hydroboration-Protonolysis (7.8B)
    * [C&+0:1]#[C&+0:2].[#1][#1]>>[*:1]=[*:2]
    * Overlap with Reactions of Alkynes: Catalytic Reduction
    * Overlap with Reactions of Alkynes: Reduction Using Na or Li Metal in NH<sub>3</sub>(*l*)
1. Reduction Using Na or Li Metal in NH<sub>3</sub>(*l*) (7.8C)
    * [C&+0:1]#[C&+0:2].[#1][#1]>>[*:1]=[*:2]
    * Overlap with Reactions of Alkynes: Catalytic Reduction
    * Overlap with Reactions of Alkynes: Hydroboration-Protonolysis

## Haloalkanes, Halogenation, and Radical Reactions

1. Chlorination and Bromination of Alkanes (8.4)
    * [C&!H0&+0:1].[Cl&+0:2]-[Cl&+0:3]>>[*:1]-[*:2].[*:3]
    * [C&!H0&+0:1].[Br&+0:2]-[Br&+0:3]>>[*:1]-[*:2].[*:3]
1. Allylic Bromination (8.6)
    * [C&+0]=[C&+0]-[C&+0:1].[O&+0]=[C&H0&+0]1-[C&H0&+0]-[C&H0&+0]-[C&H0&+0]\(=[O&+0])-[N&+0:2]1-[Br&+0:3]>>[*:1]-[*:3].[*:2]
1. Autoxidation (8.7)
    * [C&+0]=[C&+0]-[C&!H0&+0:1].[O&+0:2]=[O&+0:3]>>[*:1]-[*:2]-[*:3]
1. HBr Addition to Alkenes Under Radical Conditions
    * [C&+0:1]=[C&+0:2].[Br&H1&+0:3]>>[*:1]-[*:2]-[*:3]
    * Overlap with Reactions of Alkenes: Addition of HX

## Nucleophilic Substitution and β-Elimination

1. Nucleophilic Aliphatic Substitution: S<sub>N</sub>2 (9.3)
    * Carbon center limited based on page 388 (also has beta carbon stats)
    * Pick some leaving groups from page 389 for SN2
    * Pick some nucleophiles from page 393 & 376
    * Keep in mind the nature of substitutions on page 394 for halogens
    * Keep in mind effect of solvent from page 392
    * General comparison on Page 397
1. Nucleophilic Aliphatic Substitution: S<sub>N</sub>1 (9.3)
    * Carbon center limited based on page 388
    * Check stability on page 386
    * Keep in mind the nature of substitutions on page 394 for halogens
    * Keep in mind effect of solvent from page 391
    * Keep in mind resonance from page 386
    * General comparison on Page 397