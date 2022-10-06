# Example Operators

This is a list of sample SMARTS operators, which cover some reactions in Organic Chemistry Eighth Edition by Brown, Iverson, Anslyn, and Foote (2017).  Stereoselectivity and chirality are not included, so there is some overlap between operators.

WARNING: these have not necessarily been fully tested, but have been included to provide a reference.

## Reactions of Alkenes (Chapter 6)

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

## Alkynes (Chapter 7)
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

## Haloalkanes, Halogenation, and Radical Reactions (Chapter 8)

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

## Nucleophilic Substitution and β-Elimination (Chapter 9)

1. Nucleophilic Aliphatic Substitution: S<sub>N</sub>2 (9.3)
    * Carbon center limited based on page 388 (also has beta carbon stats)
    * Pick some leaving groups from page 389 for SN2
    * Pick some nucleophiles from page 393 & 376
    * Keep in mind the nature of substitutions on page 394 for halogens
    * Keep in mind effect of solvent from page 392
    * General comparison on Page 397, 413
1. Nucleophilic Aliphatic Substitution: S<sub>N</sub>1 (9.3)
    * Carbon center limited based on page 388
    * Check stability on page 386
    * Keep in mind the nature of substitutions on page 394 for halogens
    * Keep in mind effect of solvent from page 391
    * Keep in mind resonance from page 386
    * General comparison on Page 397, 413
1. β-Elimination: E1 (9.6, 9.7)
1. β-Elimination: E2 (9.6, 9.7)
1. Neighboring Group Participation

## Alcohols (Chapter 10)

1. Acidity of Alcohols (10.3)
1. Reaction with Active Metals (10.4)
1. Reaction with HCl, HBr, and HI (10.5A)
1. Reaction with PBr<sub>3</sub>
1. Reaction with SOCl<sub>2</sub> and SOBr<sub>2</sub> (10.5C)
1. Acid-Catalyzed Dehydration (10.6)
1. Pinacol Rearrangement (10.7)
1. Oxidation of a Primary Alcohol to a Carboxylic Acid (10.8A)
1. Oxidation of a Secondary Alcohol to a Ketone (10.8A-10.8D)
1. Oxidation of a Primary Alcohol to an Aldehyde (10.8B)
1. Oxidative Cleavage of a Glycol (10.8E)
1. Acidity of Thiols (10.9F)
1. Oxidation of Thiols to DIsulfides (10.9G)

## Ethers, Epoxides, and Sulfides (Chapter 11)

1. Williamson Ether Synthesis (11.4A)
1. Acid-Catalyzed Dehydration of Alcohols (11.4B)
1. Acid-Catalyzed Addition of Alcohols to Alkenes (11.4C)
1. Acid-Catalyzed Cleavage of Dialkyl Ethers (11.5A)
1. Reaction of Alcohols with Chloro-*tert*-butyldimethylsilane (11.6)
1. Oxidation of Alkenes with Peroxycarboxylic Acids (11.8C)
1. Synthesis of Epoxides from Halohydrins (11.8C)
1. Sharpless Asymmetric Epoxidation (11.8D)
1. Acid-Catalyzed Hydrolysis of Epoxides (11.9A)
1. Nucleophilic Ring Opening of Epoxides (11.9B)
1. Reduction of an Epoxide to an Alcohol (11.9B)
1. Oxidation of Sulfides (11.12C)

## An Introduction to Organometallic Compounds (Chapter 15)

1. Formation of Organomagnesium (Grignard) and Organolithium Compounds (15.1A)
1. Reaction of RMgX and RLi with Proton Donors (15.1B)
1. Reaction of a Grignard Reagent with an Epoxide (15.1C)
1. Formation of Gilman Reagents (15.2A)
1. Treatement of a Gilman Reagent with an Alkyl, Aryl, or Alkenyl Halide (15.2B)
1. Reaction of Dichloro- or Dibromocarbene with an Alkene (15.3B)
1. The Simmons-Smith Reaction (15.3C)

## Aldehydes and Ketones (Chapter 16)

1. Reaction with Grignard Reagents (16.5A)
1. Reaction with Organolithium Reagents (16.5B)
1. Reaction with Anions of Terminal Alkynes (16.5C)
1. Reaction with HCN to Form Cyanohydrins (16.5D)
1. The Wittig Reaction (16.6)
1. The Horner-Emmons-Wadsworth Modification of the Wittig Reaction (16.6)
1. Hydration (16.7A)
1. Addition of Alcohols to Form Hemiacetals (16.7B)
1. Addition of Alcohols to Form Acetals (16.7B)
1. Addition of Ammonia and Its Derivatives: Formation of Imines (16.8A)
1. Addition of Secondary Amines: Formation of Enamines (16.8A)
1. Addition of Hydrazine and Its Derivatives (16.8B)
1. Keto-Enol Tautomerism (16.9B)
1. Oxidation of an Aldehyde to a Carboxylic Acid (16.10A)
1. Metal Hydride Reduction (16.11A)
1. Catalytic Reduction (16.11B)
1. Reductive Amination (16.11D)
1. Clemmensen Reduction (16.11E)
1. Wolff-Kishner Reduction (16.11E)
1. Deuterium Exchange at an α-Carbon (16.12B)
1. Halogenation at an α-Carbon (16.12C)

## Carboxylic Acids (Chapter 17)

1. Acidity of Carboxylic Acids (17.4A)
1. Reaction of Carboxylic Acids with Bases (17.4B)
1. Carbonation of a Grignard Reagent (17.5)
1. Industrial Preparation of Acetic Acid by the Carbonylation of Methanol (17.6)
1. Reduction by Lithium Aluminum Hydride (17.6A)
1. Fischer Esterification (17.7A)
1. Reaction with Diazomethane (17.7B)
1. Conversion to Acid Halides (17.8)
1. Decarboxylation of β-Ketoacids (17.9A)
1. Decarboxylation of β-Dicarboxylic Acids (17.9B)

## Functional Derivatives of Carboxylic Acids (Chapter 18)

1. Acidity of Imides (18.2)
1. Acidity of Sulfonamides (18.2)
1. Hydrolysis of an Acid Chloride (18.4A)
1. Hydrolysis of an Acid Anhydride (18.4B)
1. Hydrolysis of an Ester (18.4C)
1. Hydrolysis of an Amide (18.4D)
1. Hydrolysis of a Nitrile (18.4E)
1. Reaction of an Acid Chloride with an Alcohol (18.5A)
1. Reaction of an Acid Anhydride with an Alcohol (18.5B)
1. Reaction of an Ester with an Alcohol: Transesterification (18.5C)
1. Reaction of an Acid Chloride with Ammonia or an Amine (18.6A)
1. Reaction of an Acid Anhydride with Ammonia or an Amine (18.6B)
1. Reaction of an Ester with Ammonia or an Amine (18.6C)
1. Reaction of an Acid Chloride with a Carboxylic Acid Salt (18.7)
1. Reaction of an Ester with a Grignard Reagent (18.9A)
1. Reaction of an Acid Chloride with a Lithium Diorganocuprate (18.9C)
1. Reduction of an Ester (18.10A)
1. Reduction of an Amide (18.10B)
1. Reduction of a Nitrile (18.10C)

## Enolate Anions and Enamines (Chapter 19)

1. Aldol Reaction (19.2)
1. Dehydration of the Product of an Aldol Reaction (19.2)
1. Claisen Condensation (19.3A)
1. Dieckmann Condensation (19.3B)
1. Alkylation of an Enamine Followed by Hydrolysis (19.5A)
1. Acylation of an Enamine Followed by Hydrolysis (19.5B)
1. Acetoacetic Ester Synthesis (19.6)
1. Malonic Ester Synthesis (19.7)
1. Michael Reaction (19.8A)
1. Robinson Annulation (19.8C)
1. Conjugate Addition of Lithium Diorganocopper Reagents (19.8E)
