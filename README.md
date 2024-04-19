<div id="top"></div>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->

<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

<!-- PROJECT LOGO -->

<h3 align="center">DORAnet</h3>

  <p align="center">
    This is a repo containing the DORAnet network generation framework.
    <br />
    <a href="https://github.com/wsprague-nu/doranet/issues">Report Bug</a>
    ·
    <a href="https://github.com/wsprague-nu/doranet/issues">Request Feature</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->

## About The Project

This project is a rework of both [minedatabase](https://pypi.org/project/minedatabase/) and [NetGen](https://github.com/BroadbeltLab/NetGen). It is intended to streamline and extend the implementation of chemical reaction network generation for research purposes.
The aforementioned programs possess some weaknesses in terms of a lack of transferability, customizability, and ease of use. DORAnet is intended to solve these problems using the industry-standard molecule manipulation framework RDKit (though others can be implemented), an object-oriented approach to network generation, and a flexible dependency injection scheme which can handle low-level performance tradeoffs such as speed-memory and parallelization.

<p align="right">(<a href="#top">back to top</a>)</p>

### Built With

- [Anaconda](https://www.anaconda.com/)
- [Git](https://git-scm.com/)
- [GitHub](https://github.com/)
- [Python](https://www.python.org/)
- [RDKit](https://rdkit.org/)
- [Visual Studio Code](https://code.visualstudio.com/)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- GETTING STARTED -->

## Getting Started

To get a local copy up and running follow these simple example steps.

### Prerequisites

You will need access to the Python package rdkit. The recommended way to install rdkit is via the Anaconda package/environment manager (Nix flake coming soon!). If you do not have Anaconda on your computer, install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains a minimal Anaconda package manager.

### Installation

1. Activate Anaconda in a terminal (CMD/bash/etc.) using the methods provided in their documentation.
2. Navigate to the folder where you want to download pickaxe-generic using cd. Example below.
   ```sh
   cd $USERPROFILE/GitHub-Repos
   ```
3. Clone the repo using the command below. If using Windows, use Git Bash for this step only (an alternate terminal which comes with [Git](https://git-scm.com/)). Enter your GitHub credentials for the account which has been granted permissions to access the repo. You may need a [Personal Access Token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token#using-a-token-on-the-command-line).
   ```sh
   git clone https://github.com/wsprague-nu/doranet.git
   ```
4. Enter the doranet/ folder in the Anaconda terminal.
   ```sh
   cd doranet
   ```
5. Create the rdkit/pickaxe environment using the guidelines below:

   - If developing new code, environment-dev.yml is recommended.
   - If using MINE-Database version of Pickaxe, environment-compat.yml is recommended.
   - If using only DORAnet, environment.yml is recommended.

   Command to create new environment (replace env.yml with your chosen environment file)

   ```sh
   conda env create -f env.yml
   ```

   Existing environment (replace env.yml with your chosen environment file, and env-name with the existing environment name)

   ```sh
   conda env update -n env-name -f env.yml --prune
   ```

6. Activate the environment you installed rdkit into (replace env-name with your chosen environment name). If a fresh environment was installed, the default environment name should be displayed on the screen.
   ```sh
   conda activate env-name
   ```
7. Install doranet using pip. Use a -e flag after "install" if you want your installation to update automatically when changing the files in this folder (recommended). Otherwise, simply use the command below.
   ```sh
   python -m pip install .
   ```
8. When running a Python program requiring doranet, make sure you first open your terminal and activate the relevant environment using the command from Step 6.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->

## Usage

As DORAnet is intended as to be an extensible, polymorphic network generation software, effective users should understand the basic architecture of the system. While each class has its own documentation, a bird's eye view of how the program is organized should assist in development of new code and strategies with minimal overhead.

There are three design properties which informed the abstraction of network generation into an object oriented architecture:

- Compatibility
- Efficiency
- Extensibility

-future elaboration on these properties here-

The network has its start in the Engine object, provided by create_engine. This object is provided the relevant configuration options for the network expansion, such as the number of available cores (parallelism not yet available), speed/memory tradeoffs required, type of expansion strategy etc., and in exchange provides the relevant objects which meet those criteria. The Engine is the only object which requires knowledge of the entire class hierarchy since it provides them through several straightforward interfaces.

The ObjectLibrary is in charge of data storage. It provides a generalized interface for storing generic DataUnit objects such as molecules, operators, and reactions. The ObjectLibrary is abstracted since it is future-compatible with disk and external database storage. It works as an iterable, allowing the user to iterate through whichever components it contains, and also allows for lookup much like a dict using the .uid property of the stored DataUnit objects. Much like a database, a network consists of three ObjectLibraries holding molecules, operators, and reactions respectively.

The DataUnit is the generic abstract class defining an atomic unit of data. Implementing classes provide the UID, a unique identifier used as a key in ObjectLibrary objects, as well as a `.blob` function to translate the data contained within into a compressed bytestring.

- MolDat (or Molecule Data) objects represent molecules. These are currently only RDKit molecules, but could feasibly represent anything. The RDKit subclassed versions provide, in addition to the DataUnit properties, the .rdkitmol property to access the RDKit molecule and the .smiles property to access the RDKit canonical SMILES string.
- OpDat (or Operator Data) objects represent operators. These are currently only RDKit SMARTS operators, but could feasibly represent almost anything. They provide the len() and compat(MolDat,int) methods, which give the number of arguments to the operator and the compatibility of a particular MolDat with a particular argument, used to speed up reaction generation. They are also callable, and return an iterable of possible reaction product sets.
- RxnDat (or Reaction Data) objects represent reactions. These are, at the moment, not much more than a combination of a set containing reactant uids accessible through .reactants, a set containing product uids accessible through .products, and the operator uid accessible through .operator.

The Strategy puts all these components together to generate a network. The only Strategy which has a full implementation is the CartesianStrategy. This strategy attempts to combine every operator with every combination of compatible molecules to expand the network. A number of "generations" can be specified, which represent the number of times the Cartesian product is performed, with the network expanding every time. A reaction-level filter, implemented by the user, can filter out new reactions based on particular criteria in order to restrict the growth of the network. A holistic filter, which filters out molecules based on an entire new generation, is recommended to be implemented separately by the end user, but this may change.

Work in progress. Check out example_notebook.ipynb in jupyter for examples. Be sure to first install Jupyter using "conda install jupyter" while your environment is activated.

<!--
This is an example of how pickaxe-generic may be used to obtain the heat of formation of an arbitrary molecule (for which the Benson groups exist in primary_groups).

   ```python
   import ngthermo.properties as prop

   smiles = 'CC1CC(=O)CC(=O)O1'
   Hf = prop.Hf(smiles) / 1000 # Hf provided in cal/mol
   print(f'Enthalpy of {smiles}: {Hf} kcal/mol)
   ```
-->

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CONTACT -->

## Contact

Your Name - wsprague@u.northwestern.edu

Project Link: [https://github.com/wsprague-nu/doranet](https://github.com/wsprague-nu/doranet)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->

## Acknowledgments

- [Prof. Linda Broadbelt](broadbelt.northwestern.edu)
- Manali Dhawan
- Quan Zhang
- Kevin Shebek
- Lauren Lopez
- Joseph Ni
- Jon Strutz
- [Prof. Brent Shanks](https://www.engineering.iastate.edu/people/profile/bshanks/)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->

[contributors-shield]: https://img.shields.io/github/contributors/wsprague-nu/doranet.svg?style=for-the-badge
[contributors-url]: https://github.com/wsprague-nu/doranet/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/wsprague-nu/doranet.svg?style=for-the-badge
[forks-url]: https://github.com/wsprague-nu/doranet/network/members
[stars-shield]: https://img.shields.io/github/stars/wsprague-nu/doranet.svg?style=for-the-badge
[stars-url]: https://github.com/wsprague-nu/doranet/stargazers
[issues-shield]: https://img.shields.io/github/issues/wsprague-nu/doranet.svg?style=for-the-badge
[issues-url]: https://github.com/wsprague-nu/doranet/issues
[license-shield]: https://img.shields.io/github/license/wsprague-nu/doranet.svg?style=for-the-badge
[license-url]: https://github.com/wsprague-nu/doranet/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[product-screenshot]: images/screenshot.png
