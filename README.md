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

<h3 align="center">Pickaxe-Generic</h3>

  <p align="center">
    This is a private repo providing the Pickaxe-Generic network generation framework to the Broadbelt group at Northwestern.  Licensing has not yet been sorted out and this code should not be shared with external developers.
    <br />
    <a href="https://github.com/wsprague-nu/Pickaxe-Generic/issues">Report Bug</a>
    ·
    <a href="https://github.com/wsprague-nu/Pickaxe-Generic/issues">Request Feature</a>
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
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

This project is a rework of both [minedatabase](https://pypi.org/project/minedatabase/) and [NetGen](https://github.com/BroadbeltLab/NetGen).  It is intended to streamline and extend the implementation of chemical reaction network generation for research purposes.
The aforementioned programs possess some weaknesses in terms of a lack of transferability, customizability, and ease of use.  Pickaxe-Generic is intended to solve these problems using the industry-standard molecule manipulation framework RDKit (though others can be implemented), an object-oriented approach to network generation, and a flexible dependency injection scheme which can handle low-level performance tradeoffs such as speed-memory and parallelization.

<p align="right">(<a href="#top">back to top</a>)</p>



### Built With

* [Anaconda](https://www.anaconda.com/)
* [Git](https://git-scm.com/)
* [GitHub](https://github.com/)
* [Python](https://www.python.org/)
* [RDKit](https://rdkit.org/)
* [Visual Studio Code](https://code.visualstudio.com/)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple example steps.

### Prerequisites

You will need access to the Python package rdkit.  The recommended way to install rdkit is via the Anaconda package/environment manager.  If you do not have Anaconda on your computer, install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains a minimal Anaconda package manager.

Create the rdkit/pickaxe environment using the guidelines below:

* If developing new code, environment-dev.yml is recommended.
* If using MINE-Database version of Pickaxe, environment-compat.yml is recommended.
* If using only pickaxe-generic, environment.yml is recommended.

Command to create new environment (replace env.yml with your chosen environment file)
```sh
conda create -f env-name.yml -c rdkit rdkit
```
Existing environment (replace env-name with your chosen environment name)
```sh
conda install -n env-name -c rdkit rdkit
```

### Installation

1. Activate Anaconda in a terminal (CMD/bash/etc.) using the methods provided in their documentation.
2. Navigate to the folder where you want to download pickaxe-generic using cd.  Example below.
   ```sh
   cd $USERPROFILE/GitHub-Repos
   ```
3. Clone the repo using the command below.  If using Windows, use Git Bash for this step only (an alternate terminal which comes with [Git](https://git-scm.com/)).  Enter your GitHub credentials for the account which has been granted permissions to access the repo.  You may need a [Personal Access Token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token#using-a-token-on-the-command-line).
   ```sh
   git clone https://github.com/wsprague-nu/Pickaxe-Generic.git
   ```
4. Enter the Pickaxe-Generic/ folder in the Anaconda terminal.
   ```sh
   cd Pickaxe-Generic
   ```
5. Create the rdkit/pickaxe environment using the guidelines below:

   * If developing new code, environment-dev.yml is recommended.
   * If using MINE-Database version of Pickaxe, environment-compat.yml is recommended.
   * If using only pickaxe-generic, environment.yml is recommended.

   Command to create new environment (replace env.yml with your chosen environment file)
   ```sh
   conda env create -f env.yml
   ```
   Existing environment (replace env.yml with your chosen environment file, and env-name with the existing environment name)
   ```sh
   conda env update -n env-name -f env.yml --prune
   ```
6. Activate the environment you installed rdkit into (replace env-name with your chosen environment name).  If a fresh environment was installed, the default environment name should be displayed on the screen.
   ```sh
   conda activate env-name
   ```
7. Install pickaxe-generic using pip.  Use a -e flag after "install" if you want your installation to update automatically when changing the files in this folder.  Otherwise, simply use the command below.
   ```sh
   python -m pip install pickaxe-generic
   ```
8. When running a Python program requiring pickaxe-generic, make sure you first open your terminal and activate the relevant environment using the command from Step 6.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Work in progress.  Check out example_notebook.ipynb in jupyter for examples.
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

Project Link: [https://github.com/wsprague-nu/Pickaxe-Generic](https://github.com/wsprague-nu/Pickaxe-Generic)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [Prof. Linda Broadbelt](broadbelt.northwestern.edu)
* Kevin Shebek
* Lauren Lopez
* Joseph Ni
* Jon Strutz
* [Prof. Brent Shanks](https://www.engineering.iastate.edu/people/profile/bshanks/)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/wsprague-nu/Pickaxe-Generic.svg?style=for-the-badge
[contributors-url]: https://github.com/wsprague-nu/Pickaxe-Generic/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/wsprague-nu/Pickaxe-Generic.svg?style=for-the-badge
[forks-url]: https://github.com/wsprague-nu/Pickaxe-Generic/network/members
[stars-shield]: https://img.shields.io/github/stars/wsprague-nu/Pickaxe-Generic.svg?style=for-the-badge
[stars-url]: https://github.com/wsprague-nu/Pickaxe-Generic/stargazers
[issues-shield]: https://img.shields.io/github/issues/wsprague-nu/Pickaxe-Generic.svg?style=for-the-badge
[issues-url]: https://github.com/wsprague-nu/Pickaxe-Generic/issues
[license-shield]: https://img.shields.io/github/license/wsprague-nu/Pickaxe-Generic.svg?style=for-the-badge
[license-url]: https://github.com/wsprague-nu/Pickaxe-Generic/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[product-screenshot]: images/screenshot.png