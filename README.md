# Welcome to Use Xponge!

## Introduction

``Xponge`` is a lightweight and easy-customizing python package to perform pre- and post-processing of molecular simulations.

### What can Xponge do?

Xponge includes three major categories of functionality, namely, the simulation system construction, simulation data transformation and analysis, and automated workflows for complex simulations. ``Xponge`` is mainly designed for the molecular dynamics (MD) program [SPONGE](https://onlinelibrary.wiley.com/doi/epdf/10.1002/cjoc.202100456)[1], but it can also output some general format files such as mol2 and PDB, so it may help the other molecular modelling programs too.

## Installation

Xponge can be used on all operating systems (Windows/Linux/MacOS).

### 1. pip install

```bash
pip install Xponge
```

### 2. source setup

- 2.1 download or clone the source of the gitee or github repository

    The gitee repository is [here](https://gitee.com/gao_hyp_xyj_admin/xponge).
    The github repository is [here](https://github.com/xia-yijie/xponge).

        git clone http://gitee.com/gao_hyp_xyj_admin/xponge.git
        git clone http://github.com/xia-yijie/xponge.git

- 2.2 open the directory where you download or clone the repository

- 2.3 run the command

    ```bash
    python setup.py install
    ```

### Installation check

There are some unit tests in ``Xponge``. You can do the basic test to check whether the installation is successful like this:

```bash
Xponge test --do base -o test --verbose 1
```

Here, ``Xponge`` can be replaced to ``python -m Xponge``, ``python3 -m Xponge`` and so on according to your settings of the environmental variables. Some files will be generated after the test is finished.

## Usage example

Here is a simple example.

```python
import Xponge
# Import the force field you need
import Xponge.forcefield.amber.ff14sb
# Build the molecule like this
peptide = ACE + ALA + NME
# or like this
peptide2 = NALA + ALA * 10 + CALA
# or like this
peptide3 = Xponge.Get_Peptide_From_Sequence("AAAAA")
# See the documentation for more usage!
# Save them as your favorite format
Xponge.save_pdb(peptide, "ala.pdb")
Xponge.save_mol2(peptide2, "ala12.mol2")
Xponge.save_sponge_input(peptide3, "ala5")
```

## Detailed usage and API documentation

All can be seen [here](https://spongemm.cn/xponge_doc/index.html).

## Dependencies

`Xponge` does not depend on other packages except numpy for its basic use.

However, there are some complicated functions rely on some other packages. If you do not install the dependent package, you can not use the related functions.

Here is the list of all packages which may be uesd:

| package name      | description                       | how to install                 |
| ------------------| --------------------------------- | ------------------------------ |
| XpongeLib         | c/c++ compiled library for Xponge | `pip install XpongeLib`        |
| pyscf [2-4]       | quantum chemistry                 | `pip install pyscf`            |
| geometric[5]      | geometry optimization             | `pip install geometric`        |
| rdkit[6]          | cheminformatics                   | `conda install -c rdkit rdkit` |
| MDAnalysis[7-8]   | trajectory analysis               | `pip install MDAnalysis`       |
| mindspore[9]      | AI framework for machine learning | See the [official website](https://www.mindspore.cn/install)|
| mindsponge[1]     | end-to-end differentiable MD      | See the [official website](https://www.mindspore.cn/mindscience/docs/en/master/mindsponge/intro_and_install.html)|

## References

[1] Y.-P. Huang, et al. Chinese J. Chem. (2022) DOI: [10.1002/cjoc.202100456](https://doi.org/10.1002/cjoc.202100456)

[2] Q. Sun, et al. J. Chem. Phys. (2020) DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

[3] Q. Sun, et al. Wiley Interdiscip. Rev. Comput. Mol. Sci. (2018) DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)

[4] Q. Sun, J. Comp. Chem. (2015) DOI: [10.1002/jcc.23981](https://doi.org/10.1002/jcc.23981)

[5] L.-P. Wang, C.C. Song, J. Chem. Phys. (2016) DOI: [10.1063/1.4952956](https://doi.org/10.1063/1.4952956)

[6] RDKit: Open-source cheminformatics. https://www.rdkit.org

[7] R. J. Gowers, et al. Proceedings of the 15th Python in Science Conference (2016) DOI: [10.25080/majora-629e541a-00e](https://doi.org/10.25080/majora-629e541a-00e)

[8] N. Michaud-Agrawal, et al. J. Comput. Chem. (2011) DOI: [10.1002/jcc.21787](https://10.1002/jcc.21787)

[9] MindSpore: An Open AI Framwork. https://www.mindspore.cn/
