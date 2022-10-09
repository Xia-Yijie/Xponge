# Welcome to Use Xponge!

## Introduction

``Xponge`` is a lightweight and easy to customize python package to perform pre- and post-processing of molecular simulations.

### What can Xponge do?

Xponge includes three major categories of functionality, namely, the simulation system construction, simulation data transformation and analysis, and automated workflows for complex simulations. ``Xponge`` is mainly designed for the molecular dynamics (MD) program [SPONGE](https://onlinelibrary.wiley.com/doi/epdf/10.1002/cjoc.202100456)[1], but it can also output some general format files such as mol2 and PDB, so it may help the other molecular modelling programs too.

## Installation

Xponge can be used on all operating systems (Windows/Linux/MacOS). Some functions (See [here](https://spongemm.cn/xponge_doc/dependency.html#unavalable-functions-on-windows) for the detailed list) to do the quantum chemistry calculations can not be used on Windows because `pyscf` is not available on Windows.

### 1. pip install

```bash
pip install Xponge
```

### 2. source setup

- 2.1 Download or clone the source from the gitee or github repository

    The gitee repository is [here](https://gitee.com/gao_hyp_xyj_admin/xponge).
    The github repository is [here](https://github.com/xia-yijie/xponge).

        git clone http://gitee.com/gao_hyp_xyj_admin/xponge.git
        git clone http://github.com/xia-yijie/xponge.git

- 2.2 Open the directory where you download or clone the repository

- 2.3 (Optional) Configure the environment

    It is recommended to use `conda` to configure the environment. Two `yml` files named `install_requirements.yml` and `extras_requirements.yml` are provided in the repository.

    It is recommanded to use the file `install_requirements.yml` to configure the environment. The file will only install the basic dependent packages. If a `ModuleNotFoundError` is raised when you are using `Xponge`, then install the module. This allows you to avoid installing many modules that you will never use, and also makes `Xponge` more cross-platform compatible. Here are the commands to use `install_requirements.yml`.

    ```bash
    conda env create -f install_requirements.yml
    conda activate Xponge
    ```

    All the dependent packages are listed in the [dependencies](#dependencies) section. If you don't want to install the dependent packages one by one (which can be really annoying), the file `extras_requirements.yml` can help you with the environment configuration except the packages `mindspore` and `mindsponge`. The two packages should be installed according to your device (e.g. whether the backend is CPU, GPU or Huawei Ascend Processors) and can not be simply installed by conda. Here are the commands to use `extras_requirements.yml`.

    ```bash
    conda env create -f extras_requirements.yml
    conda activate Xponge
    ```

     It is worth noting that `extras_requirements.yml` can not be used on Windows because `pyscf` is not available on Windows.

- 2.4 Run the command

    ```bash
    python setup.py install
    ```

### Installation check

There are some unit tests in ``Xponge``. You can do the basic test to check whether the installation is successful like this:

```bash
Xponge test --do base -o test --verbose 1
```

Here, ``Xponge`` can be replaced to ``python -m Xponge``, ``python3 -m Xponge`` and so on according to your settings of the environmental variables. Some files will be generated after the test is finished.

## Quickstart

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

Then we can see `ala12.mol2` in VMD:

![pic2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/2.jpg)

Here is another simple example.

```python
import Xponge
import Xponge.forcefield.amber.tip3p

box = Xponge.BlockRegion(0, 0, 0, 60, 60, 60)
region_1 = Xponge.BlockRegion(0, 0, 20, 20, 20, 40)
region_2 = Xponge.BlockRegion(0, 0, 40, 20, 20, 60)
region_3 = Xponge.BlockRegion(0, 0, 0, 20, 20, 20)
region_4 = Xponge.SphereRegion(20, 10, 30, 10)
region_5 = Xponge.BlockRegion(20, 0, 20, 60, 60, 60)
region_2or3 = Xponge.UnionRegion(region_2, region_3)
region_4and5 = Xponge.IntersectRegion(region_4, region_5)
t = Xponge.Lattice("bcc", basis_molecule=CL, scale=4)
t2 = Xponge.Lattice("fcc", basis_molecule=K, scale=3)
t3 = Xponge.Lattice("sc", basis_molecule=NA, scale=3)
mol = t.Create(box, region_1)
mol = t2.create(box, region_2or3, mol)
mol = t3.create(box, region_4and5, mol)
Xponge.Save_PDB(mol, "out.pdb")
```

Then we can see `out.pdb` in VMD:

![pic1](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/1.jpg)

## Detailed usage and API documentation

All can be seen [here](https://spongemm.cn/xponge_doc/index.html).

## Contribution Guideline

If you want to contribute to the main codebase or report some issues, see [here](https://spongemm.cn/xponge_doc/contribution_guide.html) for the guides.

## Dependencies

`Xponge` does not depend on other packages except numpy for its basic use.

However, there are some complicated functions that depend on some other packages. If you do not install the dependent package, you can not use the related functions.

Here is the list of all packages which may be uesd:

| package name      | description                       | how to install                 |
| ------------------| --------------------------------- | ------------------------------ |
| XpongeLib         | c/c++ compiled library for Xponge | `pip install XpongeLib`        |
| pyscf [2-4]       | quantum chemistry                 | `pip install pyscf`            |
| geometric[5]      | geometry optimization             | `pip install geometric`        |
| rdkit[6]          | cheminformatics                   | `conda install -c rdkit rdkit` |
| MDAnalysis[7-8]   | trajectory analysis               | `pip install MDAnalysis`       |
| matplotlib        | plot and visualization            | `pip install matplotlib`       |
| mindspore[9]      | AI framework for machine learning | See the [official website](https://www.mindspore.cn/install)|
| mindsponge[1]     | end-to-end differentiable MD      | See the [official website](https://www.mindspore.cn/mindscience/docs/en/master/mindsponge/intro_and_install.html)|

## References

[0] Y. Xia, Y. Q. Gao, *J. Open Source Softw.* (2022) DOI:[10.21105/joss.04467](https://doi.org/10.21105/joss.04467)

[1] Y.-P. Huang, et al. *Chinese J. Chem.* (2022) DOI: [10.1002/cjoc.202100456](https://doi.org/10.1002/cjoc.202100456)

[2] Q. Sun, et al. *J. Chem. Phys.* (2020) DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

[3] Q. Sun, et al. Wiley Interdiscip. *Rev. Comput. Mol. Sci.* (2018) DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)

[4] Q. Sun, *J. Comp. Chem.* (2015) DOI: [10.1002/jcc.23981](https://doi.org/10.1002/jcc.23981)

[5] L.-P. Wang, C.C. Song, *J. Chem. Phys.* (2016) DOI: [10.1063/1.4952956](https://doi.org/10.1063/1.4952956)

[6] RDKit: Open-source cheminformatics. https://www.rdkit.org

[7] R. J. Gowers, et al. Proceedings of the 15th Python in Science Conference (2016) DOI: [10.25080/majora-629e541a-00e](https://doi.org/10.25080/majora-629e541a-00e)

[8] N. Michaud-Agrawal, et al. *J. Comput. Chem.* (2011) DOI: [10.1002/jcc.21787](https://10.1002/jcc.21787)

[9] MindSpore: An Open AI Framwork. https://www.mindspore.cn/
