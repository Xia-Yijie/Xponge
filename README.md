# Welcome to Use Xponge!

## What is Xponge?

``Xponge`` is a lightweight and easy-customizing python package to perform pre- and post-processing of molecular simulations.

## What can Xponge do?

Xponge includes three major categories of functionality, namely, the simulation system construction, simulation data transformation and analysis, and automated workflows for complex simulations. ``Xponge`` is mainly designed for the molecular dynamics (MD) program [SPONGE](https://onlinelibrary.wiley.com/doi/epdf/10.1002/cjoc.202100456), but it can also output some general format files such as mol2 and PDB, so it may help the other molecular modelling programs too.

## How can I get Xponge?

Xponge can be used on both Windows and Linux operating systems

### 1. pip install

```bash
pip install Xponge
```

### 2. source setup

- 2.1 download or clone the source of the gitee or github repository

    The gitee repository is [here](https://gitee.com/gao_hyp_xyj_admin/xponge).
    The github repository is [here](https://github.com/xia-yijie/xponge).

- 2.2 open the directory where you download or clone the repository

- 2.3 run the command

    ```bash
    python setup.py install
    ```

## How can I check whether I have installed Xponge correctly?

There are some unit tests in ``Xponge``. You can do the basic test to check whether the installation is successful like this:

```bash
Xponge test -do base -o test
```

Here, ``Xponge`` can be replaced to ``python -m Xponge``, ``python3 -m Xponge`` and so on according to your settings of the environmental variables.

## Can you give me a example to use?

Here is a simple example.

```python
import Xponge
# import the force field you need
import Xponge.forcefield.amber.ff14sb
# build the molecule
peptide = ACE + ALA + NME
# save it as your favorite format
Xponge.save_pdb(peptide, "ala.pdb")
Xponge.save_mol2(peptide, "ala.mol2")
Xponge.save_sponge_input(peptide, "ala")
```

## Where can I see the complicated usage, the API documentation or the way to develop my own force field?

All can be seen [here](https://spongemm.cn/xponge_doc/index.html).
