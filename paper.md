---
title: 'Xponge: A Python package to perform pre- and post-processing of molecular simulations'
tags:
  - Python
  - molecular modelling
  - molecular simulation
  - molecular dynamics
  - pre-process and post-process
authors:
  - name: Yijie Xia
    orcid: 0000-0002-7931-4015
    affiliation: 1
  - name: Yi Qin Gao^[Corresponding author]
    affiliation: 1
affiliations:
 - name: College of Chemistry and Molecular Engineering, Peking University, China
   index: 1
date: 7 May 2022
bibliography: paper.bib
---

# Summary

<<<<<<< HEAD
Xponge is a lightweight and easy-customizing Python package to perform pre- and post-processing of molecular simulations. It is mainly disigned for the MD program SPONGE [@Huang:2022], but it can also process some common format files so it is also helpful for some other simulation packages such as GROMACS [@Abraham:2015] and LAMMPS [@Thompson:2022]. Xponge includes three major categories of functionality, namely simulation system construction, simulation data transformation and analysis, and automated workflows for complex simulations. For the construction of simulation systems, Xponge can generate 3-dimensional molecular structures or read structures downloaded from online databases such as RCSB [] and PubChem [], and then perform force field parameterization. The current force fields supported by Xponge contain CHARMM27 [], ff14SB [@Maier:2015], ff19SB [@Tian:2020] for proteins, lipid14 [] and lipid17 for lipids, GAFF [] for small organic molecules and tip3p, tip4pew, opc for water. For simulation data transformation and analysis, Xponge is combined with the Python package MDAnalysis [] for processing. Xponge now also has an integrated workflow for free energy perturbation calculations with dual topology processing, enabling users to calculate hydration and binding free energies easily and foolproofly. At the same time, Xponge itself is highly modular and easily customizable, enabling simple extensions that mimic existing modules to develop one's own force fields, data analytics, and workflows.
=======
Xponge is a lightweight and easy-customizing Python package to perform pre- and post-processing of molecular simulations. It is mainly disigned for the molecular dynamics (MD) program SPONGE [@Huang:2022], but the preprocessing performance can be helpful to other molecular simulation programs.
>>>>>>> ee18a480af2aeff01748e0f5fbb1224e4f4bda12

# Statement of need

Since molecules are the basis of the world, molecular simulations are an important tool for understanding nature. For example, the computational simulation of proteins and small organic molecules in water can help with the drug designing. In addition to the simulation calculation itself, it is also important to construct and parameterize the simulation system and to process and analyze the simulation results. Because of the complexity and variety of systems to be simulated, it is often necessary to develop new and different force fields, approximations, algorithms, or analytical methods to obtain more accurate results.

There have been many pre-processing tools for molecular simulation, such as offline antechamber [], pdb2gmx [] and psfgen [], and online LigParGen [] and CGenFF [], but they are designed for a particular form of force field and are not suitable for force field development. Xponge, on the other hand, is highly modular, making it well suited for force field development, especially the machine learning force field. At the same time, Xponge includes certain post-processing capabilities and combines specific complex pre-processing-simulation-post-processing into one workflow, making complex molecular simulations simple. Meanwhile, machine learning force field for molecular simulation is one of the hot research at present, and Xponge written by python can be better merged into the computational graph, which is more suitable for molecular simulation software like SPONGE [@Huang:2022] that can incorporate machine learning.

# Acknowledgements

The authors thank the National Key R&D Program of China (2017YFA0204702), the National Natural Science Foundation of China (21821004, 21873007 and 21927901) and CAAI-Huawei MindSpore Open Fund for financial support. This research was supported, in part, through the computational resources provided by the SZBL supercomputing center.

# References
