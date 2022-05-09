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

Xponge is a lightweight and easy-customizing Python package to perform pre- and post-processing of molecular simulations. It is mainly disigned for the MD program SPONGE [@Huang:2022], but it can also process some common format files so it is also helpful for some other simulation packages such as GROMACS [@Abraham:2015] and LAMMPS [@Thompson:2022]. Xponge includes three major categories of functionality, namely simulation system construction, simulation data transformation and analysis, and automated workflows for complex simulations. For the construction of simulation systems, Xponge can generate 3-dimensional molecular structures or read structures downloaded from online databases such as RCSB [@Berman:2000] and PubChem [@Sayers:2021], and then perform force field parameterization. The current force fields supported by Xponge contain CHARMM27 [@MacKerell:1998; @Mackerell:2004], ff14SB [@Maier:2015], ff19SB [@Tian:2020] for proteins, lipid14 [@Dickson:2014] and lipid17 for lipids, GAFF [@Wang:2004] for small organic molecules and tip3p, tip4pew, opc for water. For simulation data transformation and analysis, Xponge is combined with the Python package MDAnalysis [@Gowers:2016; @Michaud-Agrawal:2011] for processing. Xponge now also has an integrated workflow for free energy perturbation calculations with dual topology processing, enabling users to calculate hydration and binding free energies easily and foolproofly. At the same time, Xponge itself is highly modular and easily customizable, enabling simple extensions that mimic existing modules to develop one's own force fields, data analytics, and workflows.

# Statement of need

Since molecules are the basis of the world, molecular simulations are an important tool for understanding nature. For example, the computational simulation of proteins and small organic molecules in water can help with the drug designing. In addition to the simulation calculation itself, it is also important to construct and parameterize the simulation system and to process and analyze the simulation results. Because of the complexity and variety of systems to be simulated, it is often necessary to develop new and different force fields, approximations, algorithms, or analytical methods to obtain more accurate results.

There have been many pre-processing tools for molecular simulation, such as offline LEaP [@LEaP], pdb2gmx [@pdb2gmx] and psfgen [@psfgen], and online LigParGen [@LigParGen] and CGenFF [@CGenFF], but they are designed for a particular form of force field and are not suitable to modify. Xponge, on the other hand, is highly modular, making it well suited for development. At the same time, Xponge includes certain post-processing capabilities and combines specific complex pre-processing - simulation - post-processing into one workflow, making complex molecular simulations simple. Meanwhile, machine learning force field for molecular simulation is one of the hot research at present, and Xponge written by python can be better merged into the computational graph, which is more suitable for molecular simulation software like SPONGE [@Huang:2022] that can incorporate machine learning.

# Availability

Xponge is freely available and open source under the Apache License 2.0 (Apache-2.0). You can download the package and access the online documentations via https://gitee.com/gao_hyp_xyj_admin/xponge or https://github.com/Xia-Yijie/Xponge..

# Acknowledgements

The authors thank the National Key R&D Program of China (2017YFA0204702), the National Natural Science Foundation of China (21821004, 21927901, 92053202 and 22050003) for financial support.

# References
