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
date: 22 May 2022
bibliography: paper.bib
---

# Summary

Xponge is a lightweight and easy-customizing Python package to perform pre- and post-processing of molecular simulations. It is mainly designed for the MD program SPONGE [@Huang:2022], but it can also process common format files and therefore it should also be useful for other simulation packages such as GROMACS [@Abraham:2015] and LAMMPS [@Thompson:2022]. Xponge includes three major categories of functionality, namely, the simulation system construction, simulation data transformation and analysis, and automated workflows for complex simulations. For the construction of simulation systems, Xponge can generate 3-dimensional molecular structures or read structures downloaded from online databases such as RCSB [@Berman:2000] and PubChem [@Sayers:2021], and perform force field parameterization. The current force fields supported by Xponge include CHARMM27 [@MacKerell:1998; @Mackerell:2004], ff14SB [@Maier:2015], ff19SB [@Tian:2020] for proteins, lipid14 [@Dickson:2014] and lipid17 for lipids, GAFF [@Wang:2004] for small organic molecules and tip3p, tip4pew, opc for water. For simulation data transformation and analysis, Xponge is combined with the Python package MDAnalysis [@Gowers:2016; @Michaud-Agrawal:2011] for processing. Xponge has an integrated workflow for free energy perturbation calculations with dual topology construction, facilitating calculations on hyd	ration and binding free energies. At the same time, Xponge itself is highly modular and easily customizable, enabling simple extensions that mimic existing modules to develop one's own force fields, data analytics, and workflows.

# Statement of need

Molecular simulation is becoming an important and useful tool in many different research areas. For example, the computational simulation of proteins and small organic molecules in water is now widely used in drug design. In MD simulations, it is important to construct and parameterize the simulation system and to process and analyze the simulation results. Because of the complexity and variety of systems to be simulated, it is often necessary to develop new and different force fields, approximations, algorithms, or analytical methods to obtain and visualize results.
There exist a variety of pre-processing tools for molecular simulations, such as offline LEaP [@LEaP], pdb2gmx [@pdb2gmx] and psfgen [@psfgen], and online LigParGen [@LigParGen] and CGenFF [@CGenFF], but they are designed for a particular form of force field and can be difficult to be customized. Xponge, is designed to be highly modular and intended to be developer-friendly. At the same time, Xponge includes a number of post-processing capabilities and combines specific complex pre-processing - simulation - post-processing into one single workflow. Meanwhile, since machine learning force field for molecular simulation represents one of the current research frontiers, Xponge is written by python and thus can be better merged into the computational graph. Such a setup fits well with molecular simulation softwares such as SPONGE [@Huang:2022] which incorporates machine learning methods.

# Availability

Xponge is freely available and open source under the Apache License 2.0 (Apache-2.0). You can download the package and access the online documentations on [gitee](https://gitee.com/gao_hyp_xyj_admin/xponge) or [github](https://github.com/Xia-Yijie/Xponge).

# Acknowledgements

The authors thank the National Key R&D Program of China (2017YFA0204702), the National Natural Science Foundation of China (21821004, 21927901, 92053202 and 22050003) for financial support.

# References
