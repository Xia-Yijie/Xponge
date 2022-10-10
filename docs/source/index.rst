Xponge 1.2.6.11 Documentation
------------------------------------

|JOSS| |mdanalysis| |rdkit|

.. |JOSS| image:: https://joss.theoj.org/papers/10.21105/joss.04467/status.svg
    :alt: JOSS
    :target: https://doi.org/10.21105/joss.04467

.. |mdanalysis| image:: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
    :alt: Powered by MDAnalysis
    :target: https://www.mdanalysis.org

.. |rdkit| image:: https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC
    :target: https://www.rdkit.org/
    :alt: Powered by RDKit

What is Xponge?
    ``Xponge`` is a lightweight and easy to customize python package to perform pre- and post-processing of molecular simulations. 

What can Xponge do?
    ``Xponge`` is mainly designed for the molecular dynamics (MD) program `SPONGE <https://onlinelibrary.wiley.com/doi/epdf/10.1002/cjoc.202100456>`_, but it can also output some general format files such as mol2 and PDB, so it may help the other molecular modelling programs too.

How can I get Xponge?
    Xponge is available on all operating systems (Windows/Linux/MacOS). See :ref:`installation` for detailed instructions.

    ``Xponge`` now is also a part of `MindSponge <https://gitee.com/mindspore/mindscience>`_ as toolkits. You can download and use ``Xponge`` with MindSponge.


.. toctree::
    :maxdepth: 2
    :caption: START 

    installation
    dependency
    test
 
.. toctree::
    :maxdepth: 2
    :caption: BASIC USAGE 

    building_systems
    type_assignment
    structure_pre-process
    MD_post-process
    command_line_tools
    
.. toctree::
    :maxdepth: 2
    :caption: ADVANCE USAGE 

    developing_forcefield
    automated_workflow

.. toctree::
    :maxdepth: 2
    :caption: TUTORIAL

    covid bound to b6 <https://gitee.com/mindspore/mindscience/blob/master/MindSPONGE/applications/molecular_dynamics/tradition/covid_bound_to_b6/covid_bound_to_b6.ipynb>
    the folding of ALA12 <https://gitee.com/mindspore/mindscience/blob/master/MindSPONGE/applications/molecular_dynamics/tradition/the_folding_of_ala12/the_folding_of_ala12.ipynb>
    building of KALP15 <https://gitee.com/xiayijie/mindscience/blob/master/MindSPONGE/applications/molecular_dynamics/building_of_KALP15/building_of_KALP15-Copy1.ipynb>
    relative binding energy of FXR <https://gitee.com/mindspore/mindscience/blob/master/MindSPONGE/applications/molecular_dynamics/tradition/relative_binding_energy_of_FXR/relative_binding_energy_of_FXR.ipynb>
    tutorials on Gitee (in Chinese) <https://gitee.com/gao_hyp_xyj_admin/sponge_tutorial>

.. toctree::
    :maxdepth: 2
    :caption: API

    modules

.. toctree::
    :maxdepth: 2
    :caption: OTHERS
    
    contribution_guide
    official website of SPONGE <https://spongemm.cn>
    github repository of Xponge <https://github.com/Xia-Yijie/xponge>
    gitee repository of Xponge <https://gitee.com/gao_hyp_xyj_admin/xponge>
    gitee repository of MindScience <https://gitee.com/mindspore/mindscience>

Indices and tables
------------------------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
