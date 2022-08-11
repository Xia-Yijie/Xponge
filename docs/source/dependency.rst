dependency
-------------

Although the basic usage of ``Xponge`` do not depend on other packages, some complicated functions (quantum mechanics calculation for example) rely on some other packages.

When an ``ImportError`` or a ``ModuleNotFoundError`` is raised, you need to install the module yourself.

Most of the module can be installed via ``pip``, while the ``RDKit`` package should be installed via ``conda``.

Here is the list of all packages which may be used:

.. list-table:: Xponge Dependency
    :widths: 10 20 20
    :header-rows: 1
    :stub-columns: 1
    
    * - package name
      - usage description
      - how to install
    * - mindspore
      - the basic package for mindsponge
      - ``pip install mindspore-gpu``
    * - mindsponge
      - the basic package for mindsponge.toolkits
      - ``pip install mindscience-sponge-gpu``
    * - XpongeLib
      - c/c++ compiled library for Xponge
      - ``pip install XpongeLib``
    * - pyscf
      - quantum chemistry
      - ``pip install pyscf``
    * - geometric
      - geometry optimization
      - ``pip install geometric``
    * - rdkit
      - cheminformatics
      - ``conda install -c rdkit rdkit``
    * - MDAnalysis
      - trajectory analysis
      - ``pip install MDAnalysis``
