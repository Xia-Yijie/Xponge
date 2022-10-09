dependencies
-------------

Although the basic usage of ``Xponge`` does not depend on other packages, some complicated functions (quantum mechanics calculation for example) rely on some other packages.

When an ``ImportError`` or a ``ModuleNotFoundError`` is raised, you need to install the module yourself.

Most of the module can be installed via ``pip``, while the ``RDKit`` package should be installed via ``conda``.

It is recommended to use ``conda`` to configure the environment. Two ``yml`` files named ``install_requirements.yml`` and ``extras_requirements.yml`` are provided in the repository.

It is recommanded to use the file `install_requirements.yml` to configure the environment. The file will only install the basic dependent packages. If a ``ModuleNotFoundError`` is raised when you are using ``Xponge``, then install the module. This allows you to avoid installing many modules that you will never use, and also makes `Xponge` more cross-platform compatible. Here are the commands to use ``install_requirements.yml``::

    conda env create -f install_requirements.yml
    conda activate Xponge

If you don't want to install the dependent packages one by one (which can be really annoying), the file ``extras_requirements.yml`` can help you with the environment configuration except the packages ``mindspore`` and ``mindsponge``. The two packages should be installed according to your device (e.g. whether the backend is CPU, GPU or Huawei Ascend Processors) and can not be simply installed by conda. Here are the commands to use ``extras_requirements.yml``::

    conda env create -f extras_requirements.yml
    conda activate Xponge

It is worth noting that ``extras_requirements.yml`` can not be used on Windows because ``pyscf`` is not available on Windows.

Unavalable functions on Windows
=================================

Here is the list of functions which are not available on Windows.

- Xponge.assign.calculate_charge(method="RESP")
- Xponge.assign.resp
- Xponge.assign.resp.resp_fit

The list of dependent packages
=================================

Here is the list of all packages which may be used:

.. list-table:: Xponge Dependency
    :widths: 10 20 20 10
    :header-rows: 1
    :stub-columns: 1
    
    * - package name
      - usage description
      - how to install
      - reference
    * - mindspore
      - AI framework for machine learning
      - See the `official website <https://www.mindspore.cn/install>`_
      - [1]
    * - mindsponge
      - the end-to-end differentiable MD
      - See the `official website <https://www.mindspore.cn/mindscience/docs/en/master/mindsponge/intro_and_install.html>`_
      - [2]
    * - XpongeLib
      - c/c++ compiled library for Xponge
      - ``pip install XpongeLib``
      - None
    * - pyscf
      - quantum chemistry
      - ``pip install pyscf``
      - [3-5]
    * - geometric
      - geometry optimization
      - ``pip install geometric``
      - [6]
    * - rdkit
      - cheminformatics
      - ``conda install -c rdkit rdkit``
      - [7]
    * - MDAnalysis
      - trajectory analysis
      - ``pip install MDAnalysis``
      - [8-9]
    * - matplotlib
      - plot and visualization
      - ``pip install matplotlib``
      -

References for the dependent packags

[1] MindSpore: An Open AI Framwork. https://www.mindspore.cn/

[2] Y.-P. Huang, et al. Chinese J. Chem. (2022) DOI: [10.1002/cjoc.202100456](https://doi.org/10.1002/cjoc.202100456)

[3] Q. Sun, et al. J. Chem. Phys. (2020) DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

[4] Q. Sun, et al. Wiley Interdiscip. Rev. Comput. Mol. Sci. (2018) DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)

[5] Q. Sun, J. Comp. Chem. (2015) DOI: [10.1002/jcc.23981](https://doi.org/10.1002/jcc.23981)

[6] L.-P. Wang, C.C. Song, J. Chem. Phys. (2016) DOI: [10.1063/1.4952956](https://doi.org/10.1063/1.4952956)

[7] RDKit: Open-source cheminformatics. https://www.rdkit.org

[8] R. J. Gowers, et al. Proceedings of the 15th Python in Science Conference (2016) DOI: [10.25080/majora-629e541a-00e](https://doi.org/10.25080/majora-629e541a-00e)

[9] N. Michaud-Agrawal, et al. J. Comput. Chem. (2011) DOI: [10.1002/jcc.21787](https://10.1002/jcc.21787)

