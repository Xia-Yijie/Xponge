installation
-------------

To install ``Xponge``, the version of Python should be not less than 3.6.0

``Xponge``, the Python package itself is easy to install, and there are two ways to install.

.. TIP::
    
    Pay attention to the version of ``Xponge``.

.. NOTE::

    - ``Xponge`` installed from ``pip`` is a more stable version
    - ``Xponge`` in gitee or github Xponge repository is usually under development, it may have more functions but less stable.

I. Xponge as mindsponge.toolkits
========================================

Xponge now is part of MindScience as ``mindsponge.toolkits``. It can only be used in Linux at the moment, but you can install the MD package at the same time, which is very useful.

1. pip install
############################

Not Implemented yet, but coming soon!

2. source setup
############################

- 2.1 download or clone the source of the gitee repository

    The gitee repository is `MindScience <https://gitee.com/mindspore/mindscience>`_.
    
    You may need to choose a stable version instead of the head of the default branch in the repository.

- 2.2 open the directory where you download or clone the repository, and open the directory ``MindSponge``

- 2.3 run the command to build the wheel

.. code-block:: bash

    bash build.sh -e gpu -t on -j32

- 2.4 install the wheel

.. code-block:: bash

    pip install output/*.whl

II. Xponge as an independent package
========================================

Xponge can be an independent package. With the independent package, you can use it on Windows, and do not need to run a build.sh if you want to develop it.

1. pip install
############################

.. code-block:: bash

    pip install Xponge

2. source setup
############################

- 2.1 Download or clone the source from the gitee or github repository

    The gitee repository is `Xponge <https://gitee.com/gao_hyp_xyj_admin/xponge>`_
    The github repository is `here <https://github.com/xia-yijie/xponge>`_.
    
- 2.2 Open the directory where you download or clone the repository

- 2.3 (Optional) Configure the environment

    It is recommended to use `conda` to configure the environment. Two `yml` files named `install_requirements.yml` and `extras_requirements.yml` are provided in the repository.

    If you only need the basic functionality of Xponge, use the following command to configure the environment::

        conda env create -f install_requirements.yml
        conda activate Xponge

    If you want to experience all the features of Xponge, use the following command to configure the environment::

        conda env create -f extras_requirements.yml
        conda activate Xponge

- 2.4 run the command

.. code-block:: bash

    python setup.py install

.. code-block:: bash

    Xponge test -do base -o test
    
Here, ``Xponge`` can be replaced to ``python -m Xponge``, ``python3 -m Xponge`` or ``python -m mindsponge/toolkits`` and so on accordind to your settings of the environmental variables.
