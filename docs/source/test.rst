installation check
--------------------
.. include:: namespace.rst

There are some unit tests in ``Xponge``. You can do the basic test to check whether the installation is successful like this:

.. code-block:: bash

    Xponge test --do base -o test
    
Here, ``Xponge`` can be replaced to ``python -m Xponge``, ``python3 -m Xponge`` or ``python -m mindsponge.toolkits`` and so on according to your settings of the environmental variables.

You can use any visulization tool such as VMD or pymol to see the generated PDB file and the mol2 file to check the results. You can also compare all the files you get with the files in the test_standard directory in the gitee repository.

If you want to check the installation of the dependent packages, you can do the test for assignment:

.. code-block:: bash

    Xponge test --do assign -o test --verbose 1
    
Here are the list all the tests you can do:

.. list-table::
    :widths: 10 50
    :header-rows: 1
    :stub-columns: 1
    
    * - command
      - what to test
    * - base
      - the basic funtions of ``Xponge``
    * - assign
      - the assignment functions and dependency
    * - charmm27
      - the forcefield of charmm27
    * - all
      - all above

.. TIP::
    
    You can do multiple tests at the same time. For example::

        Xponge test --do base charmm27 -o test
    