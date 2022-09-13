MD post-process
--------------------
.. include:: namespace.rst

Read the Mdout File
=======================

The mdout file is the file to record the basic information of the trajectory.

Here is an example named ``mdout.txt``::

            step         time  temperature    potential           LJ          PME      nb14_LJ      nb14_EE         bond        angle     dihedral         cmap 
            1000        1.000         0.56       878.58        -2.14       -59.32         4.85         6.10       933.95        11.50        14.34    -0.939009 
            2000        2.000         0.89       719.55        -2.22       -59.41         4.92         5.63       893.44        11.66        14.03     0.230595 
            3000        3.000         0.44       662.13        -0.84       -60.37         4.89         6.97       842.41        10.17        12.95    -0.985439 
            4000        4.000         0.52       631.07        -1.16       -60.47         6.95         6.42       793.67         9.22        13.30    -1.004319 
            5000        5.000         0.18       711.92        -0.95       -57.50         4.54         7.43       769.72         7.77        14.58    -2.345700 
            6000        6.000         0.80       706.59        -0.75       -60.08         4.76         6.68       741.36         7.02        13.63     0.140813 
            7000        7.000         0.23       712.97        -2.17       -55.32         4.16         3.43       720.23         4.65        13.95     0.428522 
            8000        8.000         0.37       757.94        -1.05       -59.79         5.17         6.89       704.59         9.16        15.01    -0.965497 
            9000        9.000         0.13       746.94        -1.57       -57.58         5.54         5.33       720.50         9.80        16.33    -2.496909 
           10000       10.000         0.43       758.65        -0.63       -57.45         5.38         4.12       743.64         9.09        13.59    -1.105110 

Then we can use ``MdoutReader`` in ``Xponge.analysis`` to read it and give the numpy arrays for the information::

    from Xponge.analysis import MdoutReader
    
    mdout = MdoutReader("mdout.txt")
    
    #directly print
    print(mdout.step)

    #use matplotlib to plot
    import matplotlib.pyplot as plt
    plt.plot(mdout.time, mdout.temperature)
    plt.show()

Process the Coordinates or Trajectory File
=========================================================

You can use ``SpongeCoordinateReader`` or ``SpongeTrajectoryReader`` in ``Xponge.analysis.md_analysis`` to read the coordinates or trajectory File::

    import Xponge.analysis.md_analysis as xmda
    import MDAnalysis as mda

    topology_file = ...
    trajectory_file = ...

    u = mda.Universe(topology_file, trajectory_file, format=xmda.SpongeTrajectoryReader)
    ...

You can use ``SpongeTrajectoryReader.with_arguments(box=box_file)`` to add the box file.

You can use ``SpongeCoordinateWriter`` or ``SpongeTrajectoryWriter`` in ``Xponge.analysis.md_analysis`` to write the coordinates or trajectory File::

    import Xponge.analysis.md_analysis as xmda
    import MDAnalysis as mda
    
    u = mda.Universe(...)

    with xmda.SpongeTrajectoryWriter(args.o) as w:
        for _ in u.trajectory:
            w.write(u)

Trajectory Analysis
====================

Xponge use the MDAnalysis to do the analysis. Here take the calculation of RDF as an example::

    import Xponge.analysis.md_analysis as xmda
    from MDAnalysis import Universe

    u = Universe("test.pdb", "mdcrd.dat", box="mdbox.txt", format=xmda.SpongeTrajectoryReader)

    #See the API of MDAnalysis for the following lines
    O = u.select_atoms("resname WAT and name O")
    from MDAnalysis.analysis import rdf as RDF
    rdf = RDF.InterRDF(O, O)
    rdf.run()
    import matplotlib.pyplot as plt
    plt.plot(rdf.results.bins[1:], rdf.results.rdf[1:])
    plt.show()

The final figure looks like

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/14.png