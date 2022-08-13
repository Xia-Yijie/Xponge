
.. TIP::
    
    You may need to change the commands and codes according to your environmental variables.

    For example, you may change this if you use Xponge as a part of mindscience::
    
        import Xponge
        import Xponge.forcefield.amber.ff14sb
    
    and

    .. code-block:: bash

        Xponge test --do base
        Xponge.mdrun SPONGE -mdin mdin.txt

    to::
    
        import mindsponge.toolkits as Xponge
        Xponge.source("mindsponge.toolkits.forcefield.amber.ff14sb")
        
    and

    .. code-block:: bash
        
        python -m mindsponge.toolkits test --do base
        python -m mindsponge.toolkits.mdrun SPONGE -mdin mdin.txt
