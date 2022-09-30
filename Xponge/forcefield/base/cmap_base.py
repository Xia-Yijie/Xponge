"""
This **module** is the basic setting for the force field format of atom-specific cmap
"""
from ... import Generate_New_Bonded_Force_Type
from ...helper import Molecule, Xdict, set_global_alternative_names

# pylint: disable=invalid-name
CMapType = Generate_New_Bonded_Force_Type("cmap", "1-2-3-4-5", {"resolution": int, "parameters": list},
                                          False)

@CMapType.Set_Same_Force_Function
def cmap_same_force(_, atom_list):
    """
    This **function** is used to return the same force type for an atom list
    :param _:
    :param atom_list:
    :return:
    """
    return [atom_list]


@Molecule.Set_Save_SPONGE_Input("cmap")
def write_cmap(self):
    """
    This **function** is used to write SPONGE input file

    :param self: the Molecule instance
    :return: the string to write
    """
    bonds = []
    haved_cmaps = []
    haved_cmap_index = Xdict()
    for bond in self.bonded_forces.get("cmap", []):
        if bond.type not in haved_cmaps:
            haved_cmap_index[bond.type] = len(haved_cmaps)
            haved_cmaps.append(bond.type)
        bonds.append("%d %d %d %d %d %d" % (self.atom_index[bond.atoms[0]]
                                            , self.atom_index[bond.atoms[1]],
                                            self.atom_index[bond.atoms[2]],
                                            self.atom_index[bond.atoms[3]],
                                            self.atom_index[bond.atoms[4]], haved_cmap_index[bond.type]))

    if bonds:
        towrite = "%d %d\n" % (len(bonds), len(haved_cmaps))
        for bondtype in haved_cmaps:
            towrite += "%d " % bondtype.resolution
        towrite += "\n"
        for bondtype in haved_cmaps:
            for i, pi in enumerate(bondtype.parameters):
                towrite += "%f " % pi
                if (i + 1) % bondtype.resolution == 0:
                    towrite += "\n"
            towrite += "\n"
        bonds.sort(key=lambda x: list(map(int, x.split()[:5])))
        towrite += "\n".join(bonds)

        return towrite
    return None

set_global_alternative_names()
