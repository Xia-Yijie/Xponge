from ... import *

CMAP = Generate_New_Bonded_Force_Type("atom_specific_cmap", "1-2-3-4-5", {"resolution":int, "parameters":list}, False)

@Molecule.Set_Save_SPONGE_Input
def write_cmap(self, prefix, dirname):
    bonds = []
    haved_cmaps = []
    haved_cmap_index = {}
    for bond in self.bonded_forces["atom_specific_cmap"]:
        order = list(range(5))
        if self.atom_index[bond.atoms[order[0]]] > self.atom_index[bond.atoms[order[-1]]]:
            temp_order = order[::-1]
        else:
            temp_order = order
        if bond.type not in haved_cmaps:            
            haved_cmap_index[bond.type] = len(haved_cmaps)
            haved_cmaps.append(bond.type)
        bonds.append("%d %d %d %d %d %d"%(self.atom_index[bond.atoms[temp_order[0]]]
        , self.atom_index[bond.atoms[temp_order[1]]], self.atom_index[bond.atoms[temp_order[2]]],
        self.atom_index[bond.atoms[temp_order[3]]], self.atom_index[bond.atoms[temp_order[4]]], haved_cmap_index[bond.type]))
    
    if (bonds):
        towrite = "%d %d\n"%(len(bonds), len(haved_cmaps))
        for bondtype in haved_cmaps:
            towrite += "%d "%bondtype.resolution
        towrite += "\n"
        for bondtype in haved_cmaps:
            for i, pi in enumerate(bondtype.parameters):
                towrite += "%f "%pi
                if (i+1) % bondtype.resolution == 0:
                    towrite += "\n"
            towrite += "\n"
        bonds.sort(key = lambda x: list(map(int, x.split()[:5])))
        towrite += "\n".join(bonds)
        
        f = open(os.path.join(dirname, prefix + "_cmap.txt"),"w")
        f.write(towrite)
        f.close()
