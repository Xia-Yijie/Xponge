from ... import *
from ..BASE import CHARGE, MASS, LJ
from ..BASE import NB14, NB14_EXTRA, EXCLUDE
from ..BASE import BOND, SOFT_BOND

LJType = LJ.LJType

LJType.New_From_String("""name A  B
ZERO_LJ_ATOM-ZERO_LJ_ATOM  0  0
""")

AtomType.Add_Property({"LJtypeB": str})
AtomType.Add_Property({"subsys": int})

def _find_common_forces(forcetype, Aforces, Bforces, molB2molA):
    toret = []
    temp_map = {}
    temp_map2 = {}
    for force in Bforces:
        temp_map2[force] = True
        for fatoms in forcetype.Same_Force(force.atoms):
            temp_map["-".join(list(map(lambda atom: str(molB2molA[atom]), fatoms)))] = force
    for force in Aforces:
        tofind = "-".join(list(map(str, force.atoms)))
        Bforce = temp_map.get(tofind, None)
        toret.append([force, Bforce])
        temp_map2[Bforce] = False
    for force in Bforces:
        if temp_map2[force]:
            toret.append([None, force])
    return toret


FEP_Bonded_Force_Merge_Rule = {}


def nb14_extra_merge_rule(molR, molA, molB, forcetype, Rforces, Bforces, _lambda, molR2molA, molA2molR, molR2molB,
                          molB2molR):
    TINY = 1e-10 / 18.2223
    forcepair = _find_common_forces(forcetype, Rforces, Bforces, molB2molR)
    for fR, fB in forcepair:
        if fB is None:
            temp_charge0 = fR.atoms[0].charge if abs(fR.atoms[0].charge) > TINY else TINY
            temp_charge1 = fR.atoms[1].charge if abs(fR.atoms[1].charge) > TINY else TINY

            if fR.nb14_ee_factor is not None:
                nb14_ee_factor = fR.nb14_ee_factor
            else:
                nb14_ee_factor = fR.kee * molR2molA[fR.atoms[0]].charge * \
                                 molR2molA[fR.atoms[1]].charge
            fR.kee = nb14_ee_factor / temp_charge0 / temp_charge1

            fR.kee *= 1 - _lambda
            fR.A *= 1 - _lambda
            fR.B *= 1 - _lambda
        elif fR is None:
            fR = NB14_EXTRA.NB14Type.entity(list(map(lambda x: molB2molR[x], fB.atoms)), fB.type, fB.name)

            temp_charge0 = fR.atoms[0].charge if abs(fR.atoms[0].charge) > TINY else TINY
            temp_charge1 = fR.atoms[1].charge if abs(fR.atoms[1].charge) > TINY else TINY

            if fB.nb14_ee_factor is not None:
                nb14_ee_factor = fB.nb14_ee_factor
            else:
                nb14_ee_factor = fB.kee * fB.atoms[0].charge * \
                                 fB.atoms[1].charge
            fR.kee = nb14_ee_factor / temp_charge0 / temp_charge1

            fR.A = fB.A * _lambda
            fR.B = fB.B * _lambda
            fR.kee = f2.kee * _lambda
            molR.Add_Bonded_Force(fR)
        else:
            temp_charge0 = fR.atoms[0].charge if abs(fR.atoms[0].charge) > TINY else TINY
            temp_charge1 = fR.atoms[1].charge if abs(fR.atoms[1].charge) > TINY else TINY

            if fR.nb14_ee_factor is not None:
                nb14_ee_factor = fR.nb14_ee_factor
            else:
                nb14_ee_factor = fR.kee * molR2molA[fR.atoms[0]].charge * \
                                 molR2molA[fR.atoms[1]].charge
            fR.kee = nb14_ee_factor / temp_charge0 / temp_charge1

            if fB.nb14_ee_factor is not None:
                nb14_ee_factor = fB.nb14_ee_factor
            else:
                nb14_ee_factor = fB.kee * fB.atoms[0].charge * \
                                 fB.atoms[1].charge
            kee = nb14_ee_factor / temp_charge0 / temp_charge1
            
            fR.kee = fR.kee * (1 - _lambda) + kee * _lambda
            fR.A = fR.A * (1 - _lambda) + fB.A * _lambda
            fR.B = fR.B * (1 - _lambda) + fB.B * _lambda


FEP_Bonded_Force_Merge_Rule["nb14_extra"] = {"lambda_name": "dihedral", "merge_function": nb14_extra_merge_rule}


def bond_merge_rule(molR, molA, molB, forcetype, Rforces, Bforces, _lambda, molR2molA, molA2molR, molR2molB, molB2molR):
    forcepair = _find_common_forces(forcetype, Rforces, Bforces, molB2molR)
    for fR, fB in forcepair:
        if fB is None:
            molR.bonded_forces["bond"].remove(fR)
            fR.from_AorB = 0
            molR.Add_Bonded_Force(fR, "soft_bond")
        elif fR is None:
            fR = BOND.BondType.entity(list(map(lambda x: molB2molR[x], fB.atoms)), fB.type, fB.name)
            fR.k = fB.k
            fR.b = fB.b
            fR.from_AorB = 1
            molR.Add_Bonded_Force(fR, "soft_bond")
        elif abs(fR.b - fB.b) < 1e-5:
            fR.k = fR.k * (1 - _lambda) + fB.k * _lambda
        else:
            fR2 = BOND.BondType.entity(list(map(lambda x: molB2molR[x], fB.atoms)), fB.type, fB.name)
            fR2.k = fB.k * _lambda
            fR.b = fB.b
            fR.k *= (1 - _lambda)
            molR.Add_Bonded_Force(fR2)


FEP_Bonded_Force_Merge_Rule["bond"] = {"lambda_name": "bond", "merge_function": bond_merge_rule}

from ..BASE import ANGLE


def angle_merge_rule(molR, molA, molB, forcetype, Rforces, Bforces, _lambda, molR2molA, molA2molR, molR2molB,
                     molB2molR):
    forcepair = _find_common_forces(forcetype, Rforces, Bforces, molB2molR)
    for fR, fB in forcepair:
        if fB is None:
            fR.k *= 1 - _lambda
        elif fR is None:
            fR = ANGLE.AngleType.entity(list(map(lambda x: molB2molR[x], fB.atoms)), fB.type, fB.name)
            fR.k *= _lambda
            fR.b = fB.b
            molR.Add_Bonded_Force(fR)
        elif abs(fR.b - fB.b) < 1e-5:
            fR.k = fR.k * (1 - _lambda) + fB.k * _lambda 
        else:
            fR2 = ANGLE.AngleType.entity(list(map(lambda x: molB2molR[x], fB.atoms)), fB.type, fB.name)
            fR2.k = fB.k * _lambda
            fR2.b = fB.b
            fR.k *= (1 - _lambda)
            molR.Add_Bonded_Force(fR2)


FEP_Bonded_Force_Merge_Rule["angle"] = {"lambda_name": "angle", "merge_function": angle_merge_rule}

from ..BASE import DIHEDRAL


def dihedral_merge_rule(molR, molA, molB, forcetype, Rforces, Bforces, _lambda, molR2molA, molA2molR, molR2molB,
                        molB2molR):
    forcepair = _find_common_forces(forcetype, Rforces, Bforces, molB2molR)
    for fR, fB in forcepair:
        if fB is None:
            for i in range(fR.multiple_numbers):
                fR.ks[i] *= 1 - _lambda
        elif fR is None:
            fR2 = DIHEDRAL.ProperType.entity(list(map(lambda x: molB2molR[x], fB.atoms)), fB.type, fB.name)
            fR2.multiple_numbers = fB.multiple_numbers
            for i in range(fB.multiple_numbers):
                fR2.ks.append(fB.ks[i] * _lambda)
                fR2.phi0s.append(fB.phi0s[i])
                fR2.periodicitys.append(fB.periodicitys[i])
            molR.Add_Bonded_Force(fR2)
        else:
            sameforce = fR.multiple_numbers == fB.multiple_numbers
            if sameforce:
                check_map = {}
                for i in range(fB.multiple_numbers):
                    check_map[fB.periodicitys[i]] = fB.phi0s[i]
                for i in range(fR.multiple_numbers):
                    if abs(check_map.get(fR.periodicitys[i], float("Inf")) - fR.phi0s[i]) > 1e-5:
                        sameforce = False
                        break
            if sameforce:
                for i in range(fR.multiple_numbers):
                    fR.ks[i] = fR.ks[i] * (1 - _lambda) + fB.ks[i] * _lambda
            else:
                for i in range(fR.multiple_numbers):
                    fR.ks[i] *= (1 - _lambda)
                fR2 = DIHEDRAL.ProperType.entity(list(map(lambda x: molB2molR[x], fB.atoms)), fB.type, fB.name)
                fR2.multiple_numbers = fB.multiple_numbers
                for i in range(fB.multiple_numbers):
                    fR2.ks.append(fB.ks[i] * _lambda)
                    fR2.phi0s.append(fB.phi0s[i])
                    fR2.periodicitys.append(fB.periodicitys[i])
                molR.Add_Bonded_Force(fR2)


FEP_Bonded_Force_Merge_Rule["dihedral"] = {"lambda_name": "dihedral", "merge_function": dihedral_merge_rule}


def improper_merge_rule(molR, molA, molB, forcetype, Rforces, Bforces, _lambda, molR2molA, molA2molR, molR2molB,
                        molB2molR):
    forcepair = _find_common_forces(forcetype, Rforces, Bforces, molB2molR)
    for fR, fB in forcepair:
        if fB is None:
            fR.k *= 1 - _lambda
        elif fR is None:
            fR = DIHEDRAL.ImproperType.entity(list(map(lambda x: molB2molR[x], fB.atoms)), fB.type, fB.name)
            fR.k = fB.k * _lambda
            fR.phi0 = fB.phi0
            fR.periodicity = fB.periodicity
            molR.Add_Bonded_Force(fR)
        elif abs(fR.phi0 - fB.phi0) < 1e-5 and fR.periodicity == fB.periodicity:
            fR.k = fR.k * (1 - _lambda) + fB.k * _lambda
        else:
            fR2 = DIHEDRAL.ImproperType.entity(list(map(lambda x: molB2molR[x], fB.atoms)), fB.type, fB.name)
            fR2.k = fB.k * _lambda
            fR2.phi0 = fB.phi0
            fR2.periodicity = fB.periodicity
            fR.k *= (1 - _lambda)
            molR.Add_Bonded_Force(fR2)


FEP_Bonded_Force_Merge_Rule["improper"] = {"lambda_name": "dihedral", "merge_function": improper_merge_rule}


def Save_Hard_Core_LJ():
    Molecule.Set_Save_SPONGE_Input("LJ")(LJ.write_LJ)
    Molecule.Del_Save_SPONGE_Input("LJ_soft_core")
    Molecule.Del_Save_SPONGE_Input("subsys_division")

def Save_Soft_Core_LJ():
    Molecule.Del_Save_SPONGE_Input("LJ")

    @Molecule.Set_Save_SPONGE_Input("subsys_division")
    def write_subsys_division(self):
        towrite = "%d\n"%len(self.atoms)
        for atom in self.atoms:
            if getattr(atom, "subsys", None) is None:
                towrite += "%d\n"%0
            else:
                towrite += "%d\n"%atom.subsys
        return towrite

    @Molecule.Set_Save_SPONGE_Input("LJ_soft_core")
    def write_LJ(self):
        LJtypes = []
        LJtypemap = {}
        LJtypesB = []
        LJtypemapB = {}
        for atom in self.atoms:
            if atom.LJtype not in LJtypemap.keys():
                LJtypemap[atom.LJtype] = len(LJtypes)
                LJtypes.append(atom.LJtype)
            if atom.LJtypeB is None:
                atom.LJtypeB = atom.LJtype
            if atom.LJtypeB not in LJtypemapB.keys():
                LJtypemapB[atom.LJtypeB] = len(LJtypesB)
                LJtypesB.append(atom.LJtypeB)
        
        As, Bs = LJ.find_AB_LJ(LJtypes)
        AsA, BsB = LJ.find_AB_LJ(LJtypesB)

        checks = LJ.get_checks(LJtypes, As, Bs)
        same_type = LJ.judge_same_type(LJtypes, checks)
        real_LJtypes = LJ.get_real_LJ(LJtypes, same_type)
        real_As, real_Bs = LJ.find_AB_LJ(real_LJtypes)

        checks = LJ.get_checks(LJtypesB, AsA, BsB)
        same_typeB = LJ.judge_same_type(LJtypesB, checks)
        real_LJtypesB = LJ.get_real_LJ(LJtypesB, same_typeB)
        real_AsB, real_BsB = LJ.find_AB_LJ(real_LJtypesB)
        
        towrite = "%d %d %d\n\n" % (len(self.atoms), len(real_LJtypes), len(real_LJtypesB))
        count = 0
        for i in range(len(real_LJtypes)):
            for j in range(i + 1):
                towrite += "%16.7e" % real_As[count] + " "
                count += 1
            towrite += "\n"
        towrite += "\n"

        count = 0
        for i in range(len(real_LJtypes)):
            for j in range(i + 1):
                towrite += "%16.7e" % real_Bs[count] + " "
                count += 1
            towrite += "\n"
        towrite += "\n"

        count = 0
        for i in range(len(real_LJtypesB)):
            for j in range(i + 1):
                towrite += "%16.7e" % real_AsB[count] + " "
                count += 1
            towrite += "\n"
        towrite += "\n"

        count = 0
        for i in range(len(real_LJtypesB)):
            for j in range(i + 1):
                towrite += "%16.7e" % real_BsB[count] + " "
                count += 1
            towrite += "\n"

        towrite += "\n"
        towrite += "\n".join(
            ["%d %d" % (same_type[LJtypemap[atom.LJtype]], same_typeB[LJtypemapB[atom.LJtypeB]]) for atom in
             self.atoms])
        return towrite


def Intramolecule_NB_To_NB14(molA, perturbing_residues):
    if isinstance(perturbing_residues, Residue):
        perturbing_residues = [perturbing_residues]
    BUILD.Build_Bonded_Force(molA)
    A_Exclude = EXCLUDE.Exclude.current.Get_Excluded_Atoms(molA)
    for residue1 in perturbing_residues:
        for atomA1 in residue1.atoms:
            for residue2 in perturbing_residues:
                for atomA2 in residue2.atoms:
                    if atomA1 == atomA2:
                        continue
                    if atomA2 not in A_Exclude[atomA1]:
                        A, B = NB14_EXTRA.Get_NB14EXTRA_AB(atomA1, atomA2)
                        new_force = NB14_EXTRA.NB14Type.entity([atomA1, atomA2], NB14_EXTRA.NB14Type.types["UNKNOWNS"])
                        new_force.A = A
                        new_force.B = B
                        new_force.kee = 1
                        molA.Add_Bonded_Force(new_force)

                        atomA1.Extra_Exclude_Atom(atomA2)
                        A_Exclude[atomA1].add(atomA2)
                        A_Exclude[atomA2].add(atomA1)


def Get_Free_Molecule(molA, perturbing_residues, intra_FEP=False):
    if isinstance(perturbing_residues, Residue):
        perturbing_residues = [perturbing_residues]
    BUILD.Build_Bonded_Force(molA)
    NB14_EXTRA.NB14_To_NB14EXTRA(molA)
    molB = molA.deepcopy()
    molA2molB = {}
    for i, atomA in enumerate(molA.atoms):
        molA2molB[atomA] = molB.atoms[i]

    for residue in perturbing_residues:
        for atomA in residue.atoms:
            atom = molA2molB[atomA]
            atom.charge = 0
            atom.LJtype = "ZERO_LJ_ATOM"
            atom.subsys = 1
            atomA.subsys = 1

    return molB


def Merge_Dual_Topology(mol, ResidueA, ResidueB, AssignA, AssignB, tmcs = 60):
    
    BUILD.Build_Bonded_Force(mol)
    BUILD.Build_Bonded_Force(ResidueB)
    
    from ...assign.RDKit_tools import Assign2RDKitMol, Get_Part_Align, \
        Set_Conformer_Coordinate_From_Residue, Get_Conformer_Coordinate_To_Residue, \
        Get_Conformer_Coordinate, Insert_Atom_Type_To_RDKitMol
    from rdkit import Chem
    from rdkit.Chem import rdFMCS as MCS

    RDmolA = Assign2RDKitMol(AssignA, True)
    RDmolB = Assign2RDKitMol(AssignB, True)

    atom_type_dict = {}
    Insert_Atom_Type_To_RDKitMol(RDmolA, ResidueA, AssignA, atom_type_dict)
    Insert_Atom_Type_To_RDKitMol(RDmolB, ResidueB, AssignB, atom_type_dict)
    print("FINDING MAXIMUM COMMON SUBSTRUCTURE")

    result = MCS.FindMCS([RDmolA, RDmolB], atomCompare=MCS.AtomCompare.CompareIsotopes, completeRingsOnly=True, timeout = tmcs)
    #print(result.smartsString)
    #RDmol_mcs = Chem.MolFromSmarts(result.smarts)
    RDmol_mcs = result.queryMol
    #from rdkit.Chem import Draw
    #Draw.MolToImageFile(RDmol_mcs, "test.jpg")
    matchA = RDmolA.GetSubstructMatch(RDmol_mcs)
    matchB = RDmolB.GetSubstructMatch(RDmol_mcs)
    matchmap = {matchB[j]: matchA[j] for j in range(len(matchA))}
  
    Set_Conformer_Coordinate_From_Residue(RDmolA, ResidueA, AssignA)
    Set_Conformer_Coordinate_From_Residue(RDmolB, ResidueB, AssignB)
    print("ALIGNING TOPOLOGY AND COORDINATE")
    Get_Part_Align(RDmolA, RDmolB, matchA, matchB)

    ResidueTypeA = ResidueA.type

    if isinstance(ResidueB, Residue):
        ResidueTypeB = ResidueB.type
    elif isinstance(ResidueB, ResidueType):
        ResidueTypeB = ResidueB
    else:
        raise TypeError

    Get_Conformer_Coordinate_To_Residue(RDmolB, ResidueTypeB, AssignB)
    
    forcopy = hash(str(time.time()))
    restypeAB = ResidueTypeA.deepcopy(ResidueTypeA.name + "_" + ResidueTypeB.name, forcopy)

    extraA = []
    extraB = []
    RBmap = {value: key for key, value in matchmap.items()}
    for i in range(len(ResidueTypeA.atoms)):
        atom = ResidueTypeA.atoms[i]
        atom.x = ResidueA._name2atom[atom.name].x
        atom.y = ResidueA._name2atom[atom.name].y
        atom.z = ResidueA._name2atom[atom.name].z
        if i not in matchA:
            extraA.append(atom.copied[forcopy])            
            extraA[-1].subsys = 1
            
    for i in range(len(ResidueTypeB.atoms)):
        if i not in matchB:
            RBmap[len(restypeAB.atoms)] = i
            atom = ResidueTypeB.atoms[i]
            restypeAB.Add_Atom(atom.name + "R2", atom.type, atom.x, atom.y, atom.z)
            atom.copied[forcopy] = restypeAB.atoms[-1]
            atom.copied[forcopy].contents = {key: value for key, value in atom.contents.items()}
            atom.copied[forcopy].name = atom.name + "R2"
            extraB.append(atom.copied[forcopy])
            extraB[-1].subsys = 2
        else:
            ResidueTypeB.atoms[i].copied[forcopy] = restypeAB.atoms[matchmap[i]]

    for atomi in extraA:
        for atomj in extraB:
            atomi.Extra_Exclude_Atom(atomj)

    for atom, connect_set in ResidueTypeB.connectivity.items():
        for aton in connect_set:
            restypeAB.Add_Connectivity(atom.copied[forcopy], aton.copied[forcopy])
    

    for bond_entities in ResidueTypeB.bonded_forces.values():
        for bond_entity in bond_entities:
            tocopy = False
            for atom in bond_entity.atoms:
                if ResidueTypeB._atom2index[atom] not in matchmap.keys():
                    tocopy = True
                    break
            if tocopy:
                restypeAB.Add_Bonded_Force(bond_entity.deepcopy(forcopy))
    
                          

    for atom in ResidueTypeB.atoms:
        for key, linked_atoms in atom.copied[forcopy].linked_atoms.items():
            for aton in atom.linked_atoms.get(key, []):
                if not (ResidueTypeB._atom2index[aton] in matchmap.keys()
                        and ResidueTypeB._atom2index[atom] in matchmap.keys()):
                    atom.copied[forcopy].Link_Atom(key, aton.copied[forcopy])
    
    for atom in ResidueTypeA.atoms:
        atom.copied.pop(forcopy)

    for atom in ResidueTypeB.atoms:
        atom.copied.pop(forcopy)

    restypeBA = restypeAB.deepcopy(ResidueTypeB.name + "_" + ResidueTypeA.name)

    BUILD.Build_Bonded_Force(restypeBA)
    BUILD.Build_Bonded_Force(restypeAB)
    
    NB14_EXTRA.NB14_To_NB14EXTRA(restypeBA)
    NB14_EXTRA.NB14_To_NB14EXTRA(restypeAB)
    for i in range(len(restypeAB.atoms)):
        if i < len(ResidueTypeA.atoms):
            restypeAB.atoms[i].contents.update({key: value for key, value in ResidueTypeA.atoms[i].contents.items() if key != "name" })
        else:
            restypeAB.atoms[i].LJtype = "ZERO_LJ_ATOM"
            restypeAB.atoms[i].charge = 0
            restypeAB.atoms[i].subsys = 2
            restypeBA.atoms[i].subsys = 2

        if i in RBmap:
            restypeBA.atoms[i].contents.update({key: value for key, value in ResidueTypeB.atoms[RBmap[i]].contents.items() if key != "name" })
        else:
            restypeBA.atoms[i].LJtype = "ZERO_LJ_ATOM"
            restypeBA.atoms[i].charge = 0
            restypeAB.atoms[i].subsys = 1
            restypeBA.atoms[i].subsys = 1

    molA = Molecule(mol.name + "A")
    molB = Molecule(mol.name + "B")

    for res in mol.residues:
        if res == ResidueA:
            molA.Add_Residue(restypeAB)
            molB.Add_Residue(restypeBA)
        else:
            molA.Add_Residue(res)
            molB.Add_Residue(res)

    for reslink in mol.residue_links:
        molA.residue_links.append(reslink)
        molB.residue_links.append(reslink)

    BUILD.Build_Bonded_Force(molA)
    BUILD.Build_Bonded_Force(molB)
 


    return molA, molB


def Merge_Force_Field(molA, molB, default_lambda, specific_lambda=None, intra_FEP=False):
    if specific_lambda is None:
        specific_lambda = {}
    BUILD.Build_Bonded_Force(molA)
    BUILD.Build_Bonded_Force(molB)

    # ����NB14��ȫ��NB14_extra
    # �������ظ�������Ϊ��ԭλ�����ģ�molAB�����һ�κ�Ͷ�û��������nb14��
    NB14_EXTRA.NB14_To_NB14EXTRA(molA)
    NB14_EXTRA.NB14_To_NB14EXTRA(molB)

    assert len(molA.atoms) == len(molB.atoms)
    molA2molB = {}
    for i, atomA in enumerate(molA.atoms):
        molA2molB[atomA] = molB.atoms[i]

    molR = molA.deepcopy()

    molR2molA = {}
    molR2molB = {}
    for i, atomRet in enumerate(molR.atoms):
        molR2molA[atomRet] = molA.atoms[i]
        molR2molB[atomRet] = molB.atoms[i]
    molB2molR = {value: key for key, value in molR2molB.items()}
    molA2molR = {value: key for key, value in molR2molA.items()}

    charge_lambda = specific_lambda.get("charge", default_lambda)
    mass_lambda = specific_lambda.get("mass", default_lambda)

    for i, atom in enumerate(molR.atoms):
        atom.charge = molR2molA[atom].charge * (1 - charge_lambda) + molR2molB[atom].charge * charge_lambda 
        atom.mass = molR2molA[atom].mass * (1 - mass_lambda) + molR2molB[atom].mass * mass_lambda 
        atom.LJtypeB = molR2molB[atom].LJtype

    for forcename, Rforces in molR.bonded_forces.items():
        if forcename in FEP_Bonded_Force_Merge_Rule.keys():
            temp_lambda = specific_lambda.get(FEP_Bonded_Force_Merge_Rule[forcename]["lambda_name"], default_lambda)
            Bforces = molB.bonded_forces.get(forcename, [])
            FEP_Bonded_Force_Merge_Rule[forcename]["merge_function"](molR, molA, molB,
                                                                     GlobalSetting.BondedForcesMap[forcename],
                                                                     Rforces, Bforces, temp_lambda, molR2molA,
                                                                     molA2molR, molR2molB, molB2molR)
        elif len(Rforces) > 0:
            raise NotImplementedError(forcename + " is not supported for FEP to merge force field yet.")

    for forcename, parameters in FEP_Bonded_Force_Merge_Rule.items():
        temp_lambda = specific_lambda.get(parameters["lambda_name"], default_lambda)

    return molR
