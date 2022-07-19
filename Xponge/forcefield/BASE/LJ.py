"""
This **module** is the basic setting for the force field property of charge
"""
import numpy as np
from ... import Generate_New_Pairwise_Force_Type
from ...helper import Molecule, AtomType, GlobalSetting

AtomType.Add_Property({"LJtype": str})

LJType = Generate_New_Pairwise_Force_Type("LJ",
                                          {"epsilon": float, "rmin": float, "sigma": float, "A": float, "B": float})

LJType.Set_Property_Unit("rmin", "distance", "A")
LJType.Set_Property_Unit("sigma", "distance", "A")
LJType.Set_Property_Unit("epsilon", "energy", "kcal/mol")
LJType.Set_Property_Unit("A", "energy·distance^6", "kcal/mol·A^6")
LJType.Set_Property_Unit("B", "energy·distance^12", "kcal/mol·A^12")


@GlobalSetting.Add_Unit_Transfer_Function(LJType)
def lj_unit_transfer(self):
    """
This **function** is used to transfer the units of lj
    :param self:
    :return:
    """
    if self.A != None and self.B != None:
        if self.B == 0 or self.A == 0:
            self.sigma = 0
            self.epsilon = 0
        else:
            self.sigma = (self.A / self.B) ** (1 / 6)
            self.epsilon = 0.25 * self.B * self.sigma ** (-6)
        self.A = None
        self.B = None
    if self.sigma != None:
        self.rmin = self.sigma * (4 ** (1 / 12) / 2)
        self.sigma = None


def Lorentz_Berthelot_For_A(epsilon1, rmin1, epsilon2, rmin2):
    """
This **function** is used to calculate the A coefficient for Lorentz_Berthelot mix rule
    :param epsilon1:
    :param rmin1:
    :param epsilon2:
    :param rmin2:
    :return:
    """
    return np.sqrt(epsilon1 * epsilon2) * ((rmin1 + rmin2) ** 12)


def Lorents_Berthelot_For_B(epsilon1, rmin1, epsilon2, rmin2):
    """
This **function** is used to calculate the A coefficient for Lorentz_Berthelot mix rule
    :param epsilon1:
    :param rmin1:
    :param epsilon2:
    :param rmin2:
    :return:
    """
    return np.sqrt(epsilon1 * epsilon2) * 2 * ((rmin1 + rmin2) ** 6)


def _find_ab_lj(ljtypes, stat=True):
    """

    :param ljtypes:
    :param stat:
    :return:
    """
    As = []
    Bs = []

    for i in range(len(ljtypes)):
        LJ_i = LJType.types[ljtypes[i] + "-" + ljtypes[i]]
        if stat:
            j_max = len(ljtypes)
        else:
            j_max = i + 1
        for j in range(j_max):
            LJ_j = LJType.types[ljtypes[j] + "-" + ljtypes[j]]
            finded = False
            findnames = [ljtypes[i] + "-" + ljtypes[j], ljtypes[j] + "-" + ljtypes[i]]
            for findname in findnames:
                if findname in LJType.types.keys():
                    finded = True
                    LJ_ij = LJType.types[findname]
                    As.append(LJType.combining_method_A(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin))
                    Bs.append(LJType.combining_method_B(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin))
                    break
            if not finded:
                As.append(LJType.combining_method_A(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin))
                Bs.append(LJType.combining_method_B(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin))
    return As, Bs


def _get_checks(ljtypes, As, Bs):
    """

    :param ljtypes:
    :param As:
    :param Bs:
    :return:
    """
    checks = {}
    count = 0
    for i in range(len(ljtypes)):
        check_string_A = ""
        check_string_B = ""
        for j in range(len(ljtypes)):
            check_string_A += "%16.7e" % As[count] + " "
            check_string_B += "%16.7e" % Bs[count] + " "
            count += 1

        checks[i] = check_string_A + check_string_B
    return checks


def _judge_same_type(ljtypes, checks):
    """

    :param ljtypes:
    :param checks:
    :return:
    """
    same_type = {i: i for i in range(len(ljtypes))}
    for i in range(len(ljtypes) - 1, -1, -1):
        for j in range(i + 1, len(ljtypes)):
            if checks[i] == checks[j]:
                same_type[j] = i
    return same_type


def _get_real_lj(ljtypes, same_type):
    """

    :param ljtypes:
    :param same_type:
    :return:
    """
    real_ljtypes = []
    tosub = 0
    for i in range(len(ljtypes)):
        if same_type[i] == i:
            real_ljtypes.append(ljtypes[i])
            same_type[i] -= tosub
        else:
            same_type[i] = same_type[same_type[i]]
            tosub += 1
    return real_ljtypes


def write_LJ(self):
    """
This **function** is used to write SPONGE input file
    :param self:
    :return:
    """
    ljtypes = []
    ljtypemap = {}
    for atom in self.atoms:
        if atom.LJtype not in ljtypemap.keys():
            ljtypemap[atom.LJtype] = len(ljtypes)
            ljtypes.append(atom.LJtype)

    coefficients_a, coefficients_b = _find_ab_lj(ljtypes)

    checks = _get_checks(ljtypes, coefficients_a, coefficients_b)

    same_type = _judge_same_type(ljtypes, checks)

    real_ljtypes = _get_real_lj(ljtypes, same_type)

    real_as, real_bs = _find_ab_lj(real_ljtypes, False)

    towrite = "%d %d\n\n" % (len(self.atoms), len(real_ljtypes))
    count = 0
    for i in range(len(real_ljtypes)):
        for j in range(i + 1):
            towrite += "%16.7e" % real_as[count] + " "
            count += 1
        towrite += "\n"
    towrite += "\n"

    count = 0
    for i in range(len(real_ljtypes)):
        for j in range(i + 1):
            towrite += "%16.7e" % real_bs[count] + " "
            count += 1
        towrite += "\n"
    towrite += "\n"
    towrite += "\n".join(["%d" % (same_type[ljtypemap[atom.LJtype]]) for atom in self.atoms])
    return towrite


Molecule.Set_Save_SPONGE_Input("LJ")(write_LJ)
