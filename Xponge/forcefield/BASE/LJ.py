"""
This **module** is the basic setting for the force field property of charge
"""
import numpy as np
from ... import Generate_New_Pairwise_Force_Type
from ...helper import Molecule, AtomType, GlobalSetting, set_dict_value_alternative_name

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
    if self.A is not None and self.B is not None:
        if self.B == 0 or self.A == 0:
            self.sigma = 0
            self.epsilon = 0
        else:
            self.sigma = (self.A / self.B) ** (1 / 6)
            self.epsilon = 0.25 * self.B * self.sigma ** (-6)
        self.A = None
        self.B = None
    if self.sigma is not None:
        self.rmin = self.sigma * (4 ** (1 / 12) / 2)
        self.sigma = None


def lorentz_berthelot_for_a(epsilon1, rmin1, epsilon2, rmin2):
    """
This **function** is used to calculate the A coefficient for Lorentz_Berthelot mix rule
    :param epsilon1:
    :param rmin1:
    :param epsilon2:
    :param rmin2:
    :return:
    """
    return np.sqrt(epsilon1 * epsilon2) * ((rmin1 + rmin2) ** 12)


set_dict_value_alternative_name(globals(), lorentz_berthelot_for_a)


def lorentz_berthelot_for_b(epsilon1, rmin1, epsilon2, rmin2):
    """
This **function** is used to calculate the A coefficient for Lorentz_Berthelot mix rule
    :param epsilon1:
    :param rmin1:
    :param epsilon2:
    :param rmin2:
    :return:
    """
    return np.sqrt(epsilon1 * epsilon2) * 2 * ((rmin1 + rmin2) ** 6)


set_dict_value_alternative_name(globals(), lorentz_berthelot_for_b)


def _find_ab_lj(ljtypes, stat=True):
    """

    :param ljtypes:
    :param stat:
    :return:
    """
    coefficients_a = []
    coefficients_b = []

    for i in range(len(ljtypes)):
        lj_i = LJType.types[ljtypes[i] + "-" + ljtypes[i]]
        if stat:
            j_max = len(ljtypes)
        else:
            j_max = i + 1
        for j in range(j_max):
            lj_j = LJType.types[ljtypes[j] + "-" + ljtypes[j]]
            finded = False
            findnames = [ljtypes[i] + "-" + ljtypes[j], ljtypes[j] + "-" + ljtypes[i]]
            for findname in findnames:
                if findname in LJType.types.keys():
                    finded = True
                    lj_ij = LJType.types[findname]
                    coefficients_a.append(
                        LJType.combining_method_A(lj_ij.epsilon, lj_ij.rmin, lj_ij.epsilon, lj_ij.rmin))
                    coefficients_b.append(
                        LJType.combining_method_B(lj_ij.epsilon, lj_ij.rmin, lj_ij.epsilon, lj_ij.rmin))
                    break
            if not finded:
                coefficients_a.append(LJType.combining_method_A(lj_i.epsilon, lj_i.rmin, lj_j.epsilon, lj_j.rmin))
                coefficients_b.append(LJType.combining_method_B(lj_i.epsilon, lj_i.rmin, lj_j.epsilon, lj_j.rmin))
    return coefficients_a, coefficients_b


def _get_checks(ljtypes, coefficients_a, coefficients_b):
    """

    :param ljtypes:
    :param coefficients_a:
    :param coefficients_b:
    :return:
    """
    checks = {}
    count = 0
    for i in range(len(ljtypes)):
        check_string_a = ""
        check_string_b = ""
        for _ in range(len(ljtypes)):
            check_string_a += "%16.7e" % coefficients_a[count] + " "
            check_string_b += "%16.7e" % coefficients_b[count] + " "
            count += 1

        checks[i] = check_string_a + check_string_b
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


def write_lj(self):
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
        for _ in range(i + 1):
            towrite += "%16.7e" % real_as[count] + " "
            count += 1
        towrite += "\n"
    towrite += "\n"

    count = 0
    for i in range(len(real_ljtypes)):
        for _ in range(i + 1):
            towrite += "%16.7e" % real_bs[count] + " "
            count += 1
        towrite += "\n"
    towrite += "\n"
    towrite += "\n".join(["%d" % (same_type[ljtypemap[atom.LJtype]]) for atom in self.atoms])
    return towrite


Molecule.Set_Save_SPONGE_Input("LJ")(write_lj)
