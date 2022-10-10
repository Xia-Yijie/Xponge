"""
This **module** is used to calculate the Gasteiger charge.
"""
from rdkit.Chem import rdPartialCharges
from ..helper import Xprint, set_global_alternative_names
from ..helper.rdkit import assign_to_rdmol


def gasteiger(assign):
    """
    This **function** is used to calculate the Gasteiger charge of an assignment

    :param assign: the Assign instance
    :return: a list of charges
    """
    rdmol = assign_to_rdmol(assign)
    rdPartialCharges.ComputeGasteigerCharges(rdmol)
    return [float(atom.GetProp("_GasteigerCharge")) for atom in rdmol.GetAtoms()]


set_global_alternative_names()

Xprint("""Reference for Gasteiger charge
  Gasteiger, J., Marsili, M.
    A new model for calculating atomic charges in molecules
    Tetrahedron Letters, 1978 19(34) 3181-3184
    DOI: 10.1016/S0040-4039(01)94977-9""")
