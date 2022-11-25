"""
This **module** helps to assign bond orders
"""
from itertools import product
from . import AssignRule, Xdict, OrderedDict, set_attribute_alternative_names, deepcopy, np

bo = AssignRule("bo", pure_string=True)

@bo.add_rule("X")
def _(i, assign):
    return assign.atoms[i] in ("H", "F", "Cl", "Br", "I")


@bo.add_rule("Ccn")
def _(i, assign):
    return assign.Atom_Judge(i, "C1") and [assign.atoms[key] for key in assign.bonds[i].keys()][0] == "N"


@bo.add_rule("Cx1")
def _(i, assign):
    return assign.Atom_Judge(i, "C1")


@bo.add_rule("Co2")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "C3") and bonded_os == 2


@bo.add_rule("C")
def _(i, assign):
    return assign.atoms[i] == "C"


@bo.add_rule("Nnn1")
def _(i, assign):
    return assign.Atom_Judge(i, "N1") and [assign.atoms[key] for key in assign.bonds[i].keys()][0] == "N"


@bo.add_rule("Nx1")
def _(i, assign):
    return assign.Atom_Judge(i, "N1")


@bo.add_rule("Nnn2")
def _(i, assign):
    return assign.Atom_Judge(i, "N2") and "N" in {assign.atoms[key] for key in assign.bonds[i].keys()}


@bo.add_rule("Nx2")
def _(i, assign):
    return assign.Atom_Judge(i, "N2")


@bo.add_rule("No2")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "N3") and bonded_os == 2


@bo.add_rule("No1")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "N3") and bonded_os == 1


@bo.add_rule("Nx3")
def _(i, assign):
    return assign.Atom_Judge(i, "N3")


@bo.add_rule("Nx4")
def _(i, assign):
    return assign.Atom_Judge(i, "N4")


@bo.add_rule("On1")
def _(i, assign):
    return assign.Atom_Judge(i, "O1") and [assign.atoms[key] for key in assign.bonds[i].keys()].count("N") == 1


@bo.add_rule("Ox1")
def _(i, assign):
    return assign.Atom_Judge(i, "O1")


@bo.add_rule("Ox2")
def _(i, assign):
    return assign.Atom_Judge(i, "O2")


@bo.add_rule("Px1")
def _(i, assign):
    return assign.Atom_Judge(i, "P1")


@bo.add_rule("Px2")
def _(i, assign):
    return assign.Atom_Judge(i, "P2")


@bo.add_rule("Px3")
def _(i, assign):
    return assign.Atom_Judge(i, "P3")


@bo.add_rule("Po3")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "P4") and bonded_os == 3



@bo.add_rule("Po2")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "P4") and bonded_os == 2


@bo.add_rule("Px4")
def _(i, assign):
    return assign.Atom_Judge(i, "P4")


@bo.add_rule("Sn1")
def _(i, assign):
    return assign.Atom_Judge(i, "S1") and [assign.atoms[key] for key in assign.bonds[i].keys()].count("N") == 1


@bo.add_rule("Sx1")
def _(i, assign):
    return assign.Atom_Judge(i, "S1")


@bo.add_rule("Sx2")
def _(i, assign):
    return assign.Atom_Judge(i, "S2")


@bo.add_rule("S3")
def _(i, assign):
    return assign.Atom_Judge(i, "S3")


@bo.add_rule("So3")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "S4") and bonded_os >= 3


@bo.add_rule("So2")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "S4") and bonded_os == 2


@bo.add_rule("Sx4")
def _(i, assign):
    return assign.Atom_Judge(i, "S4")


class BondOrderAssignment:
    """
    This **class** includes the functions to assign bond orders

    :param original_penalties: the original penalties dict
    :param max_stat: the max valence stats to iterate
    :param assign: the father Assignment instance
    """
    atomic_valence = Xdict({
        "X": OrderedDict({0: 64, 1: 0, 2: 64}),
        "Cn1": OrderedDict({3: 0, 4: 1, 5: 32}),
        "Cx1": OrderedDict({3: 1, 4: 0, 5: 32}),
        "Co2": OrderedDict({4: 32, 5: 0, 6: 32}),
        "C": OrderedDict({2: 64, 3: 32, 4: 0, 5: 32, 6: 64}),
        "Nnn1": OrderedDict({2: 0, 3: 0}),
        "Nx1": OrderedDict({2: 3, 3: 0, 4: 32}),
        "Nnn2": OrderedDict({3: 1, 4: 0}),
        "Nx2": OrderedDict({2: 4, 3: 0, 4: 32}),
        "No2": OrderedDict({3: 64, 4: 32, 5: 0, 6: 32}),
        "No1": OrderedDict({3: 1, 4: 0}),
        "Nx3": OrderedDict({2: 32, 3: 0, 4: 1, 5: 2}),
        "Nx4": OrderedDict({3: 64, 4: 0, 5: 64}),
        "On1": OrderedDict({1: 0, 0: 1}),
        "Ox1": OrderedDict({2: 0, 1: 1, 3: 64}),
        "Ox2": OrderedDict({1: 32, 2: 0, 3:64}),
        "Px1": OrderedDict({2: 2, 3: 0, 4: 32}),
        "Px2": OrderedDict({2: 4, 3: 0, 4: 2}),
        "Px3": OrderedDict({3: 32, 4: 0, 5: 1, 6: 2}),
        "Po2": OrderedDict({5: 32, 6: 0, 7: 32}),
        "Po3": OrderedDict({6: 32, 7: 0}),
        "Px4": OrderedDict({3: 64, 4: 1, 5: 0, 6: 32}),
        "Sn1": OrderedDict({1: 0, 2: 1}),
        "Sx1": OrderedDict({1: 2, 2: 0, 3: 64}),
        "Sx2": OrderedDict({1: 32, 2: 0, 3: 32}),
        "Sx3": OrderedDict({3: 1, 4: 0, 5: 2, 6: 2}),
        "So2": OrderedDict({6: 0, 7: 32}),
        "So3": OrderedDict({6: 32, 7: 0}),
        "Sx4": OrderedDict({4: 4, 5: 2, 6:0}),
    })
    def __init__(self, original_penalties, max_step, max_stat, assign, debug=False):
        self.debug = debug
        self.max_step = max_step
        self.max_stat = max_stat
        self.assign = assign
        self.original_penalties = original_penalties
        self.uc = [set(filter(lambda aj: self.assign.bonds[ai][aj] == -1, self.assign.bonds[ai]))
                     for ai in range(self.assign.atom_numbers)]
        self.connected = [sum(filter(lambda x: x > 0, self.assign.bonds[ai].values()))
                     for ai in range(self.assign.atom_numbers)]
        self.valence_best = [next(iter(pi)) - self.connected[i] for i, pi in enumerate(self.original_penalties)]
        self.valence = deepcopy(self.valence_best)
        self.penalties = Xdict(not_found_method=lambda x: [])
        self._get_penalties(original_penalties)
        self.stat_position = 0
        self.current_stat = 0
        self.points = []
        self.cached = {}
        set_attribute_alternative_names(self)

    def _get_penalties(self, original_penalties):
        """

        :param original_penalties:
        :return:
        """
        for atom, pi in enumerate(original_penalties):
            for valence, penalty in pi.items():
                self.penalties[penalty].append((atom, penalty, valence))

    def _preprocess_penalties(self, n):
        if n in self.cached:
            return self.cached[n]
        if n == 1:
            toret = [[point] for point in self.penalties[n]]
        else:
            toret =  []
            have_added = set()
            for i in range(1, n // 2 + 1):
                r1 = self.cached[i]
                r2 = self.cached[n - i]
                for ri, rj in product(r1, r2):
                    if {r0[0] for r0 in ri} & {r0[0] for r0 in rj}:
                        continue
                    rij = ri + rj
                    rij.sort()
                    rij_checkstring = "+".join([f"{r0[0]}-{r0[1]}-{r0[2]}" for r0 in rij])
                    if rij_checkstring not in have_added:
                        have_added.add(rij_checkstring)
                        toret.append(rij)
            toret.extend([[point] for point in self.penalties.get(n, [])])
        self.cached[n] = toret
        return toret

    def _get_next_valence(self):
        if self.current_stat < self.max_stat:
            if self.stat_position < len(self.points):
                self.valence = deepcopy(self.valence_best)
                has_negative_value = False
                for point in self.points[self.stat_position]:
                    self.valence[point[0]] = point[2] - self.connected[point[0]]
                    if self.valence[point[0]] < 0:
                        has_negative_value = True
                        break
                self.stat_position += 1
                if has_negative_value:
                    return True
            else:
                self.current_stat += 1
                self.stat_position = 0
                self.points = self._preprocess_penalties(self.current_stat)
                if self.debug:
                    print("stat=", self.current_stat)
                    print("points=\n", self.points)
                return True


    def main(self):
        """
        This **function** is the main function to do the bond order assignment

        :return: True for success, False for failure
        """
        count = 0
        success = False
        while count < self.max_step and self.current_stat < self.max_stat and not success:
            bonds = deepcopy(self.assign.bonds)
            uc = deepcopy(self.uc)
            valence = deepcopy(self.valence)
            guess_bonds = []
            determined = False
            while not determined:
                index_sort = np.argsort([len(uci) for uci in uc]).tolist()
                no_basic_rule = True
                for i in index_sort:
                    if self.debug:
                        print(i, valence, uc)
                    if len(uc[i]) == valence[i] and valence[i] != 0:
                        while uc[i]:
                            j = uc[i].pop()
                            bonds[i][j] = 1
                            bonds[j][i] = 1
                            uc[j].remove(i)
                            valence[j] -= 1
                        valence[i] = 0
                        no_basic_rule = False
                    elif len(uc[i]) == 1:
                        j = uc[i].pop()
                        uc[j].remove(i)
                        bonds[j][i] = valence[i]
                        bonds[i][j] = valence[i]
                        valence[j] -= valence[i]
                        valence[i] -= valence[i]
                        no_basic_rule = False
                success = True
                for i in range(self.assign.atom_numbers):
                    if valence[i] != 0 or len(uc[i]) != 0:
                        success = False
                    if len(uc[i]) == 0 and valence[i] != 0 or valence[i] == 0 and len(uc[i]) != 0:
                        success = False
                        determined = True
                        break
                if success:
                    determined = True
                if not determined and no_basic_rule:
                    if not guess_bonds:
                        i = index_sort.pop()
                        j = uc[i].pop()
                        uc[j].remove(i)
                        guess_bonds = [i, j]
                        bonds[i][j] = 1
                        bonds[j][i] = 1
                        valence[i] -= 1
                        valence[j] -= 1
                    else:
                        i, j = guess_bonds
                        if bonds[i][j] == 3:
                            determined = True
                            success = False
                            break
                        else:
                            bonds[i][j] += 1
                            bonds[j][i] += 1
                            valence[i] -= 1
                            valence[j] -= 1
            if self.debug:
                print(valence, uc)
            if not success:
                count += 1
                while self._get_next_valence():
                    pass
                if self.debug:
                    print("-"*20, self.points[self.stat_position - 1])
        if success:
            self.assign.bonds = bonds
        return success

for key, value in BondOrderAssignment.atomic_valence.items():
    BondOrderAssignment.atomic_valence[key] = OrderedDict(sorted(value.items(), key=lambda t: t[1]))
