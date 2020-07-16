from os.path import exists, join, basename
from libtestset.constants import atomic_energies, Hydrogen


import csv
import numpy as np


class ReactionError(Exception):
    pass


class AtomizationEnergyError(Exception):
    pass


class DeviationError(Exception):
    pass


class Reaction(object):
    """Calculates reaction energies from testset data.

    Parses the reaction equation string and uses the systems dictionary to get
    to calculate the reaction energy.
    """

    def __init__(self, reaction, systems, reference):
        """Instantiation method.

        Inputs:
        @:param reaction: Reaction string as in the input file.
        @:param systems: Testset systems dictionary with DFTB energies.
        @:param reference: reaction's reference energy.
        """
        self.reaction = reaction
        self.ref = reference
        self._parse_reaction()
        self._energy = self._calculate_energy(systems)

    @property
    def energy(self):
        """Returns reaction energy in kcal/mol."""
        return self._energy

    @property
    def deviation(self):
        """Returns deviation from reference energy in kcal/mol."""
        return self.ref - self.energy

    def _parse_reaction(self):
        """Parses a reaction string.

        @:returns dict with system name and stochiometric factor pairs.
        """
        found_arrow = False
        reac_split = self.reaction.split()
        reac_dict = dict()
        amount = None
        system = None
        i = 0
        while i < len(reac_split):
            if reac_split[i] == "+" or reac_split[i] == "->":
                # Add entry to react_dict with correct stochiometric sign
                if found_arrow:
                    reac_dict[system] = amount
                else:
                    reac_dict[system] = -amount
                # Reset factor and molecule pair
                amount = None
                system = None
                if reac_split[i] == "->":
                    found_arrow = True
            elif not amount:  # found unsigned stochiometric factor
                try:
                    amount = int(reac_split[i])
                except ValueError:
                    amount = 1
                    system = reac_split[i]
            elif not system:  # found molecule name
                system = reac_split[i]
            else:  # error catcher
                msg = "Invalid reaction entry '%s' of reaction '%s' found!"
                raise ReactionError(msg % (reac_split[i], self.reaction))
            i += 1
        reac_dict[system] = amount  # add last entry
        self.reaction_dictionary = reac_dict.copy()

    def _calculate_energy(self, systems):
        energy = 0.0
        for sys, factor in self.reaction_dictionary.items():
            try:
                sys_energy = systems[sys]
            except KeyError:
                msg = "Energy for system '%s' of reaction '%s' not found!"
                raise ReactionError(msg % (sys, self.reaction))
            if sys == "h2":  # Corrected hydrogen molecule
                sys_energy = Hydrogen.total_energy
            energy += abs(sys_energy) * factor
        return energy


class AtomizationEnergy(object):
    """Calculates atomization energies from testset data."""

    def __init__(self, xyz, total_energy, reference):
        """Instantiation of AtomizationEnergy object.

        Inputs:
        @:param xyz: path to xyz-file of system.
        @:param total_energy: system's total electronic energy in kcal/mol.
        @:param reference: system's reference energy from testset.
        """
        self.atomics = atomic_energies
        self.total_energy = total_energy
        self.xyz = xyz
        self.ref = reference
        self._eat = self._calculate_eat()

    @property
    def eat(self):
        """Returns the atomization energy of given molecule."""
        return self._eat

    @property
    def deviation(self):
        """Returns the deviation from reference in kcal/mol."""
        return self.ref - self.eat

    def _calculate_eat(self):
        """Calculates the atomization energy.

        In case of the hydrogen molecule (H2) it uses the corrected
        DFTB3 energy.
        """
        if basename(self.xyz) == "h2.xyz":
            energy = Hydrogen.eat
        else:
            atoms = self._parse_xyz(self.xyz)
            energy = self.total_energy
            for at, amount in atoms.items():
                energy -= self.atomics[at] * amount
        return abs(energy)

    @staticmethod
    def _parse_xyz(xyz):
        """Opens object's xyz file and counts atoms and occurrences.

        Inputs:
        @:param xyz: path to xyz file.

        @:returns dictionary with atoms and their occurences in the geometry.
        """
        if not exists(xyz):
            msg = "XYZ file at %s not found!"
            from os import getcwd
            print(getcwd())
            raise AtomizationEnergyError(msg % xyz)
        with open(xyz, "r") as geom_file:
            geom = geom_file.readlines()[2:]
        atoms = dict()
        for line in geom:
            splt = line.split()
            if len(splt) < 4:
                continue
            atom = splt[0].lower()
            if atom not in atoms.keys():
                atoms[atom] = 1
            else:
                atoms[atom] += 1
        return atoms


class Deviations(object):
    """Calculates deviations from Reaction or AtomizationEnergy objects.

    Uses lists of Reaction objects or AtomizationEnergy objects of a test set as
    inputs and calculates Mean Signed Deviation (MSD), Mean Absolute Deviation
    (MAD), Root Mean Square Deviation (RMSD), and the maximum absolute deviation
    (MAX) of the test set.
    """

    def __init__(self, obj_list):
        """Instantiation method.

        Inputs:
        @:param input: list of Reaction or AtomizationEnergy objects.
        """
        signed_deviations = []
        abs_deviations = []
        square_sum = 0.0
        for entry in obj_list:
            # sanity check
            if not isinstance(entry, (AtomizationEnergy, Reaction)):
                msg = ("Each entry in the list has to be a AtomizationEnergy "
                       "or Reaction object, instead received: %s")
                raise DeviationError(msg % type(entry))
            # deviation sums
            signed_deviations.append(entry.deviation)
            abs_deviations.append(abs(entry.deviation))
            square_sum += entry.deviation**2
        n_entries = len(signed_deviations)  # number of entries
        # Calculate error values
        self._msd = np.sum(signed_deviations) / n_entries
        self._mad = np.sum(abs_deviations) / n_entries
        self._rmsd = np.sqrt(square_sum / n_entries)
        self._max = np.amax(abs_deviations)

    @property
    def msd(self):
        """Returns Mean Signed Deviation."""
        return self._msd

    @property
    def mad(self):
        """Returns Mean Absolute Deviation."""
        return self._mad

    @property
    def rmsd(self):
        """Returns Root Mean Square Deviation."""
        return self._rmsd

    @property
    def max(self):
        """Returns maximum absolute Deviation."""
        return self._max


def write_results(testsets, dftb_calcs, dtnn_calcs=None):
    dftb_deviations = dict()
    dtnn_deviations = dict()
    for set_name in testsets.keys():
        set_type = testsets[set_name]["type"]
        dftb_systems = dftb_calcs[set_name]
        if dtnn_calcs:
            dtnn_systems = dtnn_calcs[set_name]
        if set_type == "reaction":
            inputs = testsets[set_name]["reactions"]
            dftb_reacs = _get_reactions(dftb_systems, inputs)
            dftb_deviations[set_name] = Deviations(dftb_reacs)
            if dtnn_calcs:
                dtnn_reacs = _get_reactions(dtnn_systems, inputs)
                dtnn_deviations[set_name] = Deviations(dtnn_reacs)
                _write_reactions(set_name, dftb_reacs, dtnn_reacs)
            else:
                _write_reactions(set_name, dftb_reacs)
        elif set_type == "atomization":
            set_path = testsets[set_name]["path"]
            inputs = testsets[set_name]["references"]
            dftb_eats = _get_atomizations(dftb_systems, inputs, set_path)
            dftb_deviations[set_name] = Deviations(dftb_eats)
            if dtnn_calcs:
                dtnn_eats = _get_atomizations(dtnn_systems, inputs, set_path)
                dtnn_deviations[set_name] = Deviations(dtnn_eats)
                _write_atomizations(set_name, dftb_eats, dtnn_eats)
            else:
                _write_atomizations(set_name, dftb_eats)
    _write_deviations(dftb_deviations, "DFTB_deviations.csv")
    if dtnn_calcs:
        _write_deviations(dtnn_deviations, "DTNN_deviations.csv")


def _get_reactions(systems, inputs):
    reac_list = []
    for reaction in inputs:
        eq = reaction["equation"]
        ref = reaction["reference"]
        reac_list.append(Reaction(eq, systems, ref))
    return reac_list


def _get_atomizations(systems, inputs, set_path):
    eat_list = []
    for name, ref in inputs.items():
        xyz = join(set_path, "%s.xyz" % name)
        eat_list.append(AtomizationEnergy(xyz, systems[name], ref))
    return eat_list


def _write_reactions(set_name, dftb_reacs, dtnn_reacs=None):
    if dtnn_reacs:
        _write_reactions_with_dtnn(set_name, dftb_reacs, dtnn_reacs)
    else:
        categories = ["Reaction", "Reference", "DFTB3", "ΔDFTB3"]
        data = []
        for reaction in dftb_reacs:
            data.append({"Reaction": reaction.reaction,
                         "Reference": str(np.round(reaction.ref, 2)),
                         "DFTB3": str(np.round(reaction.energy, 2)),
                         "ΔDFTB3": str(np.round(reaction.deviation, 2))
                         })
        print("Writing results of reaction test set %s "
              "to file %s" % (set_name, "%s.csv" % set_name))
        with open("%s.csv" % set_name, "w") as out:
            writer = csv.DictWriter(out, fieldnames=categories)
            writer.writeheader()
            for entry in data:
                writer.writerow(entry)


def _write_reactions_with_dtnn(set_name, dftb_reacs, dtnn_reacs):
    categories = ["Reaction", "Reference", "DFTB3", "DTNN", "ΔDFTB3", "ΔDTNN"]
    data = []
    for i, dftb_reaction in enumerate(dftb_reacs):
        dtnn_reaction = dtnn_reacs[i]
        data.append({"Reaction": dftb_reaction.reaction,
                     "Reference": str(np.round(dftb_reaction.ref, 2)),
                     "DFTB3": str(np.round(dftb_reaction.energy, 2)),
                     "DTNN": str(np.round(dtnn_reaction.energy, 2)),
                     "ΔDFTB3": str(np.round(dftb_reaction.deviation, 2)),
                     "ΔDTNN": str(np.round(dtnn_reaction.deviation, 2))
                     })
    print("Writing results of reaction test set %s "
          "to file %s" % (set_name, "%s.csv" % set_name))
    with open("%s.csv" % set_name, "w") as out:
        writer = csv.DictWriter(out, fieldnames=categories)
        writer.writeheader()
        for entry in data:
            writer.writerow(entry)


def _write_atomizations(set_name, dftb_eats, dtnn_eats=None):
    if dtnn_eats:
        _write_atomizations_with_dtnn(set_name, dftb_eats, dtnn_eats)
    else:
        categories = ["System", "Reference", "DFTB3", "ΔDFTB3"]
        data = []
        for energy in dftb_eats:
            sys_name = basename(energy.xyz)[:-4]
            data.append({"System": sys_name,
                         "Reference": str(np.round(energy.ref, 2)),
                         "DFTB3": str(np.round(energy.eat, 2)),
                         "ΔDFTB3": str(np.round(energy.deviation, 2))
                         })
        print("Writing results of atomization energy test set %s "
              "to file %s" % (set_name, "%s.csv" % set_name))
        with open("%s.csv" % set_name, "w") as out:
            writer = csv.DictWriter(out, fieldnames=categories)
            writer.writeheader()
            for entry in data:
                writer.writerow(entry)


def _write_atomizations_with_dtnn(set_name, dftb_eats, dtnn_eats=None):
    categories = ["System", "Reference", "DFTB3", "DTNN", "ΔDFTB3", "ΔDTNN"]
    data = []
    for i, dftb_energy in enumerate(dftb_eats):
        dtnn_energy = dtnn_eats[i]
        sys_name = basename(dftb_energy.xyz)[:-4]
        data.append({"System": sys_name,
                     "Reference": str(np.round(dftb_energy.ref, 2)),
                     "DFTB3": str(np.round(dftb_energy.eat, 2)),
                     "DTNN": str(np.round(dtnn_energy.eat, 2)),
                     "ΔDFTB3": str(np.round(dftb_energy.deviation, 2)),
                     "ΔDTNN": str(np.round(dtnn_energy.deviation, 2)),
                     })
    print("Writing results of atomization energy test set %s "
          "to file %s" % (set_name, "%s.csv" % set_name))
    with open("%s.csv" % set_name, "w") as out:
        writer = csv.DictWriter(out, fieldnames=categories)
        writer.writeheader()
        for entry in data:
            writer.writerow(entry)


def _write_deviations(dftb_deviations, filename):
    categories = ["Set Name", "MSD", "MAD", "RMSD", "MAX"]
    data = []
    for testset in dftb_deviations.keys():
        dev = dftb_deviations[testset]
        data.append({"Set Name": testset,
                     "MSD": str(np.round(dev.msd, 2)),
                     "MAD": str(np.round(dev.mad, 2)),
                     "RMSD": str(np.round(dev.rmsd, 2)),
                     "MAX": str(np.round(dev.max, 2))
                     })
    print("Writing deviations for all test sets to file %s" % filename)
    with open(filename, "w") as out:
        writer = csv.DictWriter(out, fieldnames=categories)
        writer.writeheader()
        for entry in data:
            writer.writerow(entry)
