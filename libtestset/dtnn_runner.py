from ase.calculators import dftb, mixing
from ase.io import read
from ase.optimize import BFGS
from io import StringIO
from libtestset.constants import UnitConversion as Units
from os import chdir, environ, getcwd, scandir
from os.path import join
from pathlib import Path
from random import choice
from schnetpack.interfaces.ase_interface import SpkCalculator
from shutil import rmtree
from string import ascii_lowercase

import libtestset.constants as c
import sys
import torch


class DTNNRunnerError(Exception):
    """Error raised when DTNN calculation goes wrong."""

    def __init__(self, exec_dir, msg=None):
        """Describes how to raise the error."""
        if not msg:
            msg = ("DTNN run crashed with an unknown error!\n\n"
                   "Crashed calculation in: %s" % exec_dir)
        else:
            msg += "\nCrashed calculation in: %s" % exec_dir
        super().__init__(msg)


class CaptureSTDOUT(list):
    """Captures all STDOUT and puts it in a list rather than printing.

    When used in a 'with CaptureSTDOUT' statement it will grab all
    STDOUT prints within the statement and add it to a list that can then
    be handles separately.

    Usage example:
    with CaptureSTDOUT as outputs:
        print("Hello world!")
    print("Caught:")
    print(outputs)

    > Caught:
    > ["Hello world!"]
    """

    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # Prevent memory leakage
        sys.stdout = self._stdout


class CaptureSTDERR(list):
    """Captures all STDERR and puts it in a list rather than printing.

    When used in a 'with CaptureSTDERR()' statement it will grab all
    STDOUT prints within the statement and add it to a list that can then
    be handles separately.

    Usage example:
    with CaptureSTDERR() as outputs:
        print("Hello world!")
    print("Caught:")
    print(outputs)

    > Caught:
    > ["Hello world!"]
    """

    def __enter__(self):
        self._stderr = sys.stderr
        sys.stderr = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # Prevent memory leakage
        sys.stderr = self._stderr


class DTNNDriver(object):
    """Runs DTNN calculations.

    Object that runs DTNN calculations with DFTB+ using ASE and SchNet.

    Inputs for instantiation:
    @:param model: path to the DTNN model to use.
    @:param xyz: Path to xyz geometry file.
    @:param exe: DFTB+ executable
    @:param skf: Path to Slater-Koster parameter files.
    @:param exec_dir: Directory to run the calculation in.
    """

    def __init__(self, model, xyz, exe, skf, exec_dir):
        self._model_path = model
        self.xyz = xyz
        self.ase_atoms = read(xyz)
        self.exec_dir = exec_dir
        self.base_dir = getcwd()
        self._energy = None
        self._coords = None
        self._atoms = None
        self._model = None
        self._schnet_calc = None
        environ["DFTB_COMMAND"] = exe
        environ["DFTB_PREFIX"] = skf

    def run(self):
        """Sets up the directory and runs the DTNN calculation."""
        with CaptureSTDERR() as errors:
            self._model = torch.load(self._model_path, map_location="cuda")
            self._schnet_calc = SpkCalculator(self._model, device='cuda',
                                              energy="ErepD3", forces="FOR3")
        Path(self.exec_dir).mkdir(parents=True, exist_ok=True)
        chdir(self.exec_dir)
        ase_dftb = dftb.Dftb(
            label="dftb_calculator",
            atoms=self.ase_atoms,
            run_manyDftb_steps=True,
            Hamiltonian_SCC="Yes",
            Hamiltonian_ThirdOrderFull="Yes",
            Hamiltonian_HubbardDerivs_="",
            Hamiltonian_HubbardDerivs_C=c.hubbard_derivatives["c"],
            Hamiltonian_HubbardDerivs_H=c.hubbard_derivatives["h"],
            Hamiltonian_HubbardDerivs_N=c.hubbard_derivatives["n"],
            Hamiltonian_HubbardDerivs_O=c.hubbard_derivatives["o"],
            Hamiltonian_PolynomialRepulsive_="",
            Hamiltonian_PolynomialRepulsive_setForAll="{Yes}",
            Analysis_="",
            Analysis_CalculateForces="Yes")
        mix = mixing.SumCalculator([ase_dftb, self._schnet_calc])
        self.ase_atoms.set_calculator(mix)
        opt = BFGS(self.ase_atoms, logfile="BFGS_optimization.log")
        # Redirect opt.run() STDOUT and STERR prints to file
        with CaptureSTDOUT() as outputs, CaptureSTDERR(errors) as errors:
            opt.run(fmax=0.00005)
        opt.logfile.close()
        ev2kcal = Units.ev2au * Units.au2kcal
        self._energy = self.ase_atoms.get_total_energy()[0] * ev2kcal
        self._coords = self.ase_atoms.get_positions()
        self._atoms = list(self.ase_atoms.symbols)
        with open("BFGS_output.log", "w") as outfile:
            for line in outputs:
                outfile.write(line)
        with open("BFGS_stderr.log", "w") as errfile:
            for line in errors:
                errfile.write(line)
        chdir(self.base_dir)

    @property
    def energy(self):
        """Returns calculated system energy in kcal/mol."""
        return self._energy

    @property
    def coordinates(self):
        """Returns the optimized geometry vectors.

        @:returns numpy Nx3 ndarray.
        """
        return self._coords

    @property
    def atoms(self):
        """Returns the geometry atom list.

        @:returns list of strings.
        """
        return self._atoms


def run_testset(set_definition, model, dftbplus, skf):
    """Runs all systems in a given testset.

    Inputs:
    @:param set_definition: dictionary describing the testset as given in input
        file.
    @:param model: Path to the DTNN model file.
    @:param executable: Path to DFTB+ executable.
    @:param skf: Path to skf parameter files.

    @:returns Dictionary of systems and their total energies in kcal/mol.
    """
    set_path = join(getcwd(), set_definition["path"])
    systems = dict()
    for entry in scandir(set_path):
        if entry.is_file() and entry.name.endswith(".xyz"):
            sys_name = entry.name[:-4]
            xyz = entry.path
            exec_dir = get_random_folder(prefix="dtnn_run_")
            driver = DTNNDriver(model, xyz, dftbplus, skf, exec_dir)
            driver.run()
            systems[sys_name] = driver
            rmtree(exec_dir)
    return systems


def get_random_folder(prefix="", length=8):
    random_string = "".join(choice(ascii_lowercase) for _i in range(length))
    return prefix + random_string


if __name__ == "__main__":
    pass
