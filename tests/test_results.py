from libtestset.results import ReactionError, AtomizationEnergyError
from libtestset.results import DeviationError
from os import chdir, getcwd
from os.path import exists, join
from pathlib import Path
from shutil import rmtree

import libtestset.results as results
import unittest


class DummyDriver(object):

    def __init__(self, energy):
        self.energy = energy


class TestResults(unittest.TestCase):

    def setUp(self):
        self.base_dir = getcwd()
        self.exec_dir = "testing_dir"
        self.input_dir = "input_files/results/"
        Path(self.exec_dir).mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        chdir(self.base_dir)
        rmtree(self.exec_dir)

    def test_atomization_energy(self):
        ch4_xyz = join(self.input_dir, "samples/ch4.xyz")
        c2h6_xyz = join(self.input_dir, "samples/c2h6.xyz")
        ch4_energy = -2028.4592464933714
        c2h6_energy = -3580.085689008599
        ch4_ref = 419.7
        c2h6_ref = 711.4
        ch4 = results.AtomizationEnergy(ch4_xyz, ch4_energy, ch4_ref)
        c2h6 = results.AtomizationEnergy(c2h6_xyz, c2h6_energy, c2h6_ref)
        self.assertAlmostEqual(ch4.eat, 420.3920952628314, 3)
        self.assertAlmostEqual(c2h6.eat, 714.9320171104587, 3)
        self.assertAlmostEqual(ch4.deviation, -0.6920952628314012, 3)
        self.assertAlmostEqual(c2h6.deviation, -3.532017110458696, 3)
        with self.assertRaises(AtomizationEnergyError):
            results.AtomizationEnergy("wom.bat", 0.0, 1.0)

    def test_reaction_energy(self):
        systems = {"c2h6": DummyDriver(-3580.085689008599),
                   "ch4": DummyDriver(-2028.4592464933714),
                   "h2o": DummyDriver(-2555.0106334067514),
                   "c2h5oh": DummyDriver(-5648.755518949199),
                   "h2": DummyDriver(-420.8596485791793)}
        reaction = "c2h6 + h2o -> c2h5oh + h2"
        reac = results.Reaction(reaction, systems, -24.300)
        self.assertAlmostEqual(reac.energy, -25.481, 3)
        dev1 = reac.deviation
        reaction = "2 c2h6 + 2 h2o -> 2 c2h5oh + 2 h2"
        reac = results.Reaction(reaction, systems, -48.600)
        self.assertAlmostEqual(reac.energy / 2, -25.481, 3)
        self.assertAlmostEqual(dev1, reac.deviation/2, 3)
        reaction = "2 h2o -> 2 h2 + o2"
        with self.assertRaises(ReactionError):
            results.Reaction(reaction, systems, 123.4)
        reaction = "2 h2o -> 2 h2 o2"
        with self.assertRaises(ReactionError):
            results.Reaction(reaction, systems, 432.1)

    def test_deviation(self):
        reac_list = []
        systems = {"c2h6": DummyDriver(-3580.085689008599),
                   "ch4": DummyDriver(-2028.4592464933714),
                   "h2o": DummyDriver(-2555.0106334067514),
                   "c2h5oh": DummyDriver(-5648.755518949199),
                   "h2": DummyDriver(-420.8596485791793)}
        reaction = "c2h6 + h2o -> c2h5oh + h2"
        reac_list.append(results.Reaction(reaction, systems, -24.300))
        reaction = "2 c2h6 + 2 h2o -> 2 c2h5oh + 2 h2"
        reac_list.append(results.Reaction(reaction, systems, -48.600))
        dev = results.Deviations(reac_list)
        self.assertAlmostEqual(dev.mad, abs(dev.msd), 3)
        self.assertAlmostEqual(dev.mad, 1.772, 3)
        self.assertAlmostEqual(dev.rmsd, 1.868, 3)
        self.assertAlmostEqual(dev.max, 2.363, 1)
        reac_list.append("wombat")
        with self.assertRaises(DeviationError):
            results.Deviations(reac_list)

    def test_write_results(self):
        exec_dir = join(self.exec_dir, "write_results")
        Path(exec_dir).mkdir(parents=True, exist_ok=True)
        chdir(exec_dir)
        testsets = {
            "Sample Energies": {
                "path": join("../../", self.input_dir, "samples"),
                "type": "atomization",
                "references": {
                    "ch4": 419.7,
                    "c2h6": 711.4
                }},
            "Sample Reactions": {
                "path": join("../../", self.input_dir, "samples"),
                "type": "reaction",
                "reactions": [
                    {"equation": "c2h6 + h2o -> c2h5oh + h2",
                     "reference": -24.300},
                    {"equation": "2 c2h6 + 2 h2o -> 2 c2h5oh + 2 h2",
                     "reference": -48.600}
                ]
            }
        }
        systems = {"Sample Energies": {"c2h6": DummyDriver(-3580.085689008599),
                                       "ch4": DummyDriver(-2028.4592464933714)},
                   "Sample Reactions": {"c2h6": DummyDriver(712.5),
                                        "h2o": DummyDriver(232.4),
                                        "c2h5oh": DummyDriver(809.0),
                                        "h2": DummyDriver(109.8),
                                        "c4h8": DummyDriver(1155.2)}
                   }
        results.write_results(testsets, systems)
        assert(exists("DFTB_deviations.csv"))
        assert(exists("Sample Energies.csv"))
        assert(exists("Sample Reactions.csv"))
        chdir(self.base_dir)

    def test_write_results_with_dtnn(self):
        exec_dir = join(self.exec_dir, "write_results_with_dtnn")
        Path(exec_dir).mkdir(parents=True, exist_ok=True)
        chdir(exec_dir)
        testsets = {
            "Sample Energies": {
                "path": join("../../", self.input_dir, "samples"),
                "type": "atomization",
                "references": {
                    "ch4": 419.7,
                    "c2h6": 711.4
                }},
            "Sample Reactions": {
                "path": join("../../", self.input_dir, "samples"),
                "type": "reaction",
                "reactions": [
                    {"equation": "c2h6 + h2o -> c2h5oh + h2",
                     "reference": -24.300},
                    {"equation": "2 c2h6 + 2 h2o -> 2 c2h5oh + 2 h2",
                     "reference": -48.600}
                ]
            }
        }
        dftb_systems = {"Sample Energies": {"c2h6": DummyDriver(-3580.085689),
                                            "ch4": DummyDriver(-2028.459246)},
                        "Sample Reactions": {"c2h6": DummyDriver(712.5),
                                             "h2o": DummyDriver(232.4),
                                             "c2h5oh": DummyDriver(809.0),
                                             "h2": DummyDriver(109.8),
                                             "c4h8": DummyDriver(1155.2)}}
        dtnn_systems = {"Sample Energies": {"c2h6": DummyDriver(-3613.830),
                                            "ch4": DummyDriver(-2040.619)},
                        "Sample Reactions": {"c2h6": DummyDriver(712.5),
                                             "h2o": DummyDriver(232.4),
                                             "c2h5oh": DummyDriver(809.0),
                                             "h2": DummyDriver(109.8),
                                             "c4h8": DummyDriver(1155.2)}}
        results.write_results(testsets, dftb_systems, dtnn_systems)
        assert(exists("DFTB_deviations.csv"))
        assert(exists("Sample Energies.csv"))
        assert(exists("Sample Reactions.csv"))
        chdir(self.base_dir)


if __name__ == "__main__":
    unittest.main()
