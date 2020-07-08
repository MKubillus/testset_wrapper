from libtestset.dftbplus_runner import DFTBPlusRunnerError, XYZError
from os import chdir, getcwd
from os.path import join
from pathlib import Path
from shutil import rmtree

import libtestset.dftbplus_runner as dftb
import unittest


class TestInputParser(unittest.TestCase):

    def setUp(self):
        self.base_dir = getcwd()
        self.exec_dir = "testing_dir"
        self.input_dir = "input_files/dftbplus_runner/"
        Path(self.exec_dir).mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        chdir(self.base_dir)
        rmtree(self.exec_dir)

    def test_run_testset(self):
        exec_dir = join(self.exec_dir, "run_testset")
        Path(exec_dir).mkdir(parents=True, exist_ok=True)
        chdir(exec_dir)
        hsd = join("../../", self.input_dir, "dftb_in.hsd")
        set_definition = {"path": "../../input_files/dftbplus_runner/testset",
                          "type": "atomization"}
        systems = dftb.run_testset(set_definition, hsd, "dftb+")
        self.assertAlmostEqual(systems["ch4"], -2028.4592464933714, 3)
        self.assertAlmostEqual(systems["c2h6"], -3580.085689008599, 3)
        chdir(self.base_dir)

    def test_run_dftb(self):
        exec_dir = join(self.exec_dir, "run_dftb")
        Path(exec_dir).mkdir(parents=True, exist_ok=True)
        ch4_xyz = join(self.input_dir, "testset/ch4.xyz")
        hsd = join(self.input_dir, "dftb_in.hsd")
        driver = dftb.DFTBPlusDriver("dftb+", hsd, ch4_xyz, exec_dir)
        driver.run()
        self.assertAlmostEqual(driver.energy, -2028.4592464933714, 3)

    def test_run_dftb_sccerror(self):
        exec_dir = join(self.exec_dir, "run_dftb_sccerror")
        Path(exec_dir).mkdir(parents=True, exist_ok=True)
        ch4_xyz = join(self.input_dir, "scc_error.xyz")
        hsd = join(self.input_dir, "dftb_in.hsd")
        driver = dftb.DFTBPlusDriver("dftb+", hsd, ch4_xyz, exec_dir)
        with self.assertRaises(DFTBPlusRunnerError):
            driver.run()

    def test_run_dftb_geomerror(self):
        exec_dir = join(self.exec_dir, "run_dftb_geomerror")
        Path(exec_dir).mkdir(parents=True, exist_ok=True)
        ch4_xyz = join(self.input_dir, "testset/ch4.xyz")
        hsd = join(self.input_dir, "geometry_error.hsd")
        driver = dftb.DFTBPlusDriver("dftb+", hsd, ch4_xyz, exec_dir)
        with self.assertRaises(DFTBPlusRunnerError):
            driver.run()

    def test_xyz2gen(self):
        exec_dir = join(self.exec_dir, "xyz2gen")
        Path(exec_dir).mkdir(parents=True, exist_ok=True)
        with open(join(self.input_dir, "h2.gen"), "r") as f:
            tmp = f.readlines()
        ref = [x for x in tmp if x]
        h2_xyz = join(self.input_dir, "h2.xyz")
        h2_fail = join(self.input_dir, "h2_fail.xyz")
        dftb.DFTBPlusDriver.xyz2gen(h2_xyz, join(exec_dir, "h2.gen"))
        with open(join(exec_dir, "h2.gen"), "r") as f:
            tmp = f.readlines()
        gen = [x for x in tmp if x]
        assert(len(ref) == len(gen))
        i = 0
        while i < len(ref):
            assert(ref[i] == gen[i])
            i += 1
        with self.assertRaises(XYZError):
            dftb.DFTBPlusDriver.xyz2gen(h2_fail, join(exec_dir, "h2.gen"))


if __name__ == "__main__":
    unittest.main()
