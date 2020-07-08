from libtestset.dftbplus_runner import DFTBPlusRunnerError, XYZError
from os import chdir, getcwd, mkdir
from os.path import join
from shutil import copy2, rmtree

import libtestset.dftbplus_runner as dftb
import unittest


class TestInputParser(unittest.TestCase):

    def setUp(self):
        self.base_dir = getcwd()
        self.exec_dir = "testing_dir"
        self.input_dir = "input_files/dftbplus_runner/"
        mkdir(self.exec_dir)

    def tearDown(self):
        chdir(self.base_dir)
        rmtree(self.exec_dir)

    def test_run_dftb(self):
        exec_dir = join(self.exec_dir, "run_dftb")
        mkdir(exec_dir)
        ch4_xyz = join(self.input_dir, "testset/ch4.xyz")
        hsd = join(self.input_dir, "dftb_in.hsd")
        driver = dftb.DFTBPlusDriver("/home/mkubillus/bin/dftb+", hsd, ch4_xyz, exec_dir)
        driver.run()
        self.assertAlmostEqual(driver.energy, -2028.4592464933714, 3)

    def test_run_dftb_sccerror(self):
        exec_dir = join(self.exec_dir, "run_dftb_sccerror")
        mkdir(exec_dir)
        ch4_xyz = join(self.input_dir, "scc_error.xyz")
        hsd = join(self.input_dir, "dftb_in.hsd")
        driver = dftb.DFTBPlusDriver("/home/mkubillus/bin/dftb+", hsd, ch4_xyz, exec_dir)
        with self.assertRaises(DFTBPlusRunnerError):
            driver.run()

    def test_run_dftb_geomerror(self):
        exec_dir = join(self.exec_dir, "run_dftb_geomerror")
        mkdir(exec_dir)
        ch4_xyz = join(self.input_dir, "testset/ch4.xyz")
        hsd = join(self.input_dir, "geometry_error.hsd")
        driver = dftb.DFTBPlusDriver("/home/mkubillus/bin/dftb+", hsd, ch4_xyz, exec_dir)
        with self.assertRaises(DFTBPlusRunnerError):
            driver.run()

    def test_xyz2gen(self):
        exec_dir = join(self.exec_dir, "xyz2gen")
        mkdir(exec_dir)
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
