from libtestset.dtnn_runner import DTNNDriver, run_testset
from os import chdir, getcwd
from os.path import join
from pathlib import Path
from shutil import rmtree

import unittest


class TestDFTBPlusRunner(unittest.TestCase):

    def setUp(self):
        self.base_dir = getcwd()
        self.exec_dir = "testing_dir"
        self.input_dir = "input_files/dtnn_runner/"
        Path(self.exec_dir).mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        chdir(self.base_dir)
        rmtree(self.exec_dir)

    def test_run_testset(self):
        exec_dir = join(self.exec_dir, "run_dtnn_testset")
        Path(exec_dir).mkdir(parents=True, exist_ok=True)
        set_definition = {"path": join(self.input_dir, "testset")}
        model = join(self.input_dir, "dftbnn.dtnn")
        dftbplus = "dftb+"
        skf = "/home/mkubillus/slko/3ob-3-1/"
        systems = run_testset(set_definition, model, dftbplus, skf)
        self.assertAlmostEqual(systems["ch4"], -2040.619, 3)
        self.assertAlmostEqual(systems["c2h6"], -3613.830, 3)

    def test_run_dtnn(self):
        exec_dir = join(self.exec_dir, "run_dtnn")
        Path(exec_dir).mkdir(parents=True, exist_ok=True)
        ch4_xyz = join(self.input_dir, "ch4.xyz")
        model = join(self.input_dir, "dftbnn.dtnn")
        exe = "dftb+"
        skf = "/home/mkubillus/slko/3ob-3-1/"
        driver = DTNNDriver(model, ch4_xyz, exe, skf, exec_dir)
        driver.run()
        self.assertAlmostEqual(driver.energy, -2040.619, 3)


if __name__ == "__main__":
    unittest.main()
