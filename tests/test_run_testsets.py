from os import chdir, getcwd
from os.path import join
from pathlib import Path
from shutil import copytree, rmtree

import run_testsets
import unittest


class TestDFTBPlusRunner(unittest.TestCase):

    def setUp(self):
        self.base_dir = getcwd()
        self.exec_dir = "testing_dir"
        self.input_dir = "input_files/run_testsets/"
        Path(self.exec_dir).mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        chdir(self.base_dir)
        rmtree(self.exec_dir)

    def test_run_testsets(self):
        exec_dir = join(self.exec_dir, "run_testsets")
        copytree(self.input_dir, exec_dir)
        chdir(exec_dir)
        run_testsets.main()
        chdir(self.base_dir)


if __name__ == "__main__":
    unittest.main()
