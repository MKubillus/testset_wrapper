from libtestset import input_parser
from libtestset.input_parser import InputError
from os import chdir, getcwd, mkdir
from os.path import join
from shutil import copy2, rmtree

import unittest


class TestInputParser(unittest.TestCase):

    def setUp(self):
        self.base_dir = getcwd()
        self.exec_dir = "testing_dir"
        self.input_dir = "input_files/test_input_parser"
        mkdir(self.exec_dir)

    def tearDown(self):
        chdir(self.base_dir)
        rmtree(self.exec_dir)

    def test_write_template(self):
        exec_dir = join(self.exec_dir, "write_template")
        mkdir(exec_dir)
        chdir(exec_dir)
        input_parser.write_template("tmp.yml")
        chdir(self.base_dir)

    def test_load_errors(self):
        exec_dir = join(self.exec_dir, "load_errors")
        mkdir(exec_dir)
        copy2(join(self.input_dir, "fail1.yml"), exec_dir)
        copy2(join(self.input_dir, "fail2.yml"), exec_dir)
        copy2(join(self.input_dir, "dftbplus.hsd"), exec_dir)
        chdir(exec_dir)
        with self.assertRaises(InputError):
            input_parser.load("fail1.yml")
        with self.assertRaises(InputError):
            input_parser.load("fail2.yml")
        chdir(self.base_dir)


if __name__ == "__main__":
    unittest.main()
