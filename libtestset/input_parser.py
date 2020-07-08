from inspect import currentframe, getfile
from os.path import dirname, exists, join
from shutil import copy2
from sys import exit

import yaml


class InputError(Exception):
    """Default exception for input parsing errors."""
    pass


def load(filename):
    """Main driver function to load the input."""
    if not exists(filename):
        write_template(filename)
        exit()
    with open(filename, "r") as infile:
        settings = yaml.safe_load(infile)
    if not exists(settings["Options"]["DFTBPlusHSD"]):
        msg = "DFTBPlusHSD file at %s does not exist."
        raise InputError(msg % settings["Options"]["DFTBPlusHSD"])
    if not exists(settings["Options"]["DFTBPlusPath"]):
        msg = "DFTBPlusPath executable at %s does not exist."
        raise InputError(msg % settings["Options"]["DFTBPlusPath"])
    return settings


def write_template(filename):
    base_path = dirname(getfile(currentframe()))
    template_path = join(base_path, "input_template.yml")
    copy2(template_path, filename)
    print("No input file found. An example input file was "
          "written to %s.\nTerminating..." % filename)
    try:
        exit()
    except SystemExit:
        pass


if __name__ == "__main__":
    pass
