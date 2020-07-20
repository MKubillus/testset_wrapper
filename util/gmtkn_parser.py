"""Script that reads GMTKN55 data and creates a input template.

Might not work for earlier or future versions of the GMTKN collection.
"""

import os
from os.path import join, exists
from shutil import copy2


def parse_gmtkn_55(data_path, output_path):
    print(f"parsing from {data_path}, writing to {output_path}")
    subsets = get_subsets(data_path)
    input_base = ((f"Options:\n  # Input file to use for all calculations.\n"
                   "  # All skf paths have to be absolute and the gen input "
                   "geometry has to be set\n  # to \"in.gen\" since the "
                   "program will generate that file.\n  DFTBPlusHSD: "
                   "\"dftb_in.hsd\"\n  DFTBPlusPath: \"dftb+\"  \n"
                   "Testsets:\n"))
    for subset in subsets:
        print(f"parsing subset {subset}")
        subpath = join(data_path, subset)
        reactions = reactions_from_file(subset, subpath)
        input_base += subset_to_config(reactions, subset)
    print(f"writing config file to {join(output_path, 'testsets_config.yml')}")
    with open(join(output_path, "testsets_config.yml"), "w") as out_file:
        out_file.write(input_base)
    return input_base


def reactions_from_file(subset, subpath):
    refname = f"{subset}.ref"
    reactions = dict()
    with open(join(subpath, refname), "r") as file:
        for line in file:
            split = line.rstrip().split()
            if split[0] == "#":
                continue
            _idx = int(split[0])
            ref = float(split[-1])
            nreactants = len(split[1:-1])//2
            reactants = split[1:nreactants+1]
            stoich = split[nreactants+1:-1]
            eq = f""
            prev = 0
            for factor, reactant in sorted(zip(stoich, reactants)):
                factor = int(factor)
                if factor > 0 > prev:
                    eq += "-> "
                eq += f"{abs(factor)} {reactant} "
                prev = factor
            reactions[eq] = ref
    return reactions


def subset_to_config(reactions, subset):
    template_string = ((f"  {subset}:  # Testset name\n"
                        "    path: \"{subset}_geom\"  \n"
                        "    type: \"reaction\"\n"
                        "    reactions: \n"))
    for eq, ref in reactions.items():
        template_string += (f"      - equation: \"{eq}\"\n      "
                            f"reference: \"{ref}\" \n")
    return template_string


def get_subsets(data_path):
    subsets = []
    for f in os.listdir(data_path):
        if "backup" in f:
            continue
        else:
            subsets.append(f)
    return subsets


def write_geometries(subpath, output_path):
    subname = subpath.split("/")[-1]
    for f in os.listdir(subpath):
        if ".res" in f or "ref" in f:
            continue
        else:
            old_file = join(subpath, f, "struc.xyz")
            newpath = join(output_path, subname)
            if not exists(newpath):
                os.makedirs(newpath)
            new_file = join(newpath, f"{f}.xyz")
            copy2(old_file, new_file)


if __name__ == "__main__":
    config_string = parse_gmtkn_55("path/to/GMTKN55",
                                   "output/path/GMTKN55_parsed")
