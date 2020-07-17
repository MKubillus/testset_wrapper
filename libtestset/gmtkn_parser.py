import os
from os.path import join, exists
from shutil import copy2
import re

def parse_gmtkn_55(data_path, output_path):
    print(f"parsing from {data_path}, writing to {output_path}")
    subsets = get_subsets(data_path)
    config_string = (f"Options:\n"
                     "  # Input file to use for all calculations.\n"
                     "  # All skf paths have to be absolute and the gen input geometry has to be set\n"
                     "  # to \"in.gen\" since the program will generate that file.\n"
                     "  DFTBPlusHSD: \"dftb_in.hsd\"\n"
                     "  DFTBPlusPath: \"dftb+\"  \n"
                     "Testsets:\n")
    for subset in subsets:
        print(f"parsing subset {subset}")
        subpath = join(data_path, subset)
        geom_location = write_geometries(subpath, output_path)
        reactions = reactions_from_file(subset, subpath)
        config_string += subset_to_config(reactions, subset)
    print(f"writing config file to {join(output_path, 'testsets_config.yml')}")
    with open(join(output_path, "testsets_config.yml"), "w") as out_file:
        out_file.write(config_string)
    return config_string

def reactions_from_file(subset, subpath):
    refname = f"{subset}.ref"
    reactions = dict()
    with open(join(subpath, refname), "r") as file:
        for l in file:
            split = l.rstrip().split()
            if split[0]=="#":
                continue
            idx = int(split[0])
            ref = float(split[-1])
            nreactants = len(split[1:-1])//2
            reactants = split[1:nreactants+1]
            stoich = split[nreactants+1:-1]
            eq = f""
            prev = 0
            for factor, reactant in sorted(zip(stoich, reactants)):
                factor = int(factor)
                if factor > 0 and prev < 0:
                    eq += "-> "
                eq += f"{abs(factor)} {reactant} "
                prev = factor
            reactions[eq] = ref
    return reactions

def subset_to_config(reactions, subset):
    template_string = (f"  {subset}:  # Testset name\n"
                       f"    path: \"{subset}_geom\"  \n"
                        "    type: \"reaction\"\n"
                        "    reactions: \n")
    for eq, ref in reactions.items():
        template_string += f"      - equation: \"{eq}\"\n      reference: \"{ref}\" \n"
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
            #print(f"copying {old_file} to {new_file}")
            copy2(old_file, new_file)
        break

if __name__=="__main__":
    config_string = parse_gmtkn_55("/media/mila/HD1/GMTKN55", "/media/mila/HD1/GMTKN55_parsed")
