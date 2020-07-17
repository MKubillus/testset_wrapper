#!/bin/python3

from libtestset import dftbplus_runner, dtnn_runner
from libtestset import results, input_parser


def main():
    settings = input_parser.load("testsets_config.yml")
    dftb_calcs = dict()
    dtnn_calcs = None
    hsd = settings["Options"]["DFTBPlusHSD"]
    dftbplus = settings["Options"]["DFTBPlusPath"]
    dtnn = ("DTNN" in settings["Options"])
    if dtnn:
        dtnn_calcs = dict()
        model = settings["Options"]["DTNN"]["DTNNModel"]
        skf = settings["Options"]["DTNN"]["DTNNSkfPath"]
    for testset in settings["Testsets"]:
        set_definition = settings["Testsets"][testset]
        dftb_calcs[testset] = dftbplus_runner.run_testset(set_definition,
                                                          hsd, dftbplus)
        if dtnn:
            dtnn_calcs[testset] = dtnn_runner.run_testset(
                set_definition, model, dftbplus, skf)
    results.write_results(settings["Testsets"], dftb_calcs, dtnn_calcs)


if __name__ == "__main__":
    main()
