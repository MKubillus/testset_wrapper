from libtestset import dftbplus_runner
from libtestset import results, input_parser


def main():
    settings = input_parser.load("testsets_config.yml")
    calculations = dict()
    hsd = settings["Options"]["DFTBPlusHSD"]
    dftbplus = settings["Options"]["DFTBPlusPath"]
    for testset in settings["Testsets"]:
        set_definition = settings["Testsets"][testset]
        calculations[testset] = dftbplus_runner.run_testset(set_definition,
                                                            hsd, dftbplus)
    results.write_results(settings["Testsets"], calculations)


if __name__ == "__main__":
    main()
