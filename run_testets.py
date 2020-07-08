from libtestset import dftbplus_runner
from libtestset import results, input_parser


def main():
    settings = input_parser.load("testsets_config.yml")
    calculations = dict()
    for testset in settings["Testsets"].keys():
        set_definition = settings["Testsets"][testset]
        hsd = settings["Options"]["DFTBPlusHSD"]
        dftbplus = settings["Options"]["DFTBPlusPath"]
        calculations[testset] = dftbplus_runner.run_testset(set_definition,
                                                            hsd, dftbplus)
    results.write_results(settings, calculations)


if __name__ == "__main__":
    main()
