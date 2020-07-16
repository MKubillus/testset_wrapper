class UnitConversion:
    """Unit conversions factors."""
    au2kcal = 627.509
    kcal2au = 1 / au2kcal
    au2ev = 27.211386245988
    ev2au = 1 / au2ev


class Hydrogen:
    """Energy of the H2 molecule."""
    eat = -109.87889987999591
    total_energy = -460.8595304429359


# DFTB's atomic electronic energies in kcal/mol
atomic_energies = {"c": -906.10589010466,
                   "h": -175.49031528147,
                   "n": -1391.15578889586,
                   "o": -1971.61300945930}

# DFTB3 Hubbard Derivatives for input files
hubbard_derivatives = {
    "br": -0.0573,
    "c": -0.1492,
    "ca": -0.0340,
    "cl": -0.0697,
    "f": -0.1623,
    "h": -0.1857,
    "i": -0.0433,
    "k": -0.0339,
    "mg": -0.02,
    "n": -0.1535,
    "na": -0.0454,
    "o": -0.1575,
    "p": -0.14,
    "s": -0.11,
    "zn": -0.03,
}
