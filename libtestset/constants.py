class UnitConversion:
    """Unit conversions factors."""
    au2kcal = 627.509
    kcal2au = 1 / au2kcal


class Hydrogen:
    """Energy of the H2 molecule."""
    eat = -109.87889987999591
    total_energy = -460.8595304429359


# DFTB's atomic electronic energies in kcal/mol
atomic_energies = {"c": -906.10589010466,
                   "h": -175.49031528147,
                   "n": -1391.15578889586,
                   "o": -1971.61300945930}
