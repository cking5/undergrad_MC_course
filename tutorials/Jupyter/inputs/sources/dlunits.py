"""Constants related to physical units

Values here are to be consistent with those used in
DL-MONTE constants.f90

DL_MONTE Units for input/output

Temperature    Kelvin
Length         Angstrom
Energy         Defined by user in FIELD file
Charge         Electronic charge (?)
Pressure       kilo Atmospheres
"""

from inputs.parameter import Parameter
from inputs.util import Label

LABEL_NA = Label("N_A", "Avogadro's Number", "per mole")
LABEL_KB = Label("k_b", "Boltzmann Constant", "Joules per Kelvin")
LABEL_QE = Label("q", "Electronic charge", "Coulomb")

NA_SI = Parameter(6.022140e+23, LABEL_NA)
KB_SI = Parameter(1.3806488e-23, LABEL_KB)
QE_SI = Parameter(1.602176e-19, LABEL_QE)


def k_boltzmann(dlstr):

    """Return value of Boltzmann constant in current DL units

    Args:
    dlstr (string): "ev" | 'kJ" | "kcal" | "k" | "internal"

    Returns:
    Boltzmann constant in the relevant units
    """

    key = dlstr.lower()

    label = Label("k_b", "Boltzmann Constant", None)

    if key == "ev":
        label.units = "eV per Kelvin"
        parameter = Parameter((1.0/QE_SI)*KB_SI, label)
    elif key == "kcal":
        label.units = "kCal per mole per Kelvin"
        parameter = Parameter((1.0/4184.0)*KB_SI*NA_SI, label)
    elif key == "kj":
        label.units = "kJoules per mole per Kelvin"
        parameter = Parameter((1.0/1000.0)*KB_SI*NA_SI, label)
    elif key == "k":
        label.units = "10 Joules per mole per Kelvin"
        parameter = Parameter((1.0/10.0)*KB_SI*NA_SI, label)
    elif key == "internal":
        label.units = "10 Joules per mole per Kelvin"
        parameter = Parameter((1.0/10.0)*KB_SI*NA_SI, label)
    else:
        raise ValueError("units not recognised {}".format(dlstr))

    return parameter
