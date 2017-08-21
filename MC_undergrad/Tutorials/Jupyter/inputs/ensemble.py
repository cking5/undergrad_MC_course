"""Monte Carlo Ensemble

A classification of available Ensemble types.
"""

class Ensemble(object):

    """Abstract Ensemble class"""

    def __init__(self, name, longname):

        self.name = name
        self.longname = longname


    def __delattr__(self, name):

        msg = "Cannot delete {!s}.{!s}".format(type(self).__name__, name)
        raise AttributeError(msg)


    def __hash__(self):

        return hash(self.name) ^ hash(self.longname)


    def __eq__(self, other):

        return self.name == other.name and self.longname == other.longname


    def __str__(self):

        return "{!s} ({!s})".format(self.name, self.longname)


    def __repr__(self):

        strme = "name= {!r}, longname= {!r}".format(self.name, self.longname)
        return "{!s}({!s})".format(type(self).__name__, strme)


class EnsembleNPT(Ensemble):

    """Isothermal Isobaric ensemble (NPT)"""

    def __init__(self):

        super(EnsembleNPT, self).__init__("NPT", "Isothermal-Isobaric")


class EnsembleNVT(Ensemble):

    """Canonical Ensemble (NVT)"""

    def __init__(self):

        super(EnsembleNVT, self).__init__("NVT", "Canonical")


class EnsembleMuVT(Ensemble):

    """Grand Canonical Ensemble (muVT)"""

    def __init__(self):

        super(EnsembleMuVT, self).__init__("muVT", "Grand Canonical")


class EnsembleNVE(Ensemble):

    """Microcanonical Ensemble (NVE)"""

    def __init__(self):

        super(EnsembleNVE, self).__init__("NVE", "Microcanonical")


KEY_ENSEMBLE = {"nvt": EnsembleNVT, "npt": EnsembleNPT, "muvt": EnsembleMuVT,
                "nve": EnsembleNVE}

def from_string(dlstr):

    """Return an Ensemble"""

    key = dlstr.lower()

    if key not in KEY_ENSEMBLE:
        raise ValueError("Unrecognised Enseble: {!r}".format(dlstr))

    return KEY_ENSEMBLE[key]()
