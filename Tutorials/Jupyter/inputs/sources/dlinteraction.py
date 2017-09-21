"""DL-MONTE Interaction descriptions

Interaction types are merely to identify and to hold the parameters
associated with the different types of bonds, angles, etc.

Some of the interations may appear as either bonded or non-bonded,
other not; this is not really captured here.

No computation takes place.
"""

from collections import OrderedDict


class Interaction(object):

    """The interaction abstract class, or interface

    Attributes:
      key (string): each concrete class has a key which identifies
      it in the DL-MONTE FIELD file format, e.g., "lj" for Lennard-
      Jones.
    """

    def to_dct(self):

        """Implemented by child classes"""

        raise NotImplementedError()

    @classmethod
    def from_string(cls, dlstr):

        """Implmented by child classes as factory method"""

        raise NotImplementedError()


class InteractionLJ(Interaction):

    """Standard Lennard Jones potential (bonded or non-bonded)

    U(r) = 4 epsilon [ (sigma/r)^12 - (sigma/r)^6 ]
    """

    key = "lj"

    def __init__(self, epsilon, sigma):

        """
        Arguemnts:
          epsilon (float): the energy scale
          sigma   (float): length scale
        """

        self.key = InteractionLJ.key
        self.type = "Lennard-Jones"
        self.epsilon = epsilon
        self.sigma = sigma

    def __str__(self):

        """Return (sub-) string relevant for DL FIELD file"""

        return "{!s} {!s} {!s}".format(self.key, self.epsilon, self.sigma)


    def __repr__(self):

        """Return internal representation"""

        rep = "key={!r}, type={!r}, epsilon={!r}, sigma={!r}" \
            .format(self.key, self.type, self.epsilon, self.sigma)
        return "Interaction({!s})".format(rep)


    def to_dct(self):

        """Construct and return a dictionary description"""

        dct = OrderedDict()
        dct.update({"KEY" : self.key})
        dct.update({"EPSILON" : self.epsilon})
        dct.update({"SIGMA" : self.sigma})

        return dct


    @classmethod
    def from_string(cls, dlstr):

        """Generate and return an object from DL-MONTE style string

        Arguments:
          dlstr (string): description

        Example:
        "lj 1.0 2.0" gives epsilon = 1.0, sigma = 2.0

        """

        try:
            name, epsilon, sigma = dlstr.split()
            if name.lower() != cls.key:
                raise ValueError()

        except ValueError:
            raise ValueError("Require `lj eps sigma` not {!r}".format(dlstr))

        epsilon = float(epsilon)
        sigma = float(sigma)

        return cls(epsilon, sigma)


class InteractionNM(Interaction):

    """Generalised Lennard-Jones potential

    U = (e0/(n-m)) x [m(r0/r)^n - n(r0/r)^m]

    """

    key = "nm"

    def __init__(self, e0, n, m, r0):

        """Initialise description

        Arguments:
          e0 (float): an energy
          n  (int):   exponent (+ve term)
          m  (int):   exponent (-ve term)
          r0 (float): a length
        """

        self.key = InteractionNM.key
        self.type = "Lennard Jones General N-M"
        self.e0 = e0
        self.n = n
        self.m = m
        self.r0 = r0


    def __str__(self):

        """Return string for DL FIELD file"""

        dlstr = "{!s} {!s} {!s} {!s} {!s}"\
            .format(self.key, self.e0, self.n, self.m, self.r0)

        return dlstr


    def __repr__(self):
        me = "key={!r}, type={!r}, e0={!r}, n={!r}, m={!r}, r0={!r}" \
            .format(self.key, self.type, self.e0, self.n, self.m, self.r0)
        return "Interaction({!s})".format(me)



    def to_dct(self):
        """Construct and return a dictionary description"""

        dct = OrderedDict()
        dct.update({"KEY" : self.key})
        dct.update({"E_0" : self.e0})
        dct.update({"N" : self.n})
        dct.update({"M" : self.m})
        dct.update({"R_0" : self.r0})

        return dct


    @classmethod
    def from_string(cls, dlstr):
        """Generate and return an object from a DL-MONTE style string

        Example:
        dlstr = "nm 1.0 12 6 2.5"
        """
        name, e0, n, m, r0 = dlstr.split()
        assert name == InteractionNM.key
        e0 = float(e0)
        n = int(n)
        m = int(m)
        r0 = float(r0)
        return cls(e0, n, m, r0)


class Interaction12_6(Interaction):

    """Lennard-Jones-like potential (bonded or non-bonded)

    U = (a/r^12) - (b/r^6)
    """

    key = "12-6"

    def __init__(self, a, b):
        """
        Arguments:
          a (float):  parameter in r^12 term
          b (float):  parameter in r^6 term
        """
        self.key = Interaction12_6.key
        self.type = "Lennard-Jones-like 12-6"
        self.a = a
        self.b = b


    def __str__(self):

        """Return (sub-) string in DL FIELD style"""

        return "{!s} {!s} {!s}".format(self.key, self.a, self.b)


    def __repr__(self):
        me = "key={!r}, type={!r}, a={!r}, b={!r}" \
            .format(self.key, self.type, self.a, self.b)
        return "Interaction({!s})".format(me)


    def to_dct(self):
        """Construct and return a dictionary description"""

        dct = OrderedDict()
        dct.update({"KEY" : self.key})
        dct.update({"A" : self.a})
        dct.update({"B" : self.b})

        return dct


    @classmethod
    def from_string(cls, dlstr):
        """Generate and return an object from DL-MONTE style string

        Exmaple:
        dlstr = "12-6 1.0 2.0" gives a = 1.0 and b = 2.0
        """

        name, a, b = dlstr.split()
        assert name == cls.key
        a = float(a)
        b = float(b)
        return cls(a, b)


class InteractionBuckingham(Interaction):

    """Buckingham potential (bonded or non-bonded)

    U(r) = A exp [-r/rho] - C/r^6
    """

    key = "buck"

    def __init__(self, a, rho, c):

        """
        Arguments:
          a    (float): an energy
          rho  (float): a length
          c    (float): an energy
        """

        self.key = InteractionBuckingham.key
        self.type = "Buckingham"
        self.a = a
        self.rho = rho
        self.c = c


    def __str__(self):

        """Return a string in DL FIELD style"""

        dlstr = "{!s} {!s} {!s} {!s}" \
            .format(self.key, self.a, self.rho, self.c)
        return dlstr


    def __repr__(self):
        me = "key={!r}, type={!r}, a={!r}, rho={!r}, c={!r}" \
            .format(self.key, self.type, self.a, self.rho, self.c)
        return "Interaction({!s})".format(me)


    def to_dct(self):
        """Construct and return a dictionary description"""

        dct = OrderedDict()
        dct.update({"KEY" : self.key})
        dct.update({"A" : self.a})
        dct.update({"rho" : self.rho})
        dct.update({"C" : self.c})

        return dct


    @classmethod
    def from_string(cls, dlstr):
        """Generate and return an object from a DL-MONTE style string

        Example:
        dlstr = "buck 1.0 2.0, 3.0" gives a = 1.0, rho = 2.0, and c = 3.0
        """

        name, a, rho, c = dlstr.split()
        assert name == cls.key
        a = float(a)
        rho = float(rho)
        c = float(c)
        return cls(a, rho, c)


class InteractionHS(Interaction):

    """Hard Spheres

    U(r) = Arb. large r < sigma and zero if r >= sigma, where sigma
    is the hard sphere radius.
    """

    key = "hs"

    def __init__(self, sigma):

        """Initilaise descirption

        Arguments:
          sigma (float):  hard-sphere cut off
        """

        self.key = InteractionHS.key
        self.type = "Hard Sphere"
        self.sigma = sigma

    def __str__(self):

        """Return string in DL FIELD style"""

        return "{!s} {!s}".format(self.key, self.sigma)


    def __repr__(self):
        me = "key={!r}, type={!r}, sigma={!r}" \
            .format(self.key, self.type, self.sigma)
        return "Interaction({!s})".format(me)


    def to_dct(self):
        """Construct and return OrderedDict"""

        dct = OrderedDict()
        dct.update({"KEY" : self.key})
        dct.update({"SIGMA" : self.sigma})

        return dct

    @classmethod
    def from_string(cls, dlstr):

        """Genreate and return object from DL-style string

        Example:
        dlstr = "hs 1.5" gives sigma = 1.5
        """

        try:
            name, sigma = dlstr.split()
            if name.lower() != cls.key:
                raise ValueError()

        except ValueError:
            raise ValueError("Require `hs sigma` not {!r}".format(dlstr))

        sigma = float(sigma)
        return cls(sigma)


def from_string(dlstr):

    """Factory method taking the string as appearing in FIELD input"""

    input_key = dlstr.split()[0].lower()

    for subcls in Interaction.__subclasses__():
        if subcls.key == input_key:
            return subcls.from_string(dlstr)

    raise ValueError("No interaction available for {!r}".format(str))
