"""Containers for DL-MONTE FED Method descriptions

The FED Method is part of the FED section of the CONTROL file.

FEDMethod provides a classification of the various concrete
methods available, each of which has its own class.

FEDMethod:
  UmbrellaSampling
  WangLandau
  ExpandedEnsemble
  TransitionMatrix

A module level from_string() method is supplied to allow
creation from a well formed string, and the relevant __str__()
method should return the same string that can be used in the
CONTROL file.
"""


class FEDMethod(object):

    """Container for an FED method"""

    def __str__(self):
        """To be implemented by subclasses"""
        raise NotImplementedError("Should be implemented by subclass")

    @classmethod
    def from_string(cls, dlstr):
        """To be implemented by subclasses"""
        raise NotImplementedError("Should be implemented by subclass")


class BiasSmoother(object):

    """Container for Method bias smoother parameters"""

    def __init__(self, n_itr, i_beg, i_end, omega):

        """Arguments:

        n_itr (int):    number of bias smoothing iterations 1 <= n_itr <= 5
        i_beg (int):    start of bias grid
        i_end (int):    end of bias grid
        omega (float):  weight of central bin
        """

        self.n_itr = n_itr
        self.i_beg = i_beg
        self.i_end = i_end
        self.omega = omega

    def __repr__(self):

        """Return something readable"""

        repme = "n_itr= {!r}, i_beg= {!r}, i_end= {!r}, omega= {!r}"\
            .format(self.n_itr, self.i_beg, self.i_end, self.omega)

        return "BiasSmoother({!s})".format(repme)


    def __str__(self):

        """Return a string in DL CONTROL style"""

        # This appears at the end of the fed method line
        strme = "{} {} {} {}"\
            .format(self.n_itr, self.i_beg, self.i_end, self.omega)

        return strme


class UmbrellaSampling(FEDMethod):

    """Umbrella sampling CONTROL FILE entry"""

    key = "us"

    def __init__(self, x0, kf, n_upd):

        """Initialise container content

        Arguments:
            x0 (float):             bias minimum
            kf (float):             bias spring constant
            n_upd (int):            update frequency
        """

        self.x0 = x0
        self.kf = kf
        self.n_upd = n_upd

    def __repr__(self):

        """Return parameters as readable string"""

        repme = "x0= {!r}, kf= {!r}, n_upd= {!r}"\
            .format(self.x0, self.kf, self.n_upd)

        return "UmbrellaSampling({!s})".format(repme)


    def __str__(self):

        """Return string in DL CONTROL style"""

        strme = "fed method {} {} {} {}"\
            .format(UmbrellaSampling.key, self.x0, self.kf, self.n_upd)

        return strme


    @classmethod
    def from_string(cls, dlstr):

        """Generate object from DL CONTROL file string

        e.g., dlstr = "us 1.0 0.5 1"
        """

        try:

            key, x0, kf, n_upd = dlstr.lower().split()
            if key != UmbrellaSampling.key:
                raise ValueError()

            x0 = float(x0)
            kf = float(kf)
            n_upd = int(n_upd)

        except ValueError:
            raise ValueError("Require 'us x0 kf n_upd' not {!r}".format(dlstr))

        return UmbrellaSampling(x0, kf, n_upd)


class ExpandedEnsemble(FEDMethod):

    """Container for expanded ensemble method parameters"""

    key = "ee"

    def __init__(self, eta0, c_upd, n_upd, smooth=None):

        """Initialise container parameters

        Arguments:
            eta0 (float):           initial damping factor 0 < eta0 <= 1
            c_upd (float):          scaling factor > 1
            n_upd (int):            update frequency
            smooth (BiasSmoother):  optional bias smoother parameters
        """

        self.eta0 = eta0
        self.c_upd = c_upd
        self.n_upd = n_upd
        self.smooth = smooth

    def __repr__(self):

        """Return a readable string"""

        repme = "eta0= {!r}, c_upd= {!r}, n_upd= {!r}, smooth= {!r}"\
            .format(self.eta0, self.c_upd, self.n_upd, self.smooth)

        return "ExpandedEnsemble({!s})".format(repme)


    def __str__(self):

        """Return a DL CONTROL style string"""

        strme = "fed method {} {} {} {}"\
            .format(ExpandedEnsemble.key, self.eta0, self.c_upd, self.n_upd)
        if self.smooth:
            strme = "{!s} {!s}".format(strme, self.smooth)

        return strme


    @classmethod
    def from_string(cls, dlstr):

        """Return instance from DL CONTROL file string

        e.g., dlstr = "ee eta0 c_upd n_upd [smooth]"
        """

        smooth = None

        try:
            tokens = dlstr.lower().split()

            if tokens[0] != ExpandedEnsemble.key:
                raise ValueError()

            eta0 = float(tokens[1])
            c_upd = float(tokens[2])
            n_upd = int(tokens[3])

            try:
                n_itr = int(tokens[4])
                i_beg = int(tokens[5])
                i_end = int(tokens[6])
                omega = float(tokens[7])
                smooth = BiasSmoother(n_itr, i_beg, i_end, omega)

            except IndexError:
                # assume optional arguments not present
                pass

        except (IndexError, ValueError):
            msg = "Expect 'ee eta0 c_upd u_upd []'; got {!r}".format(dlstr)
            raise ValueError(msg)

        return ExpandedEnsemble(eta0, c_upd, n_upd, smooth)


class WangLandau(FEDMethod):

    """Container for Wang Landau meothd parameters"""

    key = "wl"

    def __init__(self, delta0, c_upd, n_upd, smooth=None):

        """Initialise container parameters

        Arguemnts:
            delta0 (float):
            c_upd (float):
            n_upd (int):
            smooth (BiasSmoother):  optional bias smoother parameters
        """

        self.delta0 = delta0
        self.c_upd = c_upd
        self.n_upd = n_upd
        self.smooth = smooth


    def __repr__(self):

        """Return a readable string"""

        repme = "delta0= {!r}, c_upd= {!r}, n_upd= {!r}, smooth= {!r}"\
            .format(self.delta0, self.c_upd, self.n_upd, self.smooth)

        return "WangLandau({!s})".format(repme)


    def __str__(self):

        """Return a well-formed DL CONTROL string"""

        strme = "fed method {} {} {} {}"\
            .format(WangLandau.key, self.delta0, self.c_upd, self.n_upd)
        if self.smooth:
            strme = "{} {}".format(strme, self.smooth)

        return strme


    @classmethod
    def from_string(cls, dlstr):

        """Generate instance from DL CONTROL file string

        e.g. dlstr = "1.0 2.0 1 [bias smoother parameters]"
        """

        smooth = None

        try:

            tokens = dlstr.lower().split()

            if tokens[0] != WangLandau.key:
                raise ValueError

            delta0 = float(tokens[1])
            c_upd = float(tokens[2])
            n_upd = int(tokens[3])

            try:
                n_itr = int(tokens[4])
                i_beg = int(tokens[5])
                i_end = int(tokens[6])
                omega = float(tokens[7])
                smooth = BiasSmoother(n_itr, i_beg, i_end, omega)
            except IndexError:
                # assume optional argument not present
                pass

        except (IndexError, ValueError):
            msg = "Expect 'WL delta0 c_upd n_upd []'; got {!r}".format(dlstr)
            raise ValueError(msg)

        return WangLandau(delta0, c_upd, n_upd, smooth)


class TransitionMatrix(FEDMethod):

    """Container for transition matrix method parameters"""

    key = 'tm'

    def __init__(self, nout, n_upd, mode="new"):

        """Initialise content

        Arguments:
            nout (int):       bias freq
            n_upd (int):      output freq
            mode (string):    "new" or "resume"
        """

        self.nout = nout
        self.n_upd = n_upd
        self.mode = mode

    def __repr__(self):

        """Return content as string"""

        repme = "nout= {!r}, n_upd= {!r}, mode= {!r}"\
            .format(self.nout, self.n_upd, self.mode)

        return "TransitionMatrix({!s})".format(repme)


    def __str__(self):

        """Return well-formed DL CONTROL string"""

        strme = "fed method {} {} {} {}"\
            .format(TransitionMatrix.key, self.nout, self.n_upd, self.mode)

        return strme


    @classmethod
    def from_string(cls, dlstr):

        """Generate instance from DL CONTROL string"""

        mode = "new" # the default

        try:
            tokens = dlstr.lower().split()

            if tokens[0] != TransitionMatrix.key:
                raise ValueError

            nout = int(tokens[1])
            n_upd = int(tokens[2])

            try:
                mode = str(tokens[3])
            except IndexError:
                # assume optional argument not present
                pass

        except (IndexError, ValueError):
            usage = "fed method tm nout n_upd [mode]"
            raise ValueError("Expected {!r}: got {!r}".format(usage, dlstr))

        return TransitionMatrix(nout, n_upd, mode)


def from_string(dlstr):

    """Return a FEDMethod object from DL CONTROL file string"""

    methods = {UmbrellaSampling.key: UmbrellaSampling,
               ExpandedEnsemble.key: ExpandedEnsemble,
               WangLandau.key: WangLandau,
               TransitionMatrix.key: TransitionMatrix}

    try:
        # Remove the "fed method" and pass on the rest
        keystr = dlstr.split(None, 2)[2]
        key = keystr.split()[0].lower()

        instance = methods[key].from_string(keystr)

    except KeyError:
        raise ValueError("Unrecognised fed method: {!r}".format(dlstr))

    return instance
