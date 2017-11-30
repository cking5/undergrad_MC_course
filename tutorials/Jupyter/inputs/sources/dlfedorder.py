"""Continaer for DL CONTROL FED order parameter descriptions

The order parameter is part of the fed section of the CONTROL
file.

There could be a different class for each order parameter type
but, at the moment, the differences are small enough that one
class only is provided.

The notable exception is "com2", which requires a potentially
complex block in the CONTROL input.

TODO: the "com1" and "com2" lines of the com2 order parameter
are not parsed at the moment.
"""

from collections import OrderedDict


PARAMETERS = ("ps", "psmc", "hardps", "volume", "temp", "beta", "com2")
"""Valid order parameter choices"""


class FEDOrderParameter(object):

    """Container for order parameters"""

    def __init__(self, name, ngrid, xmin, xmax, npow=None):

        """Initialise order parameter"""

        self.name = name
        self.ngrid = ngrid
        self.xmin = xmin
        self.xmax = xmax
        self.npow = npow


    def __str__(self):

        """Return a well-formed DL CONTROL file string"""

        strme = "fed order param {}".format(self.name)
        strme = "{!s} {} {} {}".format(strme, self.ngrid, self.xmin, self.xmax)

        if self.npow is not None:
            strme = "{!s} {}".format(strme, self.npow)

        return strme


    def __repr__(self):

        """Return the attritubes in readable form"""

        repme = "name= {!r}".format(self.name)
        repme = "{}, ngrid= {!r}, xmin= {!r}, xmax= {!r}"\
            .format(repme, self.ngrid, self.xmin, self.xmax)

        if self.npow is not None:
            repme = "{}, npow= {!r}".format(repme, self.npow)

        return "FEDOrderParameter({})".format(repme)


    @staticmethod
    def parse(dlstr):

        """Parse and return fed order parameter content

        Argument:
        dlstr (string):      single fed order parameter line

        The dlstr must be of the form

        fed order [[param]eter] name ngrid ibeg iend [npow]
        """

        tokens = dlstr.lower().split()

        try:

            if tokens[0] != "fed" or tokens[1] != "order":
                raise ValueError("Missing 'fed order'? {!r}".format(dlstr))

            # If we have an optional "parameter" token, lose it
            if tokens[2].startswith("param"):
                tokens.pop(2)

            name = tokens[2]
            if name not in PARAMETERS:
                raise ValueError("No order parameter; {!r}.".format(name))

            ngrid = int(tokens[3])
            xmin = float(tokens[4])
            xmax = float(tokens[5])

            try:
                npow = int(tokens[6])
            except IndexError:
                # assume optional arg not present
                npow = None

        except (ValueError, IndexError) as err:
            msg = "{} Expected 'fed order [param[eter]] ...'; got {!r}"\
                .format(err, dlstr)
            raise ValueError(msg)

        return name, ngrid, xmin, xmax, npow


    @classmethod
    def from_string(cls, dlstr):

        """Genreate instance from DL CONTROL line"""

        name, ngrid, xmin, xmax, npow \
            = FEDOrderParameter.parse(dlstr)

        if name == "com2":
            raise ValueError("com2 should be OrderCentreOfMass2")

        return FEDOrderParameter(name, ngrid, xmin, xmax, npow)


class OrderCentreOfMass2(FEDOrderParameter):

    """Centre of Mass order parameter container"""

    def __init__(self, name, ngrid, xmin, xmax, npow=None, com1=None,
                 com2=None, ncorrect=None):

        """Initialise container content

        Arguments:
        name (string):            the name "com2" (for superclass __init__)
        xmin (float):
        xmax (float):
        npow (integer):           optional
        com1 (OrderedDict):       "com1" line
        com2 (OrderedDict):       "com2" line
        ncorrect (integer):       optional sampling correction
        """

        super(OrderCentreOfMass2, self).__init__(name, ngrid, xmin, xmax, npow)

        self.com1 = OrderedDict()
        self.com2 = OrderedDict()
        self.ncorrect = ncorrect

        if com1 is not None:
            self.com1 = com1
        if com2 is not None:
            self.com2 = com2


    def __str__(self):

        """Return string appropriate for CONTROL block"""

        listme = []
        listme.append(super(OrderCentreOfMass2, self).__str__())
        listme.append(self.com1["com1"])
        listme.append(self.com2["com2"])
        if self.ncorrect is not None:
            listme.append("com sampling correction {}".format(self.ncorrect))

        listme.append("fed order param done")

        return "\n".join(listme)


    @classmethod
    def from_string(cls, dlstr):

        """Generate instance from DL CONTROL string:

        fed order [[param]eter] com2 ngrid ibeg iend [npow]
           com1 mol[elcules] <set> [atoms <set>]
           com2 mol[elcules] <set> [atoms <set>]
           [com sampling correction ncorrect]
        fed order parameter done
        """

        lines = dlstr.splitlines()

        name, ngrid, xmin, xmax, npow = FEDOrderParameter.parse(lines[0])

        try:
            com1 = OrderedDict({"com1": lines[1]})
            com2 = OrderedDict({"com2": lines[2]})
            com_sampling = lines[3]

        except IndexError:
            msg = "Expected structured order parameter block; got {!r}"\
                .format(dlstr)
            raise ValueError(msg)


        # The third line could be sampling, or just the "fed order done"

        try:
            tokens = com_sampling.lower().split()
            ncorrect = None
            if tokens[0] == "com":
                ncorrect = int(tokens[3])
        except (IndexError, ValueError):
            usage = "com sampling correction <ncorrect>"
            msg = "com2 expects {!r}; got {!r}".format(usage, com_sampling)
            raise ValueError(msg)

        return OrderCentreOfMass2(name, ngrid, xmin, xmax, npow, \
                                      com1, com2, ncorrect)


def from_string(dlstr):

    """Generate an instance from DL CONTROL entry"""

    line = dlstr.splitlines()[0]
    name = FEDOrderParameter.parse(line)[0]

    # One special case
    if name == "com2":
        return OrderCentreOfMass2.from_string(dlstr)

    return FEDOrderParameter.from_string(dlstr)

