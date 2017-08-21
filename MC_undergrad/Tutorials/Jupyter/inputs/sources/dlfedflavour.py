"""Container for DL-MONTE FED flavour option parameters

The class structure is:

FEDFlavour
  Generic
  PhaseSwitch

Each concrete class provides a class method from_string() method to
generate an instance from the appropriate DL CONTROL file entry,
while the __str__() method returns a valid string of the same form.

The DL-MONTE internal representation is in fed_interface_type.f90
"""

from collections import OrderedDict


FLAVOURS = ("gen", "generic", "ps")
"""The available flavours, key strings thereof"""


class FEDFlavour(object):

    """Abstract container for DL-MONTE FED flavour"""

    def __init__(self, nfreq=None, keys=None):

        """Initialise container content.

        Arguments:

        nfreq (integer):       frequency of fed (if specified on use fed line)
        keys (OrderedDict):    key/values describing further fed structure 
        """

        self.nfreq = None
        self.keys = OrderedDict()

        if nfreq is not None:
            self.nfreq = nfreq
        if keys is not None:
            self.keys = keys


    def __str__(self):

        """Implementeted by subclasses"""

        raise NotImplementedError("Should be implemented by subclass")


    @classmethod
    def from_string(cls, dlstr):

        """Implementated by subclasses"""

        raise NotImplementedError("Should be implemented by subclass")


    def __repr__(self):

        """Return current state"""

        repme = "nfreq= {!r}".format(self.nfreq)

        for key in self.keys:
            repme += ", {}= {!r}".format(key, self.keys[key])

        return "{}({})".format(type(self).__name__, repme)


    @staticmethod
    def _parse_use_fed_line(dlstr):

        """Parse: 'use fed <flavour> [nfreq]' and return flavour, nfreq"""

        try:
            tokens = dlstr.lower().split()
            flavour = tokens[2]

            if tokens[0] != "use" or tokens[1] != "fed":
                raise ValueError()

            if flavour not in FLAVOURS:
                raise ValueError()

            try:
                nfreq = int(tokens[3])
            except IndexError:
                # assume optional argument not present
                nfreq = None

        except (ValueError, IndexError):
            usage = "use fed <flavour> [nfreq]"
            raise ValueError("Expected {!r}; got {!r}".format(usage, dlstr))

        return flavour, nfreq


class Generic(FEDFlavour):

    """Generic flavour FED container"""

    _defaults = {"nfreq": 1}

    def __str__(self):

        """Return the DL-MONTE CONTROL file string form"""

        strme = "use fed generic"

        if self.nfreq is not None:
            strme = "{} {}".format(strme, self.nfreq)

        return strme


    @classmethod
    def from_string(cls, dlstr):

        """Genrete an instance form DL CONTROL line"""

        flavour, nfreq = FEDFlavour._parse_use_fed_line(dlstr)

        if flavour != "gen" and flavour != "generic":
            usage = "use fed gen[eric] [nfreq]"
            raise ValueError("Expected {!r}; got {!r}".format(usage, dlstr))

        return Generic(nfreq)


class PhaseSwitch(FEDFlavour):

    """Phase Switch container object following psmc_control_type.f90"""

    # Here's a dict of allowed keywords (with default values)

    _defaults = {"nfreq": 1, \
                 "switchfreq": 0, \
                 "initactive": 1, \
                 "datafreq": 100, \
                 "meltcheck": True, \
                 "meltthresh": 10, \
                 "meltfreq": 1000}

    def __str__(self):

        """Returns a well-formed DL-CONTROL file entry"""

        listme = []
        if self.nfreq is None:
            listme.append("use fed ps")
        else:
            listme.append("use fed ps {}".format(self.nfreq))

        for key in self.keys:
            # "meltcheck" appears without a value
            if key == "meltcheck":
                listme.append("  meltcheck")
            else:
                listme.append("  {} {}".format(key, self.keys[key]))

        listme.append("ps done")

        return "\n".join(listme)


    @classmethod
    def from_string(cls, dlstr):

        """Generate instance from DL CONTROL file block

        Arguments:
        dlstr (string):   lines with blank lines and comments removed,

        which should look like:

        use fed ps [nfreq]
          keyword1 value1
          keyword2 value2
          ...
        ps done
        """

        lines = dlstr.splitlines()
        line = lines.pop(0)

        flavour, nfreq = FEDFlavour._parse_use_fed_line(line)
        keys = OrderedDict()

        if flavour != "ps":
            usage = "use fed ps [nfreq]"
            raise ValueError("Expected {}; got {!r}".format(usage, line))

        done = False

        try:

            while not done:
                line = lines.pop(0)
                tokens = line.lower().split()

                if tokens[0] == "ps" and tokens[1] == "done":
                    done = True
                    break

                key = tokens[0]

                if key == "switchfreq":
                    item = {"switchfreq": int(tokens[1])}
                elif key.startswith("init"):
                    item = {"initactive": int(tokens[1])}
                elif key == "datafreq":
                    item = {"datafreq": int(tokens[1])}
                elif key == "meltcheck":
                    item = {"meltcheck": True}
                elif key == "meltthresh":
                    item = {"meltthresh": float(tokens[1])}
                elif key == "meltfreq":
                    item = {"meltfreq": int(tokens[1])}
                else:
                    # Get out of this loop and fail
                    break

                keys.update(item)

        except (ValueError, IndexError):
            raise ValueError("Parsing failed at line {!r}".format(line))

        if not done:
            msg = "Unrecognised PSMC keyword encountered before 'ps done'"
            raise ValueError("{}: {!r}".format(msg, line))

        return PhaseSwitch(nfreq, keys)

