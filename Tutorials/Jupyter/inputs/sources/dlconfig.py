"""A represenation of the DL-MONTE CONFIG file

This is simple enough that one class is defined to take all
the content (including molcules/atoms).

TODO. Atoms are only stored as strings; a more structured
format could be imagined.

A module level convenience methof from_file() is supplied
as a factory method.
"""

from collections import OrderedDict
import json
import numpy

import inputs.sources.dlutil as dlutil

DLMONTE_FORMAT_FRACTIONAL = 0
DLMONTE_FORMAT_CARTESIAN = 1
DLPOLY_LEVEL_ZERO = 0
DLPOLY_LEVEL_ONE = 1
DLPLOY_LEVEL_TWO = 2

class CONFIG(object):

    """Container class for CONFIG"""

    def __init__(self, title=None, level=None, dlformat=None, vcell=None,
                 nummol=None, molecules=None):

        """Representation of CONFIG input file
        """

        self.title = title
        self.level = level
        self.dlformat = dlformat
        self.vcell = vcell

        self.nummol = []
        self.molecules = []

        if nummol is not None:
            self.nummol = nummol
        if molecules is not None:
            self.molecules = molecules


    def __repr__(self):

        """Return a readable format"""

        me1 = "title= {!r}, level= {!r}, dlformat= {!r}" \
            .format(self.title, self.level, self.dlformat)
        me2 = "vcell= {!r}, nummol= {!r}".format(self.vcell, self.nummol)

        return "CONFIG({!s}, {!s})".format(me1, me2)


    def __str__(self):

        """Return string format which is a CONFIG file"""

        lines = []
        lines.append(self.title)
        lines.append("{} {}".format(self.level, self.dlformat))

        # Cell vectors
        for dim in range(3):
            x, y, z = self.vcell[dim]
            lines.append("{} {} {}".format(x, y, z))

        nmaxstr = " ".join(str(n) for n in self.nummol)
        lines.append("NUMMOL {}".format(nmaxstr))

        # Molecules and constituent lines
        for mol in self.molecules:
            lines.append("MOLECULE {} {} {}"\
                         .format(mol["name"], mol["natom"], mol["nmaxatom"]))

            for line in mol["atoms"]:
                lines.append(line)

        return "\n".join(lines)


    def to_dct(self):

        """Render as dict in keeping with DL-Monte style"""

        dct = OrderedDict()
        dct.update({"TITLE" : self.title})
        dct.update({"FORMAT" : self.dlformat})
        dct.update({"LEVEL" : self.level})
        dct.update({"LATVECTOR" : self.vcell})

        return dct


    def to_json(self):

        """Return JSON"""

        dct = self.to_dct()
        return json.dumps(dct, indent=2)


    @property # a property for backward compatibility reasons
    def natom(self):

        """Return the total number of atoms"""

        natom = 0

        for mol in self.molecules:
            natom += mol["natom"]

        return natom


    def volume(self):

        """Return the volume of the system based on lattice vectors"""

        a = self.vcell[0]
        b = self.vcell[1]
        c = self.vcell[2]

        return numpy.dot(a, numpy.cross(b, c))


    @classmethod
    def from_file(cls, filename, dlfield=None):

        """Return an instance from file"""

        lines = dlutil.load_ascii(filename)
        config = CONFIG.from_dlstr("\n".join(lines), dlfield)

        return config


    @classmethod
    def from_dlstr(cls, dlstr, dlfield=None):

        """Generate instance from content of DL CONFIG file.

        It is assumed comments and blank lines have been stripped out.

        dlstr (string):          the content
        dlfield (FIELD):         If present, the content can be checked
                                 against the field description

        Returns: new instance of CONFIG
        """

        lines = dlstr.splitlines()

        # "title"
        # "level format"

        try:
            title = lines.pop(0)
            level, dlformat = [int(n) for n in lines.pop(0).split()]
        except ValueError:
            raise

        # System initial cell vectors a, b, c - each (x,y,z)

        vcell = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        try:
            vcell[0] = list(float(x) for x in lines.pop(0).split())
            vcell[1] = list(float(x) for x in lines.pop(0).split())
            vcell[2] = list(float(x) for x in lines.pop(0).split())
        except ValueError:
            raise

        # Molecule information
        # NUMMOL ntot n1 n2 ...
        # ntot is total number of molecules in config file
        # n1 is maximum number of molecules moltype 1 allowed,
        # n2 is maximum number of molcules  moltype 2 allowed,
        # and so on,

        try:
            line = lines.pop(0).split()
            if line[0].lower() != "nummol":
                raise ValueError("Expecting NUMMOL ntot n1 n2 ...")

            ntot = int(line[1])
            nmolmax = list(int(n) for n in line[2:])
        except (IndexError, ValueError):
            raise

        # For each molecule in the configuration read:
        # "molecule name natompresent natommax"
        #   and for each atom present read
        #     "type"
        #     "x y z"

        # At the moment, only one line per atom record (x,y,z) is expected
        assert level == DLPOLY_LEVEL_ZERO

        molecules = []

        for _ in range(ntot):

            mol = OrderedDict()
            try:
                line = lines.pop(0)
                tokens = line.split()
                mol.update({"name": str(tokens[1])})
                mol.update({"natom": int(tokens[2])})
                mol.update({"nmaxatom": int(tokens[3])})
            except (IndexError, ValueError):
                raise ValueError("Could not parse MOLECULE: {!r}|".format(line))

            atoms = []
            for n in range(mol["natom"]):
                # The atoms are just the string descriptions
                atoms.append(lines.pop(0))
                atoms.append(lines.pop(0))

            mol.update({"atoms": atoms})
            molecules.append(mol)


        if dlfield is not None:
            assert len(nmolmax) == len(dlfield.moltypes)

            for m in range(len(dlfield.moltypes)):
                assert molecules[m]["name"] == dlfield.moltypes[m].name


        nummol = list([ntot])
        for nmol in nmolmax:
            nummol.append(nmol)

        return CONFIG(title, level, dlformat, vcell, nummol, molecules)


def from_file(filename, dlfield=None):

    """A convenience method to generate a config"""

    return CONFIG.from_file(filename, dlfield)
