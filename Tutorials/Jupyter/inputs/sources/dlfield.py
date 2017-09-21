"""DL-MONTE FIELD file decription

This module provides access to a FIELD class, which is intended
to describe the DL-MONTE FIELD input file.

The FIELD class is an aggregation of various AtomTypes, MolTypes
interactions, classes for which are also provided to impart structure.

The ability to read DL-MONTE style FIELD files is provided.
The ability to write the content in various formats, including
the DL-MONTE style itself, is also provided.

TODO: Not all interaction types are handled owing to lack of examples
in the DL MONTE regrssion tests (e.g., Three-Body interactions).

"""

from collections import OrderedDict
import json

import inputs.sources.dlutil as dlutil
import inputs.sources.dlinteraction as dlinteraction

DLMONTE_VERSION = 2.01

class Atom(object):

    """FIELD file Atom entry"""

    def __init__(self, name, atype, rpos=None):

        """Initialise container content

        Arguments:
          name (string):                  descriptive name
          atype (string):                 core|shell
          rpos ([float, float, float]):   [x,y,x] relative position
        """

        self.name = name
        self.type = atype

        if rpos is not None:
            self.rpos = rpos


    def __repr__(self):

        """Return internal representation"""

        atom = "name= {!r}, type= {!r}"\
            .format(self.name, self.type)

        if hasattr(self, "rpos"):
            atom += ", rpos= {!r}".format(self.rpos)

        return "Atom({!s})".format(atom)


    def __str__(self):

        """Return string intended to reproduce DL FIELD file string"""

        atom = "{!s} {!s}".format(self.name, self.type)

        if hasattr(self, "rpos"):
            atom = "{!s} {!s} {!s} {!s}"\
                .format(atom, self.rpos[0], self.rpos[1], self.rpos[2])

        return atom


    def to_dct(self):

        """Return a dict in keeping with DLFIELD style"""

        dct = OrderedDict()
        dct.update({"NAME" : self.name})
        dct.update({"TYPE" : self.type})
        dct.update({"RELPOS" : self.rpos})

        return dct


class AtomType(object):

    """FIELD file AtomType entry"""

    def __init__(self, name, atype, mass, charge):

        """Initialise content

        Arguments:
          name (string):       descriptive
          atype (string):      ["core"|"shell"|"metal"]
          mass (int/float):    mass
          charge (float):      charge

        Note mass for truly atomic species is generally an integer,
        while can be float. We keep the destinction.
        """

        self.name = name
        self.type = atype
        self.mass = mass
        self.charge = charge


    def __repr__(self):

        """Return internal representation"""

        atom = "name= {!r}, type= {!r}, mass= {!r}, charge= {!r}" \
            .format(self.name, self.type, self.mass, self.charge)
        return "AtomType({!s})".format(atom)


    def __str__(self):

        """Return string for DL FIELD file style output"""

        str1 = "{!s} {!s}".format(self.name, self.type)
        str2 = "{!s} {!s}".format(self.mass, self.charge)
        return "{} {}".format(str1, str2)


    def to_dct(self):

        """Return dict in the DL FIELD style"""

        dct = OrderedDict()
        dct.update({"NAME": self.name})
        dct.update({"TYPE": self.type})
        dct.update({"MASS": self.mass})
        dct.update({"CHARGE": self.charge})
        return dct


    @staticmethod
    def from_string(dlstr):

        """Generate AtomType from string"""

        try:
            name, atype, strmass, charge = dlstr.split()

            # Parse mass as float if a decimal point appears in the string
            if "." in strmass:
                mass = float(strmass)
            else:
                mass = int(strmass)
            charge = float(charge)
            atomtype = AtomType(name, atype, mass, charge)

        except (ValueError, IndexError):
            raise ValueError("Failed to parse atom type: {!r}".format(dlstr))

        return atomtype


class MolType(object):

    """FIELD file MOLTYPE entry"""

    # Structureless types have only "MAXATOM"
    # Structures have an ATOM list as well

    # Note that the bonds are specfied as an integer triple
    # relating to the order in which atoms and bond types
    # appear in the FIELD file.

    # As both python and json lists maintain order, we can
    # use the same mechanism reliably here provided the
    # file is parsed in the correct order.

    def __init__(self, name, nmaxatom):

        """Initialise content

        Arguments:
          name       string
          nmaxatom   integer  maximum number of atoms allowed in molecule
        """

        self.name = name
        self.nmaxatom = nmaxatom
        self.atoms = []
        self.bonds = []
        self.exc_coul_ints = False
        self.rigid = False


    def __repr__(self):

        """Return internal representation"""

        str1 = "name={!r}, nmaxatom={!r}, atoms={!r}, bonds={!r}" \
            .format(self.name, self.nmaxatom, self.atoms, self.bonds)
        str2 = "exc_coul_ints={!r}".format(self.exc_coul_ints)
        str3 = "rigid= {!r}".format(self.rigid)
        return "MolType({!s}, {!s}, {!s})".format(str1, str2, str3)


    def __str__(self):

        """Return string in keeping with the FIELD style"""

        listme = []
        listme.append(str(self.name))

        if not self.atoms:
            listme.append("MAXATOM {!s}".format(self.nmaxatom))
        else:
            atoms = "ATOMS {!s} {!s}".format(len(self.atoms), self.nmaxatom)
            listme.append(atoms)
            for atom in self.atoms:
                listme.append(str(atom))

        if self.bonds:
            listme.append("BONDS {!s}".format(len(self.bonds)))
            for bond in self.bonds:
                listme.append(" ".join(map(str, bond)))

        if self.exc_coul_ints:
            listme.append("EXCLUDE")
        if self.rigid:
            listme.append("RIGID")

        return "\n".join(listme)


    def to_dct(self):

        """Return a dict in keeping with DLFIELD style"""

        dct = OrderedDict()
        atoms = []
        for atom in self.atoms:
            atoms.append(atom.to_dct())

        dct.update({"NAME" : self.name})
        dct.update({"MAXATOM" : self.nmaxatom})
        dct.update({"ATOMS" : atoms})
        dct.update({"BONDS" : self.bonds})

        dct.update({"EXCLUDE": self.exc_coul_ints})
        dct.update({"RIGID": self.rigid})

        return dct



class BondTwoBody(object):

    """2-Body bonded interaction"""

    def __init__(self, atom1, atom2, interaction):

        """Initialise content

        Arguments:
          atom1 (Atom):     one end
          atom2 (Atom):     other end
          interation   DLInteration description
        """

        self.atom1 = atom1
        self.atom2 = atom2
        self.interaction = interaction


    def __repr__(self):

        """Return internal represetation"""

        bond = "atom1= {!r}, atom2= {!r}, interaction= {!r}" \
            .format(self.atom1, self.atom2, self.interaction)

        return "BondTwoBody({!s})".format(bond)


    def __str__(self):

        """Return string in DL FIELD style"""

        at1 = "{!s} {!s}".format(self.atom1.name, self.atom1.type)
        at2 = "{!s} {!s}".format(self.atom2.name, self.atom2.type)

        return "{!s} {!s} {!s}".format(at1, at2, self.interaction)


    def to_dct(self):

        """Translate to dict in keeping with FIELD style"""

        dct = OrderedDict()
        # Could use the AtomType object here, but don't replicate mass etc
        # at the moment
        at1 = self.atom1
        at2 = self.atom2
        dct.update({"ATOM1" : at1.name})
        dct.update({"TYPE1" : at2.type})
        dct.update({"ATOM2" : at2.name})
        dct.update({"TYPE2" : at2.type})
        dct.update({"INTERACTION" : self.interaction.to_dct()})

        return dct

    @staticmethod
    def from_string(dlstr):

        """Generate Bond2Body from a FIELD file string"""

        try:
            atstr1, type1, atstr2, type2, key = dlstr.split(None, 4)

        except ValueError:
            raise ValueError("Cannot parse bond: {!r}".format(dlstr))

        at1 = Atom(atstr1, type1)
        at2 = Atom(atstr2, type2)
        interaction = dlinteraction.from_string(key)

        bond = BondTwoBody(at1, at2, interaction)

        return bond


class VDW(object):

    """A VDW (aka 2-body non-bonded van der Waals) interaction"""

    def __init__(self, atom1, atom2, interaction):

        """Initialise description

        Arguments:
          atom1 (Atom):                at one end
          atom2 (Atom):                at other end
          interaction (DLInteration):  DLInteraction description
        """

        self.atom1 = atom1
        self.atom2 = atom2
        self.interaction = interaction


    def __repr__(self):

        """Return internal representation"""

        vdw = "atom1={!r}, atom2={!r}, interaction={!r}" \
            .format(self.atom1, self.atom2, self.interaction)

        return "VDW({!s})".format(vdw)


    def __str__(self):

        """Return form appropriate for DL FIELD file"""

        atom1 = "{!s} {!s}".format(self.atom1.name, self.atom1.type)
        atom2 = "{!s} {!s}".format(self.atom2.name, self.atom2.type)

        return "{!s} {!s} {!s}".format(atom1, atom2, self.interaction)


    def to_dct(self):

        """Translate to dict in keeping with FIELD style"""

        dct = OrderedDict()
        # Could use the AtomType object here, but don't replicate mass etc
        # at the moment
        at1 = self.atom1
        at2 = self.atom2

        dct.update({"ATOM1" : at1.name})
        dct.update({"TYPE1" : at2.type})
        dct.update({"ATOM2" : at2.name})
        dct.update({"TYPE2" : at2.type})
        dct.update({"INTERACTION" : self.interaction.to_dct()})

        return dct

    @staticmethod
    def from_string(dlstr):

        """Return instance from FIELD file string"""

        try:

            # Split off the interaction key string
            atstr1, type1, atstr2, type2, key = dlstr.split(None, 4)

        except ValueError:
            raise ValueError("Could not parse VDW: {!r}".format(dlstr))

        at1 = Atom(atstr1, type1)
        at2 = Atom(atstr2, type2)
        interaction = dlinteraction.from_string(key)

        vdw = VDW(at1, at2, interaction)

        return vdw


class FIELD(object):

    """FIELD file container object"""

    def __init__(self):

        """Initialise an empty container"""

        self.description = "Field file Header"
        self.cutoff = 1.0
        self.units = "eV"
        self.nconfigs = 1
        self.atomtypes = []
        self.moltypes = []
        self.vdw = []
        self.bonds2body = []

    def __str__(self):

        """Intended to reproduce a well-formed DL FIELD file"""

        listme = []
        listme.append(self.description)
        listme.append("CUTOFF {!s}".format(self.cutoff))
        listme.append("UNITS {!s}".format(self.units))
        listme.append("NCONFIGS {!s}".format(self.nconfigs))

        listme.append("ATOMS {!s}".format(len(self.atomtypes)))
        for atom in self.atomtypes:
            listme.append(str(atom))

        listme.append("MOLTYPES {!s}".format(len(self.moltypes)))
        for mol in self.moltypes:
            listme.append(str(mol))

        listme.append("FINISH")

        listme.append("VDW {!s}".format(len(self.vdw)))
        for vdw in self.vdw:
            listme.append(str(vdw))

        if self.bonds2body:
            listme.append("BONDS {!s}".format(len(self.bonds2body)))
            for bond in self.bonds2body:
                listme.append(str(bond))

        listme.append("CLOSE")

        return "\n".join(listme)


    def __repr__(self):

        """Return a string"""

        str1 = "description={!r}, cutoff={!r}, units={!r}, nconfigs={!r}" \
            .format(self.description, self.cutoff, self.units, self.nconfigs)
        str2 = "atomtypes={}, moltypes={}, vdw={}, bonds2body={}"\
            .format(self.atomtypes, self.moltypes, self.vdw, self.bonds2body)

        return "FIELD({!s}, {!s})".format(str1, str2)


    def to_dct(self):

        """Translate to dict in keeping with FIELD style"""

        dct = OrderedDict()

        dct.update({"DESCRIPTION": self.description})
        dct.update({"CUTOFF": self.cutoff})
        dct.update({"UNITS" : self.units})
        dct.update({"NCONFIGS" : self.nconfigs})

        atomtypes = []
        for atom in self.atomtypes:
            atomtypes.append(atom.to_dct())

        dct.update({"ATOMTYPES" : atomtypes})

        moltypes = []
        for mol in self.moltypes:
            moltypes.append(mol.to_dct())

        dct.update({"MOLTYPES" : moltypes})

        vdwlist = []
        for vdw in self.vdw:
            vdwlist.append(vdw.to_dct())

        dct.update({"VDW" : vdwlist})

        b2list = []
        for bond in self.bonds2body:
            b2list.append(bond.to_dct())

        dct.update({"BONDS" : b2list})

        return dct


    def get_atomtype(self, name, atype):

        """Return AtomType member if identifed, or None"""

        for atomtype in self.atomtypes:
            if atomtype.name == name and atomtype.type == atype:
                return atomtype

        return None


    def to_json(self):

        """Return JSON representation"""

        dct = self.to_dct()
        return json.dumps(dct, indent=2)


def from_string(dlstr):

    """Generate a FIELD object from the content of FIELD file

    Arguments:
    dlstr (string):  string with blank lines and comments removed
    """

    dlfield = FIELD()

    # Format; (1,2,3 are trivial; 4 and 5 have a complex structure)
    # 1. cutoff
    # 2. units
    # 3. nconfigs
    # 4. one or more species
    # 5. one or more potentials
    # 6. "close"
    #
    # This could do with being broken up, as it's rather long...

    lines = dlstr.splitlines()

    try:
        dlfield.description = lines.pop(0)

        # cutoff
        line = lines.pop(0).split()
        if line[0].lower() != "cutoff":
            raise ValueError("Excpected 'cutoff'")
        dlfield.cutoff = float(line[1])

        # units
        line = lines.pop(0).split()
        if line[0].lower() != "units":
            raise ValueError("Expected 'units'")
        dlfield.units = str(line[1])

        # nconfigs
        line = lines.pop(0).split()
        if line[0].lower() != "nconfigs":
            raise ValueError("Expected 'nconfigs'")
        dlfield.nconfigs = int(line[1])

    except (ValueError, IndexError):
        raise ValueError("Failed to parse: {!r}".format(line))

    # Species ("atom types")

    try:
        # "atom" or "atomtypes" or "atom types"
        line = lines.pop(0).split()
        if (line[0].lower()).startswith("atom"):
            natomtypes = int(line[-1])
        else:
            raise ValueError("Expected 'atom ...'")

    except (ValueError, IndexError):
        raise ValueError("Failed to parse atom types line: {!r}".format(line))

    # Read atoms

    for _ in range(natomtypes):
        line = lines.pop(0)
        atomtype = AtomType.from_string(line)
        dlfield.atomtypes.append(atomtype)


    # "mol" or "moltypes" or "mol types" or "molspec" or "mol spec"

    line = lines.pop(0).split()
    assert line[0].lower()[0:3] == "mol"
    nmoltypes = int(line[-1])

    for _ in range(nmoltypes):
        line = lines.pop(0).split()
        assert len(line) == 1
        symbol = line[0]

        # Next line: "maxatom" or "atoms ..."
        # The logic here is:
        #   if (maxatom)
        #     read maxatom and go to next moltype
        #   else if (atoms)
        #     read atoms
        #     then read zero or more bonds, angles, etc
        #     until "finish" encountered upon which next moltype

        line = lines.pop(0).split()

        if line[0].lower() == "maxatom":
            nmaxatom = int(line[-1])
            mol = MolType(symbol, nmaxatom)
            dlfield.moltypes.append(mol)

        elif line[0].lower() == "atoms":

            # Process "atoms..."

            assert len(line) == 3
            natom = int(line[1])
            nmaxatom = int(line[2])
            assert natom <= nmaxatom
            mol = MolType(symbol, nmaxatom)

            for _ in range(natom):
                nextline = lines.pop(0).split()
                assert len(nextline) == 5
                name = nextline[0]
                atype = nextline[1]

                # Make sure this entry is in the list of atomtypes
                awanted = None
                for a in dlfield.atomtypes:
                    if a.name == name and a.type == atype:
                        awanted = a

                assert awanted is not None

                x = float(nextline[2])
                y = float(nextline[3])
                z = float(nextline[4])
                atom = Atom(name, atype, [x, y, z])
                mol.atoms.append(atom)

            # Still with this mol type, look for bonds angles ...

            while True:
                nextline = lines.pop(0).split()

                if nextline[0].lower() == "finish":
                    dlfield.moltypes.append(mol)
                    # Next moltype loop iteration
                    break

                # Optional "exclude" key
                if nextline[0].lower()[0:5] == "exclu":
                    mol.exc_coul_ints = True

                # Optional "rigid"
                if nextline[0].lower() == "rigid":
                    mol.rigid = True

                if nextline[0].lower() == "bonds":
                    # Have "bonds number"
                    assert len(nextline) == 2

                    nbonds = int(line[1])
                    # Bond sepcification:
                    # "atomindex1 atomindex2 bondtypeindex"

                    for _ in range(nbonds):
                        assert len(lines[0].split()) == 3
                        ia1, ia2, ibond = map(int, lines.pop(0).split())
                        mol.bonds.append([ia1, ia2, ibond])


    # Potentials

    closed = False

    while not closed:
        line = lines.pop(0).split()

        if line[0].lower() == "close":
            # discard any lines after close statement
            closed = True
            break

        key = line[0].lower()

        # "vdw" or "vdw types"
        # "bond" or "bonds" or "bond types" etc

        if key == "vdw":
            nvdw = int(line[-1])

            for _ in range(nvdw):
                line = lines.pop(0)
                vdw = VDW.from_string(line)
                dlfield.vdw.append(vdw)

        elif key.startswith("bond"):
            nbond = int(line[-1])

            for _ in range(nbond):
                line = lines.pop(0)
                bond = BondTwoBody.from_string(line)
                dlfield.bonds2body.append(bond)

        elif key == "finish":
            # There is actually a possibility that the "finish"
            # line from the moltypes section has not been
            # consumed following a "maxtaom" key. It must
            # then be discarded here.
            pass

        else:
            # Haven't considered other potential types!
            raise NotImplementedError(line)


    return dlfield


def from_file(filename):

    """Load FIELD ASCII file via full pathname and return FIELD object"""

    lines = dlutil.load_ascii(filename)
    field = from_string("\n".join(lines))

    return field
