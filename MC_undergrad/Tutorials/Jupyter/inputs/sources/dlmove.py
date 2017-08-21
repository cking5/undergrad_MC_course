"""Containers for DL CONTROL file MC move type descriptions

Moves are part of the DL CONTROL file input. Each type of
move gets a class here.

The classification of the available Move types is as follows:

Move
  MCMove
    AtomMove
    MoleculeMove
    ... [others, some of which are untested in regression tests]
  VolumeMove (aka NPT move)
    VolumeVectorMove
    VolumeOrthoMove
    VolumeCubicMove

The dlmove.from_string(dlstr) function provides a factory method to
generate the appropriate Move from DL CONTROL style input.

"""


def parse_atom(dlstr):

    """Atoms are expected as 'Name core|shell' only"""

    try:
        tokens = dlstr.split()
        atom = {"id": "{} {}".format(tokens[0], tokens[1])}
    except IndexError:
        raise ValueError("Unrecognised Atom: {!s}".format(dlstr))

    return atom

def print_atom(atom):

    """Return atom record for DL CONTROL"""

    return atom["id"]

def parse_molecule(dlstr):

    """Molecules are 'Name" only"""

    tokens = dlstr.split()
    molecule = {"id": tokens[0]}

    return molecule

def print_molecule(molecule):

    """Return molecule record for DL CONTROL"""

    return molecule["id"]

def parse_atom_swap(dlstr):

    """Swaps are e.g., 'atom1 core atom2 core' """

    try:
        tok = dlstr.split()
        swap = {"id1": "{} {}".format(tok[0], tok[1]), \
                    "id2": "{} {}".format(tok[2], tok[3])}
    except IndexError:
        raise ValueError("Unrecognised atom swap: {!r}".format(dlstr))

    return swap


def print_atom_swap(swap):

    """Return atom swap string for DL CONTROL"""

    return "{} {}".format(swap["id1"], swap["id2"])


def parse_molecule_swap(dlstr):

    """Swap records have two tokens"""

    try:
        tokens = dlstr.split()
        swap = {"id1": tokens[0], "id2": tokens[1]}
    except IndexError:
        raise ValueError("Unrecognised molecule swap: {!r}".format(dlstr))

    return swap


def print_molecule_swap(swap):

    """Return a swap for DL CONTROL output"""

    return "{} {}".format(swap["id1"], swap["id2"])


def parse_atom_gcmc(dlstr):

    """GCMC Atoms include a chemical potential/partial pressure"""

    try:
        tok = dlstr.split()
        atom = {"id": "{} {}".format(tok[0], tok[1]), "molpot": float(tok[2])}
    except (IndexError, TypeError):
        raise ValueError("Unrecognised GCMC Atom: {!r}".format(dlstr))

    return atom


def parse_molecule_gcmc(dlstr):

    """Grand Canonical MC includes chemical potential/partial pressure"""

    try:
        tok = dlstr.split()
        molecule = {"id": tok[0], "molpot": float(tok[1])}
    except (IndexError, TypeError):
        raise ValueError("Unrecognised GCMC Molecule: {!r}".format(dlstr))

    return molecule


def print_gcmc(gcmc):

    """Return string version of chemical potential records for DL CONTROL"""

    return "{} {}".format(gcmc["id"], gcmc["molpot"])


class Move(object):

    """This classifies all moves"""

    key = None

    @classmethod
    def from_string(cls, dlstr):

        """To be implemented by MCMove or VolumeMove"""

        NotImplementedError("Should be implemented by subclass")



class MCMove(Move):

    """MC moves involve atoms or molecules"""

    parse_mover = staticmethod(None)
    print_mover = staticmethod(None)

    def __init__(self, pfreq, movers):

        """pfreq (int):      percentage probability of move per step"""

        self.pfreq = pfreq
        self.movers = movers

    def __str__(self):

        """Return well-formed DL CONTROL block"""

        strme = []
        move = "move {} {} {}".format(self.key, len(self.movers), self.pfreq)
        strme.append(move)

        for mover in self.movers:
            strme.append(self.print_mover(mover))

        return "\n".join(strme)


    def __repr__(self):

        """Return a readable form"""

        repme = "pfreq= {!r}, movers= {!r}".format(self.pfreq, self.movers)

        return "{}({})".format(type(self).__name__, repme)


    @classmethod
    def from_string(cls, dlstr):

        """Generate an instance from a DL CONTROL entry"""

        lines = dlstr.splitlines()
        line = lines.pop(0)
        pfreq = MCMove._parse_move_statement(line)[2]

        movers = []
        for line in lines:
            mover = cls.parse_mover(line)
            movers.append(mover)

        return cls(pfreq, movers)


    @staticmethod
    def _parse_move_statement(dlstr):

        """Parse move line"""

        try:
            tokens = dlstr.lower().split()
            if tokens[0] != "move":
                raise ValueError("Expected 'move' statement")

            mtype, nmove, pfreq = tokens[1], int(tokens[2]), int(tokens[3])
        except IndexError:
            raise ValueError("Badly formed 'move' statement?")

        return mtype, nmove, pfreq


class GCMove(Move):

    """Grand Canonical insertions have an extra minimum insertion
    distance parameter cf standard MCMove types
    """

    parse_mover = staticmethod(None)
    print_mover = staticmethod(None)

    def __init__(self, pfreq, rmin, movers):

        """Initalise GCMove parameters

        Arguments:
            pfreq (integer):    percentage
            rmin (float):       grace distance
            movers ():
        """

        self.pfreq = pfreq
        self.rmin = rmin
        self.movers = movers


    def __str__(self):

        """Return well-formed DL CONTROL block"""

        strme = []
        move = "move {} {} {} {}".format(self.key, len(self.movers),
                                         self.pfreq, self.rmin)
        strme.append(move)

        for mover in self.movers:
            strme.append(self.print_mover(mover))

        return "\n".join(strme)


    def __repr__(self):

        """Return a GCMove (subsclass) represetation"""

        repme = "pfreq= {!r}, rmin= {!r}, movers= {!r}"\
            .format(self.pfreq, self.rmin, self.movers)

        return "{}({})".format(type(self).__name__, repme)


    @classmethod
    def from_string(cls, dlstr):

        """Generate instance form well-formed DL CONTROL string"""

        lines = dlstr.splitlines()
        line = lines.pop(0)
        pfreq, rmin = GCMove._parse_move_statement(line)[2:]

        movers = []
        for line in lines:
            mover = cls.parse_mover(line)
            movers.append(mover)

        return cls(pfreq, rmin, movers)


    @staticmethod
    def _parse_move_statement(dlstr):

        """Parse GC move line"""

        try:
            tokens = dlstr.lower().split()
            if tokens[0] != "move":
                raise ValueError("Expected 'move' statement")

            mtype, nmove, pfreq, rmin = \
                tokens[1], int(tokens[2]), int(tokens[3]), float(tokens[4])
        except IndexError:
            raise ValueError("Badly formed 'move' statement?")

        return mtype, nmove, pfreq, rmin



class VolumeMove(Move):

    """Container for volume (NPT) moves"""

    def __init__(self, pfreq, sampling=None):

        """Initialise continaer

        Arguemnts:
            pfreq (integer):      percentage
            sampling (string):    description
        """

        self.pfreq = pfreq
        self.sampling = sampling


    def __str__(self):

        """Return well-formed DL CONTROL file string"""

        if self.sampling is not None:
            strme = "move volume {} {} {}"\
                .format(self.key, self.sampling, self.pfreq)
        else:
            strme = "move volume {} {}".format(self.key, self.pfreq)

        return strme


    def __repr__(self):

        """Return a readable string"""

        repme = "pfreq= {!r}, sampling= {!r}".format(self.pfreq, self.sampling)

        return "{}({})".format(type(self).__name__, repme)


    @classmethod
    def from_string(cls, dlstr):

        """E.g., 'move volume vector|ortho|cubic [sampling-type] pfreq' """

        tokens = dlstr.split()

        try:
            sampling = None
            pfreq = int(tokens[-1])
            # sampling-type is an optional one or two (string) tokens

            if len(tokens) == 5:
                sampling = tokens[3]
            if len(tokens) == 6:
                sampling = "{} {}".format(tokens[3], tokens[4])
        except (IndexError, TypeError):
            raise ValueError("VolumeMove: unrecognised: {!r}".format(dlstr))

        return cls(pfreq, sampling)



class AtomMove(MCMove):

    """Concrete class for atom moves"""

    key = "atom"
    parse_mover = staticmethod(parse_atom)
    print_mover = staticmethod(print_atom)

class MoleculeMove(MCMove):

    """Concrete class for molecule moves"""

    key = "molecule"
    parse_mover = staticmethod(parse_molecule)
    print_mover = staticmethod(print_molecule)

class RotateMoleculeMove(MCMove):

    """Concrete class for rotate molecule moves"""

    key = "rotatemol"
    parse_mover = staticmethod(parse_molecule)
    print_mover = staticmethod(print_molecule)

class SwapAtomMove(MCMove):

    """Concrete class for swap atom moves"""

    key = "swapatoms"
    parse_mover = staticmethod(parse_atom_swap)
    print_mover = staticmethod(print_atom_swap)

class SwapMoleculeMove(MCMove):

    """Concrete class for swap molecule moves"""

    key = "swapmols"
    parse_mover = staticmethod(parse_molecule_swap)
    print_mover = staticmethod(print_molecule_swap)

class InsertAtomMove(GCMove):

    """Concrete class for Grand Canonical atom moves"""

    key = "gcinsertatom"
    parse_mover = staticmethod(parse_atom_gcmc)
    print_mover = staticmethod(print_gcmc)

class InsertMoleculeMove(GCMove):

    """Concrete class for Grand Canonical molecule moves"""

    key = "gcinsertmol"
    parse_mover = staticmethod(parse_molecule_gcmc)
    print_mover = staticmethod(print_gcmc)

class SemiWidomAtomMove(MCMove):

    """No exmaples are available"""

    key = "semiwidomatoms"
    # Format needs to be confirmed

class SemiGrandAtomMove(MCMove):

    """No examples are available"""

    key = "semigrandatoms"
    # Format needs to be confirmed

class SemiGrandMoleculeMove(MCMove):

    """NO examples are available"""

    key = "semigrandmol"
    # Format needs to be confirmed


class VolumeVectorMove(VolumeMove):

    """Concrete class for vector volume moves"""

    key = "vector"

class VolumeOrthoMove(VolumeMove):

    """Concrete class for ortho volume moves"""

    key = "ortho"

class VolumeCubicMove(VolumeMove):

    """Concrete class for cubic volume moves"""

    key = "cubic"


def from_string(dlstr):

    """Factory method to return an instance from a well-formed
    DL CONTROL file move statement (a block of 1 plus n lines)

    Argument:
    dlstr (string)          move statement plus atom/molecule
                            descriptions
    """

    moves_volume = {"vector": VolumeVectorMove,
                    "ortho": VolumeOrthoMove,
                    "cubic": VolumeCubicMove}

    moves_mc = {"atom": AtomMove,
                "molecule": MoleculeMove,
                "rotatemol": RotateMoleculeMove,
                "gcinsertatom": InsertAtomMove,
                "gcinsertmol": InsertMoleculeMove}

    lines = dlstr.splitlines()
    tokens = lines[0].lower().split()
    if tokens[0] != "move" or len(tokens) < 4:
        raise ValueError("Expected: 'move key ...': got {!r}".format(lines[0]))

    key = tokens[1]

    # We need to allow for possible DL key abbreviations
    if key.startswith("atom"):
        key = "atom"
    if key.startswith("molecu"):
        key = "molecule"
    if key.startswith("rotatemol"):
        key = "rotatemol"

    inst = None
    if key == "volume":
        subkey = tokens[2]
        if subkey in moves_volume:
            inst = moves_volume[subkey].from_string(dlstr)
    else:
        if key in moves_mc:
            inst = moves_mc[key].from_string(dlstr)

    if inst is None:
        raise ValueError("Move unrecognised: {!r}".format(dlstr))

    return inst
