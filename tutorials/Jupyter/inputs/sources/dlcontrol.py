"""A representation of the DL-MONTE CONTROL file

The overall intention of the CONTROL class is to be able
to read in the contents of te CONTROL file, manipulate the
contents as required, as write out again a file suitable for
use by the main DL-MONTE code.

The module level function from_file() can be used to generate
an instance of CONTROL, e.g.,

ctrl = dlcontrol.from_file(file_name)


The CONTROL file has a potentially complex structure which is
here broken dewn into smaller structured parts.

The structure is as follows. Each of the subblocks has a
corresponding container class (except the title, which is
just a string):

- title:
- UseBlock:
  - use_statements: {}
  - fed_block:
    - FEDFlavour:
    - FEDMethod:
    - FEDOrderParameter:
- MainBlock
  - statements: {}
  - moves: []
  - samples: {}

"""

from collections import OrderedDict

import inputs.sources.dlutil as dlutil
import inputs.sources.dlmove as dlmove
import inputs.sources.dlfedflavour as fedflavour
import inputs.sources.dlfedmethod as fedmethod
import inputs.sources.dlfedorder as fedorder

class FEDBlock(object):

    """FED Container object for CONTROL input.

    Attributes:
    flavour:     FEDFlavour description
    method:      FEDMethod description
    orderparam:  FEDOrderParameter description

    See, e.g.,  DL-MONTE source fed_interface_type.f90
    """

    def __init__(self, flavour=None, method=None, orderparam=None):

        """Initialise an FED container"""

        self.flavour = flavour
        self.method = method
        self.orderparam = orderparam


    def __str__(self):

        """Returns a well-formed DL CONTROL FED file block"""

        listme = [str(self.flavour), str(self.method), str(self.orderparam)]
        listme.append("fed done")

        return "\n".join(listme)


    def __repr__(self):

        """Returns full information"""

        repme = "flavour= {!r}, method= {!r}, orderparam= {!r}"\
            .format(self.flavour, self.method, self.orderparam)

        return "FEDBlock({})".format(repme)


    @classmethod
    def from_string(cls, dlstr):

        """Parse the fed section of the CONTROL file
        """

        # The "use fed" is a single line, or a block ("ps")...
        # Some phaff is required to form the right single string
        # for the "ps" case ...

        # Likewise "fed order ..." is a single line, except in the
        # case of "fed order com2..." where it is a structured block

        flavour, method, orderparameter = None, None, None

        try:
            lines = dlstr.splitlines()
            line = lines[0].lower()
            tokens = line.split()

            if tokens[0] != "use" or tokens[1] != "fed":
                raise ValueError()

            if tokens[2] not in fedflavour.FLAVOURS:
                raise ValueError()

        except (ValueError, IndexError):
            msg = "Expecting 'use fed' line or block; got {!r}".format(line)
            raise ValueError(msg)

        if tokens[2] == "gen" or tokens[2] == "generic":
            flavour = fedflavour.Generic.from_string(line)
        elif tokens[2] == "ps":
            psblock = FEDBlock._join_pslist(lines)
            flavour = fedflavour.PhaseSwitch.from_string(psblock)
        else:
            # Something in fedflavor.flavours is not implemented
            raise NotImplementedError(line)

        lines.pop(0)

        # Process "fed method" lines and "fed order parameter" section
        # until "fed done" is encountered

        while True:
            # If lines are exhausted, it suggests a badly formed file

            line = lines[0].lower()
            tokens = line.split()

            if tokens[0] == "fed" and tokens[1] == "done":
                break
            elif tokens[1] == "method":
                method = fedmethod.from_string(line)
            elif tokens[1] == "order":
                # Watch out for optional "param" at position 2
                if tokens[2] == "com2" or tokens[3] == "com2":
                    opblock = FEDBlock._join_opblock(lines)
                else:
                    opblock = line
                orderparameter = fedorder.from_string(opblock)
            else:
                raise ValueError(line)

            lines.pop(0)

        return FEDBlock(flavour, method, orderparameter)


    @staticmethod
    def _join_pslist(lines):

        """Pull out a single string with a structured  ps block"""

        psblock = []

        # Leave lines[0] ("fed ps ...") in the list
        psblock.append(lines[0].lower())

        try:
            done = False
            while not done:
                line = lines.pop(1).lower()
                tokens = line.split()
                psblock.append(line)
                if tokens[0] == "ps" and tokens[1] == "done":
                    done = True
        except Exception:
            raise ValueError("Badly formed fed ps block")

        return "\n".join(psblock)


    @staticmethod
    def _join_opblock(lines):

        """Pull out a single string with a structured  order parameter block"""

        psblock = []
        psblock.append(lines[0].lower())

        try:
            done = False
            while not done:
                line = lines.pop(1).lower()
                tokens = line.split()
                psblock.append(line)
                if tokens[0] == "fed" and tokens[1] == "order" and tokens[-1] == "done":
                    done = True
        except Exception:
            raise ValueError("Badly formed fed order parameter block")

        return "\n".join(psblock)


class UseBlock(object):

    """Container for the structured 'use' block in the CONTROL file"""

    # Acceptable keys (for the CONTROL input) and the associated
    # internal DL variables (as values, just for reference).
    # "use fed" is excluded here, as it is dealt with separately.

    # The list is used to check the input; this is ok as long as
    # no abbreviations are allowed.

    # Some keys (repexch) have arguments requiring further
    # processing.

    _use_keys = ["gaspressure",
                 "repexch",
                 "rotquaternion",
                 "seqmove",
                 "seqmolmove",
                 "seqmolrot"]


    def __init__(self, use_statements=None, fed_block=None):

        """Initialise a CONTREOL file 'use block'

        Arguemnts:
        use_statements (OrderedDict): list of active use keywords (strings)
        fed_block (FEDBlock):         fed container (may be None)
        """

        if use_statements is None:
            self.use_statements = OrderedDict()
        else:
            self.use_statements = use_statements

        self.fed_block = fed_block


    def __str__(self):

        """Return well formed DL CONTROL use block"""

        lines = []
        for key in self.use_statements:
            if self.use_statements[key] is None:
                lines.append("use {}".format(key))
            else:
                values = [v for _, v in self.use_statements[key].items()]
                values = " ".join(map(str, values))
                lines.append("use {} {}".format(key, values))

        if self.fed_block is not None:
            lines.append(str(self.fed_block))

        lines.append("finish use-block")

        return "\n".join(lines)


    def __repr__(self):

        """Return full information"""

        repme = "use_statements= {!r}, fed_block= {!r}"\
            .format(self.use_statements, self.fed_block)

        return "UseBlock({})".format(repme)


    @classmethod
    def from_string(cls, dlstr):

        """Generate and return an instance from DL CONTROL block

        dlstr (string):    a new-line separated string with
                           comments and blank lines removed.

        The use block must be terminated with the line
        'finish [use-block]'
        """

        # Extract use block and [optional] fed block

        useblock, fedblock = UseBlock._split_blocks(dlstr)

        use_statements = UseBlock._use_from_string(useblock)
        fed = None

        if fedblock is not None:
            fed = FEDBlock.from_string(fedblock)

        return cls(use_statements, fed)


    @staticmethod
    def _use_from_string(dlstr):

        """Parse use statements (with any fed content removed)"""

        statements = OrderedDict()
        lines = dlstr.lower().split("\n")

        try:
            done = False
            for line in lines:
                tokens = line.split()

                if tokens[0] == "finish":
                    done = True
                    break

                key = tokens[1]
                if key in UseBlock._use_keys:
                    if key == "repexch":
                        nr = int(tokens[2])
                        dt = float(tokens[3])
                        n = int(tokens[4])
                        args = {"nr": nr, "deltat": dt, "nstep": n}
                        statements.update({key: args})
                    else:
                        statements.update({key: None})
                else:
                    msg = "Unrecognised use statement:"
                    raise ValueError("{}: {!r}".format(msg, line))

            if not done:
                raise ValueError("End of use block without 'finish'?")

        except IndexError:
            msg = "Failed to parse use block at: "
            raise ValueError("{} {!r}".format(msg, line))

        return statements


    @staticmethod
    def _split_blocks(dlstr):

        """Split use statements and structured fed block (if present)"""

        # The use block will contain at least a "finish"

        useblock = ""
        fedblock = None

        lines = dlstr.lower().split("\n")
        infedblock = False

        try:
            for line in lines:
                tokens = line.split()

                if tokens[0] == "use" and tokens[1] == "fed":
                    fedblock = ""
                    infedblock = True

                if infedblock:
                    fedblock += line + "\n"
                else:
                    useblock += line + "\n"

                if tokens[0] == "fed" and tokens[1] == "done":
                    infedblock = False

                if tokens[0] == "finish" and infedblock:
                    raise ValueError("'finish' in use fed block?")

        except IndexError:
            msg = "Could not split fed and use blocks"
            raise ValueError(msg)

        return useblock, fedblock


class MainBlock(object):

    """Container for CONTROL 'main' block information"""

    def __init__(self, statements=None, moves=None, samples=None):

        self.statements = OrderedDict()
        self.moves = []
        self.samples = OrderedDict()

        if statements is not None:
            self.statements = statements
        if moves is not None:
            self.moves = moves
        if samples is not None:
            self.samples = samples

    def __str__(self):

        """Return a string with a well-formed main section"""

        block = []
        for key in self.statements:
            try:
                values = [v for k, v in self.statements[key].items()]
                values = " ".join(map(str, values))
            except AttributeError:
                values = self.statements[key]
            block.append("{} {}".format(key, values))

        for sample in self.samples:
            values = [v for k, v in self.samples[sample].items()]
            values = " ".join(map(str, values))
            block.append("sample {} {}".format(sample, values))

        for move in self.moves:
            block.append(str(move))

        block.append("start simulation")

        return "\n".join(block)

    def __repr__(self):

        """Return a MainBlock representation"""

        block = "statements= {!r}, moves= {!r}, samples= {!r}"\
            .format(self.statements, self.moves, self.samples)

        return "MainBlock({})".format(block)


    @classmethod
    def from_string(cls, dlstr):

        """Parse 'main' block of CONTROL input"""

        statements = OrderedDict()
        moves = []
        samples = OrderedDict()

        lines = dlstr.splitlines()

        while True:

            line = lines[0]
            tokens = line.lower().split()
            key = tokens[0]

            # A number of special cases:
            #   "start" denotes the "end of input"
            #   various move types are possible...
            #   various sample types are possible...

            if key == "start":
                break
            elif key[0:4] == "move":
                move = MainBlock._move_from_string(lines)
                moves.append(move)
            elif key == "sample":
                sample = MainBlock._sample_from_string(line)
                samples.update(sample)
            else:
                # Other single line statements
                item = MainBlock._parse_single_line_statement(line)
                statements.update(item)

            lines.pop(0)

        return cls(statements, moves, samples)


    @staticmethod
    def _parse_single_line_statement(dlstr):

        """Return a dict with a keyword and appropirate values"""

        # This is rather ugly, as every key value is potentially
        # a special case requiring particular action.

        # However, there are a number of cases which are a
        # simple key value pair where no abbreviations are
        # possible and parsing can occur via...

        _simple_key_values = {"acceptatmmoveupdate": float,
                              "acceptmolmoveupdate": int,
                              "acceptmolrotupdate": int,
                              "check": int,
                              "maxatmdist": float,
                              "maxmolrot": float,
                              "maxnonbondnbrs": int,
                              "print": int,
                              "revconformat": str,
                              "stack": int,
                              "steps": int,
                              "yamldata": int}

        try:
            tokens = dlstr.lower().split()
            key = tokens[0]

            # These 'if' cases are in alphabetical order. The temporary
            # variable names correspond to the DL Fortran variable.

            # The final token[-1] is used to get a value in many cases
            # to avoid parsing optional intermediary key words.

            if key in _simple_key_values:
                value = _simple_key_values[key](tokens[1])
                item = {key: value}
            elif key == "distewald":
                item = {key: None}
            elif key.startswith("equil"):
                item = {key: int(tokens[-1])}
            elif key.startswith("ewald") and tokens[1].startswith("prec"):
                item = {"ewald precision": float(tokens[2])}
            elif key.startswith("ewald") and  tokens[1].startswith("sum"):
                # arguments are real int int int
                alpha = float(tokens[2])
                kmax1 = int(tokens[3])
                kmax2 = int(tokens[4])
                kmax3 = int(tokens[5])
                ewald = OrderedDict([("alpha", alpha), ("kmax1", kmax1), \
                                         ("kmax2", kmax2), ("kmax3", kmax3)])
                item = {"ewald sum": ewald}
            elif key == "nbrlist":
                # argument is string "auto" or int
                if tokens[1] == "auto":
                    item = {key: "auto"}
                else:
                    item = {key: int(tokens[-1])}
            elif key == "noewald":
                # argument is "all" or integer box number
                if tokens[1] == "all":
                    item = {key: "all"}
                else:
                    item = {key: int(tokens[1])}
            elif key.startswith("pres"):
                item = {"pressure": float(tokens[-1])}
            elif key == "ranseed":
                item = {key: None}
            elif key == "seeds":
                # four integers
                keys = ["seed0", "seed1", "seed2", "seed3"]
                vals = list(map(int, tokens[1:5]))
                seeds = OrderedDict(zip(keys, vals))
                item = {key: seeds}
            elif key.startswith("stat"):
                item = {key: int(tokens[-1])}
            elif key.startswith("temp"):
                item = {key: float(tokens[-1])}
            elif key.startswith("toler"):
                item = {key: float(tokens[-1])}
            else:
                # This may be a bad statement, or it may be a valid
                # statement not handled above ...
                raise NotImplementedError()

        except (IndexError, ValueError):
            msg = "MainBlock: error parsing line: {!r}".format(dlstr)
            raise ValueError(msg)
        except NotImplementedError:
            msg = "MainBlock: unrecognised statement: {!r}".format(dlstr)
            raise ValueError(msg)

        return item


    @staticmethod
    def _move_from_string(lines):

        """Generate a move object from DL CONTROL input

        Args:
        lines(list): the entire move command from the CONTROL file
                     as a list of lines
        """

        # This is slightly ugly. We need to preserve the "move"
        # statement itself in the input list, but we remove
        # the rest of the block (depending on how many lines
        # of input are expected).

        # The block is re-formed as a single string to generate
        # an object from the dlmove factory method

        try:
            line = lines[0]
            tokens = line.lower().split()
            key = tokens[1]

            move_block = []
            move_block.append(line)

            if tokens[1] == "volume":
                nlines = 0
            else:
                nlines = int(tokens[2])

            for _ in range(nlines):
                move_block.append(lines.pop(1))

        except (IndexError, TypeError):
            raise ValueError("Could not parse move block")


        return dlmove.from_string("\n".join(move_block))


    @staticmethod
    def _sample_from_string(dlstr):

        linedict = OrderedDict()

        try:

            tokens = dlstr.lower().split()
            key = tokens[1]

            if key.startswith("coord"):
                name = "coords"
                linedict.update({"nfreq": int(tokens[2])})
                # optional only <molecules>
            elif key.startswith("ener"):
                name = "energies"
                linedict.update({"nfreq": int(tokens[2])})
            elif key.startswith("rdf"):
                name = "rdfs"
                linedict.update({"ngrid": int(tokens[2])})
                linedict.update({"rcut": float(tokens[3])})
                linedict.update({"nfreq": int(tokens[4])})
            elif key.startswith("zden"):
                name = "zdensity"
                linedict.update({"ngrid": int(tokens[2])})
                linedict.update({"nfreq": int(tokens[3])})
            elif key.startswith("vol"):
                name = "volume"
                linedict.update({"mfreq": int(tokens[2])})
                linedict.update({"deltav": float(tokens[3])})
            else:
                raise NotImplementedError("Unrecognised sample: {!r}"\
                                              .format(dlstr))

        except (IndexError, ValueError):
            raise ValueError("Cannot parse sample: {!r}".format(dlstr))

        return {name: linedict}



class CONTROL(object):

    """CONTROL file description

    The full contents of the CONTROL structure are rather unwieldy.
    See DL-MONTE control_type.f90

    Broadly, the file has the following structure:

      Title (one-line comment)

      Zero or more "use directives"
      finish use-block

      Zero or more "main directives"
      start simulation
      EOF
    """

    def __init__(self, title="", use_block=None, main_block=None):

        """Initialise:

        title (string):          descriptive title (one line)
        use_block(UseBlock):     use statements and fed
        main_block (MainBlock):  main section
        """

        self.title = title
        self.use_block = use_block
        self.main_block = main_block


    def __str__(self):

        """Return a well-formed DL-CONTROL file as a string"""

        listme = []

        listme.append(self.title)
        listme.append(str(self.use_block))
        listme.append(str(self.main_block))

        return "\n".join(listme)


    def __repr__(self):

        """Return an internal represetation"""

        repme = "title= {!r}, use_block= {!r}, main_block= {!r}"\
            .format(self.title, self.use_block, self.main_block)

        return "CONTROL({})".format(repme)


    def ensemble(self):

        """Return the MC ensemble by interrogating the moves

        Returns a string ['nvt' | 'npt' | 'muvt']
        """

        # The logic is currently as follows:
        # 1. Assume NVT
        # 2. If a volume move is present => NPT
        # 3. If a GCMove move is present => muVT
        #
        # No examples of NVE are available at the time of writing
        # but could be identified by lack of temperature in main block(?)

        ensemble = "nvt"
        for move in self.main_block.moves:
            if isinstance(move, dlmove.VolumeMove):
                ensemble = "npt"
            if isinstance(move, dlmove.GCMove):
                ensemble = "muvt"

        return ensemble


    @classmethod
    def from_file(cls, filename="CONTROL"):

        """Return an instance from file"""

        lines = dlutil.load_ascii(filename)
        dlstr = "\n".join(lines)

        return cls.from_string(dlstr)


    @classmethod
    def from_string(cls, dlstr):

        """Return an instance from string representation

        Arguments:

        dlstr (string):   newline separateed string with comments
        and blank lines removed.

        Returns:
        a newly created instance

        Exceptions:
        A poorly formed string will raise a ValueError
        """

        # Split title, use block and main block

        title, use_str, main_str = CONTROL._split_input(dlstr)

        use_block = UseBlock.from_string(use_str)
        main_block = MainBlock.from_string(main_str)

        return cls(title, use_block, main_block)


    @staticmethod
    def _split_input(dlstr):

        use_lines = []
        main_lines = []

        try:
            lines = dlstr.splitlines()
            title = lines.pop(0)

            done = False
            inuseblock = True

            for line in lines:
                tokens = line.lower().split()

                if inuseblock:
                    use_lines.append(line)
                else:
                    main_lines.append(line)

                # End of input is "start [simulation]"
                if tokens[0] == "start":
                    done = True
                    break

                # Check for end of use block "finish [use-block]"
                if tokens[0] == "finish":
                    inuseblock = False

            if not done:
                raise ValueError("ControlInput: no 'start' statement")
            if inuseblock:
                raise ValueError("ControlInput: no 'finish [use-block]'")

        except IndexError:
            raise ValueError("Unexpected end of input")
        except Exception:
            raise

        use_str = "\n".join(use_lines)
        main_str = "\n".join(main_lines)

        return title, use_str, main_str


def from_file(filename="CONTROL"):

    """Module level factory method for CONTTROL object"""

    return CONTROL.from_file(filename)

