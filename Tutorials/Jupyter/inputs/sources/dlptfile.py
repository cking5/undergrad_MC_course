"""Utility to handle DL-MONTE PTFILE output

The PTFILE class holds the data (in a form based on the yaml
representation) and provides a number of methods to access
various time series.

A module level load() function is provided to read standard
PTFILE.000 file, an return a PTFILE object.

Another module level function load_yaml() reads the equivalent
YAMLDATA.000, and returns a YAMLDATA (sub) class object.

The format expected of the PTFILE is:

  1: [timestep]
  2: [Total  Enphalty RecipSpace RealSpace Nonbonded Threebody]
  3: [Pair   Angles   Fourbody   Manybody  External  Volume   ]
  4: [Virial L1       L2         L3        Lcos1     Lcos2    ]
  5: [Lcos3  Total^2  Enthalpy^2                              ]
  6: [number of atoms for 1...number of atom types]
  7: [number of molecules for 1 .. number of molcule types]
  8: [Per molecule type stats cf line 2]
  9: [Per molecule type stats cf line 3]
 10: [Per molecule type stats cf line 4]
 11: [Per molecule type stats cf line 5]

The last four lines are repeated for each molecule type, so
7 + 4*(number of molecule types) lines are expected per time
step.

Per molecule values (lines 8-11) are not stored at the moment.

YAML
Each YAMLDATA file contains two yaml documents:
1. a metadata doucment (no analogue in the standard PTFILE),
2. a per timestep data document (cf. the PTFILE).

"""

from collections import OrderedDict

import os
import yaml

import numpy

# Order is important in the following, hence tuples
# and OrderedDict()

ENERGY = (("energy", "Total Energy"),
          ("enthalpytot", "Enthalpy"),
          ("energyrcp", "Ewald Reciprocal Space Energy"),
          ("energyreal", "Ewald Real Space Energy"),
          ("energyvdw", "Non-bonded Energy"),
          ("energythree", "Three-body Energy"),
          ("energypair", "Pair Energy"),
          ("energyang", "Angle Energy"),
          ("energyfour", "Four-body Energy"),
          ("energymany", "Many-body Energy"),
          ("energyext", "External Energy"))

OTHERS = (("volume", "Volume"),
          ("virial", "Virial"),
          ("L1", "Lattice Vector (a)"),
          ("L2", "Lattice Vector (b)"),
          ("L3", "Lattice Vector (c)"),
          ("Lcos1", "Lattice angle (lcos1)"),
          ("Lcos2", "Lattice angle (lcos2)"),
          ("Lcos3", "Lattice amgle (lcos3)"),
          ("energytotsq", "Total Energy (Squared)"),
          ("enthalpytotsq", "Enthalpy (Squared)"))

KEY_NATOM = "natom"
KEY_NMOL = "nmol"
KEY_TIMESTEP = "timestamp"
KEYS = OrderedDict(ENERGY + OTHERS)

"""PTFILE_KEYS describes the order in which the various quantities
appear in the file in terms of the keys used in the YAML represenation
"""

class PTFILE(object):

    """PTFILE data container"""

    def __init__(self, data):

        """Initialise the yaml-like data"""

        self.data = data

    @property
    def keys(self):

        """Return a current list of keys in the data"""

        try:
            keys = self.data[0].keys()
        except (IndexError, TypeError):
            keys = []

        return keys

    def time_steps(self):

        """Return a time series of time steps

        Returns:
        data (list of int):    the timestep values appearing in PTFILE
        """

        data = []
        for step in self.data:
            data.append(step[KEY_TIMESTEP])

        return data


    def time_series(self, key):

        """Return a time series for given key

        Arguments:
        key (string):          identify what's wanted (from KEYS)

        Returns:
        data (list of float):  time series for quantity "key"
        """

        data = []
        for step in self.data:
            data.append(step[key])

        return data


    def natom(self):

        """ Return the number of atom types, and a set of time series;
        one for each atom type

        Returns:
        natamtypes (integer):  number of natomtypes
        data (list of list):   time series for each atom type
        """

        # Internal data need to be 'transposed' to give a list
        # of time series appropriate for each atom type

        natomtypes = len(self.data[0][KEY_NATOM])
        data = [[] for atom in range(natomtypes)]

        for step in self.data:
            for atom in range(natomtypes):
                data[atom].append(step[KEY_NATOM][atom])

        return natomtypes, data


    def nmol(self):

        """Return the number of molecule types, and a set of time series;
        one for each molecule type

        Returns:
        nmoltypes (integer):   number of molecule types
        data (list of list):   time series for each molecule type
        """

        # These data need to be 'transposed' to give a list
        # of time series appropriate for each molecule type

        nmoltypes = len(self.data[0][KEY_NMOL])
        data = [[] for mol in range(nmoltypes)]

        for step in self.data:
            for mol in range(nmoltypes):
                data[mol].append(step[KEY_NMOL][mol])

        return nmoltypes, data


    def allclose(self, other, rtol=1.0e-05, atol=1.0e-08, equal_nan=False):

        """Return True if all data are equal to within a tolerance

        Difference in numerical content is offloaded to numpy.allcose()
        so the default arguemnts are the same as this.

        Arguments:
        other (PTFILE):         to be compared
        rtol (float):           relative tolerance
        atol (float):           absolute tolerance
        equal_nan (Bool):       Are Nans equal?

        Returns:
        True or False

        Exceptions:
        If the yaml represetations are not comparable, a number of
        exceptions could arise.
        """

        same = True

        try:

            for data1, data2 in zip(self.data, other.data):

                if data1[KEY_TIMESTEP] != data2[KEY_TIMESTEP]:
                    raise ValueError("Timesteps do not match")

                for key in data1:
                    same = same and numpy.allclose(data1[key], data2[key],
                                                   rtol, atol, equal_nan)

        except (IndexError, KeyError):
            # Looks like data are not comparable
            raise

        return same


class YAMLDATA(PTFILE):

    """YAMLDATA file container"""

    def __init__(self, data, metadata):

        """Adds metadata cf PTFILE

        Arguments:
        data (yaml):        the yaml data: list of dict
        metadata (yaml):    metadata
        """

        super(YAMLDATA, self).__init__(data)
        self.metadata = metadata


def load_yaml(directory=os.curdir, yamlfile="YAMLDATA.000"):

    """Load YAMLDATA version of PTFILE

    Arguments:
    directory (string):      location
    yamlfile (string):       file name

    Returns:
    data (YAMLDATA):         data container
    """

    # safe_load_all() returns a generator, so care with
    # python2/3 compatiblity

    filename = os.path.join(directory, yamlfile)

    with open(filename, "r") as stream:
        document = yaml.safe_load_all(stream)
        metadata = next(document)
        data = next(document)

    return YAMLDATA(data, metadata)


def load(directory=os.curdir, ptfile="PTFILE.000"):

    """Load the contents of the PTFILE

    Arguemnts:
    directory (string):     location of file
    ptfile (string):        file name

    Returns:
    data (PTFILE):          data object
    """

    filename = os.path.join(directory, ptfile)

    # Read ahead to find out how many molceule types are present
    # and hence the number of lines we expect to read

    with open(filename, "r") as filecontext:

        for nblock in range(7):
            block = filecontext.readline()

        # Here is the number of lines per timestep
        nblocksize = 7 + 4*len(block.split())


    # Now read the full file

    data = []

    with open(filename, "r") as filecontext:

        nblock = 0
        block = []

        for line in filecontext:

            nblock += 1
            block.append(line)

            if nblock % nblocksize == 0:

                # Process a complete block of nblocksize lines
                # Note natom and nmol appear as float but parsed as int()

                values = block[1] + block[2] + block[3] + block[4]
                energies = [float(val) for val in values.split()]

                values = block[5]
                natom = [int(float(val)) for val in values.split()]

                values = block[6]
                nmol = [int(float(val)) for val in values.split()]

                # Per molecule data (not used)...
                values = block[7] + block[8] + block[9] + block[10]
                _ = [float(val) for val in values.split()]

                values = dict(zip(list(KEYS.keys()), energies))
                values.update({KEY_TIMESTEP: int(block[0])})
                values.update({KEY_NATOM: natom})
                values.update({KEY_NMOL: nmol})
                data.append(values)

                # Start next block
                block = []

    return PTFILE(data)
