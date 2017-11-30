"""Extends ObservableData for DL-MONTE data

The class DLMonteData extends the ObservableData class to allow
convenient access to DL-MONTE output.

The object is created by supplying the location of the directory
holding the standard DL-MONTE input and output files.
"""

# To ignore pylint numpy "no-member" errors...
# pylint: disable=E1101

import os
import subprocess

import numpy

import inputs.obs
import inputs.util
from inputs.util import Label
from inputs.util import Observable

import inputs.ensemble
import inputs.sources.dlfield as dlfield
import inputs.sources.dlconfig as dlconfig
import inputs.sources.dlptfile as dlptfile
import inputs.sources.dlcontrol as dlcontrol


class DLMonteInput(object):

    """Utility class to manage aggregate input for DL MONTE simulation"""

    def __init__(self, field, config, control, nconfigs=1):

        """Initialise with standard inputs

        Argumemts:
        field (FIELD):      the FIELD file container
        config (CONFIG):    the CONFIG file container
        control (CONTROL):  the CONTROL file container
        nconfigs (integer): number of configs expected
        """

        self.field = field
        self.config = config
        self.control = control

        self.nconfigs = nconfigs


    def __repr__(self):

        """Return details"""

        inputs = "field= {!r}, config= {!r}, control= {!r}"\
            .format(self.field, self.config, self.control)


        return "DLMonteInput({!s})".format(inputs)


    @staticmethod
    def from_directory(directory=os.curdir, nconfigs=1):

        """Read inputs from given the directory

        Arguments:
        directory (string):      existing location
        nconfigs (integer):      number of nconfig files

        Returns:
        input (DLMoneteInput):   object
        """

        filename = os.path.join(directory, "FIELD")
        field = dlfield.from_file(filename)

        filename = os.path.join(directory, "CONFIG")
        config = dlconfig.from_file(filename)

        filename = os.path.join(directory, "CONTROL")
        control = dlcontrol.from_file(filename)

        return DLMonteInput(field, config, control, nconfigs)


    def to_directory(self, directory=os.curdir):

        """Write DL MONTE FIELD, CONFIG. CONTROL files

        Arguments:
        directory (string):       directory location
        """

        filename = os.path.join(directory, "FIELD")
        with open(filename, "w") as ctxt:
            ctxt.write(self.field.__str__())

        filename = os.path.join(directory, "CONFIG")
        with open(filename, "w") as ctxt:
            ctxt.write(self.config.__str__())

        filename = os.path.join(directory, "CONTROL")
        with open(filename, "w") as ctxt:
            ctxt.write(self.control.__str__())


class DLMonteOutput(object):

    """Utility class for aggregate DL MONTE output"""

    def __init__(self, ptfile, yamldata, nconfigs=1):

        """Initialise data

        Arguments:
        ptfile (PTFILE):       PTFILE data in yaml representation
        yamldata (PTFILE):     YAMLDATA data
        nconfigs (integer):    number of configs

        Other data outputs (ARCHIVE, etc) are not represented at
        the moment.
        """

        self.ptfile = ptfile
        self.yamldata = yamldata
        self.nconfigs = nconfigs


    @staticmethod
    def load(directory=os.curdir, nconfigs=1):

        """Attempt to load all relevant output data from DL MONTE simulation

        Arguments:
        directory (string):      location
        nconfigs (integer):      number of configs expected

        Returns:
        data (DLMonteOutput):    new container object with available data
        """

        ptfile = None
        if os.path.exists(os.path.join(directory, "PTFILE.000")):
            ptfile = dlptfile.load(directory)

        yamldata = None
        if os.path.exists(os.path.join(directory, "YAMLDATA.000")):
            yamldata = dlptfile.load_yaml(directory)

        return DLMonteOutput(ptfile, yamldata, nconfigs)


class DLMonteRunner(object):

    """A utility class to help run a DL MONTE simulation"""

    def __init__(self, executable, directory=os.curdir):

        """Initialise

        Arguments:
        executable (string):      full or relative path of executable
        directory (string):       directory where input files reside
        """

        self.executable = executable
        self.directory = directory
        self.input = None
        self.output = None
        self.stderr = None


    def clone_input(self, source_dir=os.curdir):

        """Copy existing input files from source_dir to run location"""

        self.input = DLMonteInput.from_directory(source_dir)
        self.input.to_directory(self.directory)


    def execute(self, stderrfile="STDERR.000"):

        """Execute DL Monte via a subprocess"""

        try:
            self.stderr = os.path.join(self.directory, stderrfile)
            handle = open(self.stderr, "w")
        except IOError:
            # Will just have to do without ...
            handle = None

        test = subprocess.Popen(self.executable,
                                stderr=handle,
                                cwd=self.directory)
        test.wait()

        if handle is not None:
            handle.close()

        self.output = DLMonteOutput.load(self.directory)


    def remove_output(self):

        """Remove output files from the working directory"""

        files = ["ARCHIVE.000", "OUTPUT.000", "REVCON.000", "REVIVE.000",
                 "RSTART.000", "PTFILE.000", "YAMLDATA.000"]

        for name in files:
            filename = os.path.join(self.directory, name)
            if os.path.exists(filename):
                os.remove(filename)

        if self.stderr and os.path.exists(self.stderr):
            os.remove(self.stderr)


    def remove_input(self):

        """Remove input files from the working directory"""

        files = ["CONFIG", "CONTROL", "FIELD"]

        for name in files:
            filename = os.path.join(self.directory, name)
            if os.path.exists(filename):
                os.remove(filename)

    def cleanup(self):

        """Remove both inputs and outputs"""

        self.remove_input()
        self.remove_output()




class DLMonteData(htk.obs.ObservableData):

    """Container for DL-NONTE simulation data used for HTK"""

    def __init__(self, directory=os.curdir, results_dir=None):

        """Create an observable data set from DL-MONTE output

        Arguments:
        directory (string):  path to the DL-MONTE data

        The user has to specify the ensemble explicitly.
        """

        dlinput = DLMonteInput.from_directory(directory)
        estr = dlinput.control.ensemble()
        ensemble = htk.ensemble.from_string(estr)

        super(DLMonteData, self).__init__(ensemble)

        if results_dir is None:
            results_dir = directory

        self.field = dlinput.field
        self.config = dlinput.config
        self.control = dlinput.control

        self._load_data(results_dir)
        self._add_ensemble_parameters()


    def _load_data(self, directory=os.curdir):

        """Load the contents simulation

        Arguments:
        directory (string):     location

        If there is a YAMLDATA, use that. Otherwise fall back to
        the PTFILE.
        """

        output = DLMonteOutput.load(directory)
        self.data_source = directory

        if output.yamldata is not None:
            self.data_type = "DL-MONTE YAMLDATA Format"
            self._add_times(output.yamldata)
            self._add_volume_yaml(output.yamldata)
            self._add_energy_yaml(output.yamldata)
            self._add_mols_yaml(output.yamldata)
        elif output.ptfile is not None:
            self.data_type = "DL-MONTE PTFILE Format"
            self._add_times(output.ptfile)
            self._add_energies(output.ptfile)
            self._add_volume_pt(output.ptfile)
            self._add_lattice(output.ptfile)
            self._add_atoms(output.ptfile)
        else:
            # No data!
            raise IOError("No data found")


    def _add_times(self, ptfile):

        """Make ptfile time series data an observable"""

        data = numpy.array(ptfile.time_steps())
        label = Label("t", "DLM steps", "time steps")
        tobs = Observable(data, label)
        self.add_observable(tobs, independent_variable=True)


    def _add_energies(self, ptfile):

        """Look at energy data and move to obsevable"""

        energies = dict(dlptfile.ENERGY)
        eunits = self.field.units

        for key in energies:
            data = numpy.array(ptfile.time_series(key))

            if numpy.any(data):
                label = Label(key, dlptfile.KEYS[key], eunits)
                eobs = Observable(data, label)
                self.add_observable(eobs)


    def _add_volume_pt(self, ptfile):

        """In PTFILE, check to see if volume is changing"""

        data = numpy.array(ptfile.time_series("volume"))
        if numpy.any(data):
            # If volume is fixed, move to a parameter
            _, nunique = numpy.unique(data, return_counts=True)
            if nunique.size == 1:
                self._add_volume_parameter()
            else:
                self._add_volume_observable(data)


    def _add_volume_yaml(self, yamldata):

        """In YAML, volume always treated as significant if present"""

        if "volume" in yamldata.keys:
            data = numpy.array(yamldata.time_series("volume"))
            self._add_volume_observable(data)
        else:
            self._add_volume_parameter()

    def _add_volume_observable(self, vdata):

        """Add volume as observable time series"""

        label = Label("volume", "Volume", "Angstrom^3")
        self.add_observable(Observable(vdata, label))


    def _add_volume_parameter(self):

        """Add volume as a parameter"""

        label = Label("volume", "Volume", "Angstrom^3")
        volume = self.config.volume()
        self.add_parameter(volume, label)

    def _add_energy_yaml(self, yamldata):

        """In YAML, energy always treated as significant if present"""

        if "energy" in yamldata.keys:
            data = numpy.array(yamldata.time_series("energy"))
            self._add_energy_observable(data)

    def _add_energy_observable(self, data):

        """Add energy as observable time series"""

        label = Label("energy", "Energy", self.field.units)
        self.add_observable(Observable(data, label))


    def _add_energy_parameter(self):

        """Add energy as a parameter"""

        label = Label("energy", "Energy", "kT")
        self.add_parameter(energy, label)


    def _add_lattice(self, ptfile):

        """If lattice vectors fixed, ignore. Or treat together"""

        data1 = numpy.array(ptfile.time_series("L1"))
        _, nlv1 = numpy.unique(data1, return_counts=True)

        data2 = numpy.array(ptfile.time_series("L2"))
        _, nlv2 = numpy.unique(data2, return_counts=True)

        data3 = numpy.array(ptfile.time_series("L3"))
        _, nlv3 = numpy.unique(data3, return_counts=True)

        if nlv1.size > 1 or nlv2.size > 1 or nlv3.size > 1:
            label = Label("la", "Lat. vector 1", "Angstrom")
            self.add_observable(Observable(data1, label))
            label = Label("lb", "Lat. vector 2", "Angstrom")
            self.add_observable(Observable(data2, label))
            label = Label("lc", "Lat. vector 3", "Angstrom")
            self.add_observable(Observable(data3, label))

        data1 = numpy.array(ptfile.time_series("Lcos1"))
        _, nlcos1 = numpy.unique(data1, return_counts=True)

        data2 = numpy.array(ptfile.time_series("Lcos2"))
        _, nlcos2 = numpy.unique(data2, return_counts=True)

        data3 = numpy.array(ptfile.time_series("Lcos3"))
        _, nlcos3 = numpy.unique(data3, return_counts=True)

        if nlcos1.size > 1 or nlcos2.size > 1 or nlcos3.size > 1:
            label = Label("lcos1", "Angle 1", "TBC")
            self.add_observable(Observable(data1, label))
            label = Label("lcos2", "Angle 2", "TBC")
            self.add_observable(Observable(data2, label))
            label = Label("lcos3", "Angle 3", "TBC")
            self.add_observable(Observable(data3, label))


    def _add_atoms(self, ptfile):

        """Make observable from number of atoms"""

        # Units are really "None"

        natomtypes, atomdata = ptfile.natom()

        for natom in range(natomtypes):
            data = numpy.array(atomdata[natom])
            _, nau1 = numpy.unique(data, return_counts=True)
            if nau1.size > 1:
                latom = self.field.atomtypes[natom].name
                label = Label("natom" + latom, "No. atoms " + latom, None)
                self.add_observable(Observable(data, label))


    def _add_mols_yaml(self, ptfile):

        """Make observable"""

        nmoltypes, moldata = ptfile.nmol()

        for mol in range(nmoltypes):

            data = numpy.array(moldata[mol])
            lmol = self.field.moltypes[mol].name
            label = Label("nmol" + lmol, "No. molecules " + lmol, None)
            self.add_observable(Observable(data, label))


    def _add_ensemble_parameters(self):

        """Check we have parameters for the Ensemble

        Temperature: always from control file
        Volume: if constant, should be set from data
        Pressure: if required, is from control file
        N: number of particles
        """

        try:
            systemp = self.control.main_block.statements["temperature"]
            syspres = self.control.main_block.statements["pressure"]
        except KeyError:
            pass

        nlabel = Label("N", "No. atoms", None)
        tlabel = Label("systemp", "Temperature", "K")

        if isinstance(self.ensemble, htk.ensemble.EnsembleNVT):
            self.add_parameter(self.config.natom, nlabel)
            self.add_parameter(systemp, tlabel)
        elif isinstance(self.ensemble, htk.ensemble.EnsembleNPT):
            # Add pressure, temperature
            plabel = Label("syspres", "Pressure", "kAtmos.")
            self.add_parameter(self.config.natom, nlabel)
            self.add_parameter(syspres, plabel)
            self.add_parameter(systemp, tlabel)
        elif isinstance(self.ensemble, htk.ensemble.EnsembleMuVT):
            self.add_parameter(systemp, tlabel)
            # Activities
            # NEED TO CHECK MOVE TYPE
            for mover in self.control.main_block.moves[0].movers:
                activity = mover["molpot"]
                label = Label("z" + mover["id"], "Activity", "Volume^-1")
                self.add_parameter(activity, label)
        else:
            # No microcanonical examples available
            assert 0

