"""Read data associated with test Ising model in ising.py"""

# To ignore pylint numpy "no-member" errors...
# pylint: disable=E1101

import re
import numpy

import inputs.obs
import inputs.util
from inputs.util import Label
from inputs.util import Observable
from inputs.ensemble import EnsembleNVT
from inputs.parameter import Parameter

#from htk.histogram import BetaReweighter
#from htk.histogram import KTReweighter
#from htk.histogram import Reweighter

class IsingModelData(inputs.obs.ObservableData):

    """DataSet object implementation for NVT Ising model data

    It is expected that the data are read from a single file
    which may be provided to the contructor.
    """

    def __init__(self, filename=None):
        """
        Create an Ising model data set (in the NVT ensemble)
        """

        super(IsingModelData, self).__init__(EnsembleNVT())

        if filename is not None:
            self.load(filename)


    def load(self, filename=None):

        """Load data from file"""

        # Parameters

        f = open(filename, "r")

        with f:
            line = f.readline()
            line = f.readline()
            match = re.search(r" (\d+)$", line)
            n = int(match.group(0))
            v = 1.0*n*n
            line = f.readline()
            match = re.search(r" (\w+.\w+)$", line)
            j = float(match.group(0))
            line = f.readline()
            match = re.search(r" (\w+.\w+)$", line)
            h = float(match.group(0))
            line = f.readline()
            match = re.search(r" (\w+.\w+)$", line)
            kT = float(match.group(0))

        # Load the parameters

        self.add_parameter(n*n, Label("N", "Number of spins", None))
        self.add_parameter(kT, Label("kT", "Temperature", "k_bT"))
        self.add_parameter(v, Label("V", "Volume", "sites"))
        self.add_parameter(j, Label("J", "Coupling constant", "k_bT"))
        self.add_parameter(h, Label("H", "External field", "k_bT"))

        self.data_source = filename
        self.data_type = "Ising Model (2d) " + str(n) + "x" + str(n)

        # Load the observable data
        data = numpy.loadtxt(filename, skiprows=9)

        tdata = data[:, 0]
        sdata = data[:, 1]
        mdata = data[:, 2]

        # Form the total energy (per site)
        edata = sdata.copy()
        edata[:] = - j*sdata[:] - h*mdata[:]

        tobs = Observable(tdata, Label("t", "Time", "MC Sweeps"))
        sobs = Observable(sdata, Label("S", "Interaction Energy", "k_bT/site"))
        mobs = Observable(mdata, Label("M", "Magnetisation", "k_bT/site"))
        eobs = Observable(edata, Label("E", "Total Energy", "k_bT/site"))

        self.add_observable(tobs, independent_variable=True)
        self.add_observable(sobs)
        self.add_observable(mobs)
        self.add_observable(eobs)

        # Reweighters
        # Reweighting always takes place via the total energy
        # (system, not per site), so introduce a factor of the
        # volume

        vparam = self.parameter("v")
        tparam = self.parameter("kt")
        hparam = self.parameter("h")

        beta = Parameter(1.0/kT, Label("beta", "Inverse Energy", "1/k_bT"))
        rbeta = BetaReweighter("beta", beta, vparam, eobs)
        rkt = KTReweighter("kt", tparam, vparam, eobs)

        # To reweight wrt external field, a factor of v/kT is
        # required as we have magnetistation per site

        alpha = Parameter(v/kT, Label("a", "Volume/k_bT", "sites/k_bT"))
        rh = Reweighter("h", hparam, alpha, mobs)

        self.add_reweighter(rbeta)
        self.add_reweighter(rkt)
        self.add_reweighter(rh)


    def reweight_cv(self, ktnew):

        """A convenience to reweight C_V to a series of new temperatures

        Arguments:
        ktnew (float or numpy.ndarray):  the new temperatures

        Returns:
        Specfic heat capacity
        """

        volume = self.parameter("V")
        e0 = self.observable("e").data
        e1 = volume*e0[:]
        e2 = e1[:]*e1[:]

        try:
            cvnew = []
            for kt in ktnew:
                e1r = self.reweighter("kt").reweight_obs(e1, kt)
                e2r = self.reweighter("kt").reweight_obs(e2, kt)
                cv = htk.util.nvt_cv(e1r, e2r, kt, volume)
                cvnew.append(cv)

            cvnew = numpy.array(cvnew)

        except TypeError:
            e1r = self.reweighter("kt").reweight_obs(e1, ktnew)
            e2r = self.reweighter("kt").reweight_obs(e2, ktnew)
            cvnew = htk.util.nvt_cv(e1r, e2r, ktnew, volume)

        return cvnew
