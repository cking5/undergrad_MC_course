"""Observable abstract class

Dataset provided to the histogram toolkit should extend this class.
"""

from collections import OrderedDict

import numpy

import inputs.util as util
from inputs.parameter import Parameter

class ObservableData(object):

    """Top-level container for data sets
    """

    def __init__(self, ensemble):

        """Initialise an empty obsevable data set"""

        self.ensemble = ensemble
        self.params = OrderedDict()
        self.observables = []
        self.reweighters = []
        self.independent_variable = None

        self.data_source = ""
        self.data_type = ""
        self.fig_xsize = 6.0
        self.fig_ysize = 3.0


    def parameter(self, key):

        """Return parameter for given key"""

        lid = key.lower()

        return self.params[lid]


    def observable(self, key):

        """Return observable from key name"""

        keyl = key.lower()

        for obs in self.observables:
            if keyl == obs.id():
                return obs

        t = self.independent_variable
        if keyl == t.id():
            return t

        raise ValueError("No observable with label id {}".format(key))


    def reweighter(self, key):

        """Return reweight with key name"""

        keyl = key.lower()

        for reweighter in self.reweighters:
            if keyl == reweighter.rid():
                return reweighter

        raise ValueError("Non reweighter with id {}".format(key))


    def to_table(self):

        """Print information about the data set in a tabular form
        """

        fmt0 = ""
        fmt1 = "{!s:<12}"
        fmt2 = "{!s:<12} {!s:<20}"
        fmt3 = "{!s:<12} {!s:<20} {!s:<12}"
        fmt4 = "{!s:<12} {!s:<20} {!s:<12} {!s:>12}"

        tstr = []

        tstr.append(fmt2.format("Source:", self.data_source))
        tstr.append(fmt2.format("Type:", self.data_type))
        tstr.append(fmt2.format("Ensemble: ", self.ensemble))

        # parameters:
        tstr.append(fmt0)
        tstr.append(fmt1.format("Parameters:"))
        tstr.append(fmt4.format("Label", "Description", "Units", "Value"))

        for key in self.params:
            p = self.params[key]
            pstr = fmt4.format(p.label.id, p.label.name, p.label.units, p)
            tstr.append(pstr)

        tstr.append(fmt0)
        tstr.append(fmt1.format("Observables:"))
        tstr.append(fmt4.format("Label", "Description", "Units", "No. Obs."))

        t = self.independent_variable
        assert t is not None
        tstr.append(fmt4.format(t.label.id, t.label.name, t.label.units, \
                                   t.data.size))

        for obs in self.observables:
            tstr.append(obs.to_table(fmt=fmt4))

        return "\n".join(tstr)


    def add_parameter(self, value, label):

        """Add a parameter with unique label id

        Arguments:
        value (int or float):     the value of the parameter
        label (Label):            associated Label
        """

        lid = label.id.lower()

        if lid in self.params:
            raise ValueError("Label {} already exists".format(lid))

        self.params.update({lid: Parameter(value, label)})


    def add_observable(self, obs, independent_variable=False):

        """Add observable data.

        Arguemnts:
        obs (Observable):                data object (with unique label.id)
        independent_variable (Boolean):  True if this is 'time'
        """

        if independent_variable:
            if not self.independent_variable:
                self.independent_variable = obs
            else:
                raise ValueError("Independent variable already set")

        else:
            for exist in self.observables:
                key = exist.label.id.lower()
                if key == obs.label.id.lower():
                    raise ValueError("Obserable already exists: {}".format(key))
            self.observables.append(obs)


    def add_reweighter(self, re):

        """Register a reweighter"""

        for exist in self.reweighters:
            if re.rid() == exist.rid():
                raise ValueError("Reweighter exists {}".format(re.rid()))

        self.reweighters.append(re)


    def summary_time_series(self, pyplot, key, t0=0, t1=-1):

        """
        Produce a time series for an observable

        Arguments:
        pyplot (matplotlib.pyplot):       Active pyplot module
        id (string):                      Parameter key
        t0 (integer):                     First time step
        t1 (integer):                     Last time step
        """

        # We assume equally spaced observations
        # but we must allow that the time step is not unity,
        # and that the first step is not t=0.

        time = self.independent_variable
        obs = self.observable(key)

        data_t0 = time.data[0]
        dt = time.data[1] - data_t0
        n0 = 0
        if t0 > data_t0:
            n0 = int(round((t0 - data_t0)/dt))

        n1 = -1
        if t1 != -1:
            n1 = int(round((t1 - data_t0)/dt))

        fig = pyplot.figure(figsize=(self.fig_xsize, self.fig_ysize))
        self._pyplot_figure_info(fig)

        tlabel = time.label.name + " (" + time.label.units + ")"

        pyplot.xlabel(tlabel)
        pyplot.xlim(xmin=t0, xmax=time.data[n1])
        pyplot.ylabel(obs.label.name)
        pyplot.plot(time.data[n0:n1], obs.data[n0:n1], '.')

        pyplot.show()


    def summary_autocorrelation(self, pyplot, key, tlag=100):

        """
        Produce a plot of lag correlation

        Arguments:
        pyplot (matplotlib.pyplot):        Active pyplot module
        id (string):                       Observable key
        tlag (integer):                    Maximum lag (time steps)
        """

        # Construct the right time interval for the x-axis

        time = self.independent_variable
        obs = self.observable(key)

        data_t0 = time.data[0]
        dt = time.data[1] - data_t0
        nlag = int(tlag/dt)
        t = numpy.linspace(start=0.0, stop=tlag, num=nlag)

        fig = pyplot.figure(figsize=(self.fig_xsize, self.fig_ysize))
        self._pyplot_figure_info(fig, observable=obs.label.name)

        # Compute the lag correlation
        x = obs.data
        y = util.autocorrelation(x, nlag)

        pyplot.xlabel("Lag time t'")
        label = "<" + obs.id() + "(t) " + obs.id() + "(t+t')>"
        pyplot.ylabel(label)
        pyplot.plot(t, y)

        pyplot.show()


    def summary_histogram(self, pyplot, key, nbins=100):

        """
        Plot a histrogram of observable

        Arguments:

        pyplot (matplotlib.pyplot):    Active pyplot model
        id (string):                   Observable key
        nbins (integer):               Number of histogram bins
        """

        obs = self.observable(key)

        fig = pyplot.figure(figsize=(self.fig_xsize, self.fig_ysize))
        self._pyplot_figure_info(fig)

        pyplot.xlabel(obs.label.name)
        label = "P(" + obs.id() + ")"
        pyplot.ylabel(label)
        pyplot.hist(obs.data, bins=nbins, normed=True)

        pyplot.show()


    def _pyplot_figure_info(self, figure, observable=None):

        x0 = 0.95
        y0 = 0.90
        dy = 0.06
        ny = 1

        fmts = "{:<16} {:20s}"
        fmti = "{:<16} {:<20d}"

        figure.text(x0, y0 - dy*ny, fmts.format("Source:   ", self.data_source),
                    family="monospace")
        ny += 1
        figure.text(x0, y0 - dy*ny, fmts.format("Type:     ", self.data_type),
                    family="monospace")
        ny += 1
        figure.text(x0, y0 - dy*ny,
                    fmts.format("Ensemble: ", str(self.ensemble)),
                    family="monospace")
        ny += 1
        if observable is not None:
            figure.text(x0, y0 - dy*ny, fmts.format("Observable:", observable),
                        family="monospace")
            ny += 1

        figure.text(x0, y0 - dy*ny,
                    fmti.format("Measurements: ", self.nmeasure()),
                    family="monospace")
        ny += 1
        figure.text(x0, y0 - dy*ny,
                    fmti.format("Step (dt):  ", self.dt()),
                    family="monospace")
        ny += 1


    def dt(self):
        """
        Return the time interval between measurements based on the
        independnent variable
        """

        time = self.independent_variable

        assert time is not None
        assert time.data.size > 1

        dt = time.data[1] - time.data[0]
        return int(dt)


    def nmeasure(self):

        """
        Return the number of meansurements based on the independent
        variable
        """

        return self.independent_variable.data.size


    def autocorrelation_time(self, key, nmaxlag=1024):

        """
        Compute the autocorrelation time for an obsevaable
        """

        obs = self.observable(key)

        # phi(t) can become negative, or bumpy, as t becomes large.
        # To get a reasonable figure for the autocorrelation time,
        # we progressively increase the upper limit of the discrete
        # sum.

        dt = self.dt()
        phi = util.autocorrelation(obs.data, nmaxlag)

        ta0 = 0.0
        nt = 1
        while nt < nmaxlag:
            ta = numpy.sum(phi[1:nt])
            if ta < ta0: break
            ta0 = ta
            nt *= 2

        return (dt, dt*ta0)
