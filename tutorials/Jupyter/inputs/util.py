"""Utilities used by the Histogram Toolkit

Label       Container for a meaningful label and physical units
"""

import numpy

class Label(object):

    """A continer for a description for a physical quantity.

    There is no notion of what the quantity is, conversion between units,
    or anything more sophisticated at the moment.
    """

    def __init__(self, id, name, units):

        """Initialise a label for some quantity.

        E.g.,
        a = Label("K", "Absolute Temperature", "Kelvin")
        b = Label("U", "Interaction energy", "eV")

        id    (string): a short form
        name  (string): a longer form
        units (string): may be None, interpreted as "Dimensionless"
        """

        self.id = id
        self.name = name
        self.units = units

    def __str__(self):

        me = "id={!r}, name={!r}, units={!r}"  \
            .format(self.id, self.name, self.units)

        return "Label({!s})".format(me)

    def __repr__(self):

        return self.__str__()


class Observable(object):

    """A container for an observable data set."""

    def __init__(self, data, label):

        """Initialise observable data.

        data (numpy.ndarray):  the observations
        label (Label):         a description
        """

        if not label: raise ValueError("Must have a label")

        self.data = data
        self.label = label

    def __str__(self):
        me = "label={!r}, data={!r}".format(self.label, self.data)
        return "Observable({})".format(me)

    def id(self):
        return self.label.id.lower()


    def to_table(self, fmt = "{} {} {} {}"):

        str = fmt.format(self.label.id, self.label.name, self.label.units, \
                             self.data.size)

        return str


def autocorrelation(a, nmaxt):
    """
    Compute the normalised autocorrelation function.

    The lag correlator phi(t) is:
    phi(t) = ( < a(0) a(t) > - <a>^2 ) / ( <a^2> - <a>^2 )
    The normalisation is computed so that phi(0) = unity,
    and no end-effects are included.

    Arguments:
        a (numpy.ndarray): one-dimensional array of values
        nmaxt (int):       maximum lag time (nmaxt < a.size)
    Returns:
        (numpy.ndarray):   result of length nmaxt
    """

    nmax = a.size - nmaxt
    assert nmax > 0

    a1 = (1.0/nmax)*numpy.sum(a[0:nmax])
    a2 = (1.0/nmax)*numpy.sum(a[0:nmax]**2)
    phi = numpy.zeros(nmaxt, dtype = numpy.float)

    for n in range(nmaxt):
        assert (nmax+n) < a.size
        ax = (1.0/nmax)*numpy.sum(a[0:nmax]*a[n:nmax + n])
        phi[n] = (ax - a1*a1)/(a2 - a1*a1)

    return phi

def expectation_value(f):

    """
    Return the expection value of a series of measurements.
    """

    assert f.size > 0, "Zero length observable!"

    return (1.0/f.size)*numpy.sum(f[:])



def nvt_cv(e1, e2, kt, volume = 1.0):

    """Compute (specific) heat capacity in NVT ensemble.

    C_V = 1/kT^2 . ( <E^2> - <E>^2 )
    """

    cv = (1.0/(volume*kt**2))*(e2 - e1*e1)

    return cv



def multiple_free_energies(e, kT, w, nmaxit = 20):
    """
    Compute free energy offsets for multiple reweighting

    Multiple reweighting in the NVT ensemble requires the computation
    of one free energy per set of measurements. This must be done via
    an iterative method.

    Arguements:

    Returns:
        f  (numpy.ndarray):  free energy for each data set

    """

    # Look at the arguments

    nrun = len(kT)
    assert len(e) == nrun, "Check e and kT"
    assert len(w) == nrun, "Check weights and kT"

    ndata = numpy.zeros(nrun, dtype = 'int')
    fold = numpy.ones(nrun)
    f = numpy.zeros(nrun)
    beta = numpy.zeros(nrun)
    
    for irun1 in range(nrun):
        assert kT[irun1] > 0.0, "Check kT > 0"
        ndata[irun1] = e[irun1].size
        beta[irun1] = 1.0/kT[irun1]


    for n in range(nmaxit):
        # Begin iteration
        for i in range(nrun):
            f[i] = 0.0
        
            # f = free_energy(...)
            for irun1 in range(nrun):
                for idata in range(ndata[irun1]):
                    # Accumulate the denominator
                    sum = 0.0
                    for irun2 in range(nrun):
                        arg = -beta[irun2]*e[irun1][idata] + fold[irun2]
                        sum += 1.0*ndata[irun2]*w[irun2]*numpy.exp(arg)
                
                    # Accumulate term
                    arg = -beta[i]*e[irun1][idata]
                    f[i] += (1.0/sum)*w[irun1]*numpy.exp(arg)

        fold[:] = -numpy.log(f[:])

    return fold


def multiple_reweight_observable_nvt(e, obs, kT, fe, w, ktnew):

    """
    Observable reweighting based on more than one Monte Carlo data set

    Arguments:
    e     (list of numpy.ndarray):  list of energy observations
    obs   (list of numpy.ndarray):  list of observations
    kT    (list of double):         list of temperatures
    fe    numpy.ndarray:            free energy related to data set i
    w     (list of double):         list of weights
    ktnew (double):                 target temperature for reweighting in NVT
    """

    # Look at the arguments

    nrun = len(kT)
    assert len(e) == nrun, "Check e and kT"
    assert len(obs) == nrun, "Check obs and kT"
    assert len(w) == nrun, "Check weights and kT"

    ndata = numpy.zeros(nrun, dtype = 'int')
    f = numpy.zeros(nrun)
    beta = numpy.zeros(nrun)
    
    for irun1 in range(nrun):
        assert e[irun1].size == obs[irun1].size, "Check e and obs data"
        assert kT[irun1] > 0.0, "Check kT > 0"
        ndata[irun1] = e[irun1].size
        beta[irun1] = 1.0/kT[irun1]

    # Free energy at the new beta

    betanew = 1.0/ktnew
    fnew = 0.0

    for irun1 in range(nrun):
        for idata in range(ndata[irun1]):
            sum = 0.0
            for irun2 in range(nrun):
                arg = -beta[irun2]*e[irun1][idata] + fe[irun2]
                sum += 1.0*ndata[irun2]*w[irun2]*numpy.exp(arg)
                                          
            arg = -betanew*e[irun1][idata]
            fnew += (1.0/sum)*w[irun1]*numpy.exp(arg)

    fnew = -numpy.log(fnew)

    # The free energies are determined up to a constant, so subtract
    # out the target free energy
    f[:] = fe[:] - fnew
    fnew = 0.0

    # Reweight observable

    robs = 0.0
    for irun1 in range(nrun):
        for idata in range(ndata[irun1]):
            sum = 0.0
            for irun2 in range(nrun):
                arg = -beta[irun2]*e[irun1][idata] + f[irun2]
                sum += 1.0*ndata[irun2]*w[irun2]*numpy.exp(arg)
                
            arg = -betanew*e[irun1][idata] - fnew
            robs += (1.0/sum)*obs[irun1][idata]*w[irun1]*numpy.exp(arg)

    return robs
