"""A two-dimensional Ising Model implementation"""

import sys
import numpy
import numpy.random

class IsingModel(object):

    r"""
    A two-dimensional Ising Model

    The Hamiltonian is
    H  = - J \sum_ij s_i s_j - mu H \sum_i s_i =  -J S - mu H M

    Here, we have coupling contant j and external magnetic field
    h  (susceptibility mu is set to unity). The thermal energy
    scale is kT.
    Spins are s_i = +/- 1.

    This is two dimensions: \sum_ij is over four nearest neighbours.
    The boundaries are periodic.
    """

    def __init__(self, nlen, j, h, kT, seed, init = 'hot'):

        """
        A square system of nlen by nlen is created with parameters
        j and h, and temperature kT. A random number seed is
        supplied to numpy.random().
        """

        if numpy.mod(nlen, 2): raise ValueError("Please use even nlen")

        self.nlen = nlen
        self.j = j
        self.h = h
        self.kT = kT
        self.s = numpy.ndarray((nlen, nlen), dtype = numpy.int)
        self.seed = seed
        numpy.random.seed(seed)
        self.init(init)
        self.verbose = False
        self._file = None

    def init(self, init):

        """
        Initialise (or re-initialise) the system.
        Arguments:
        init (string) -- init is either 'cold' for all spins -1 or
                         'hot' for random initialisation with no net spin
        """

        linit = init.lower()
        if not (linit == 'hot' or linit == 'cold'):
            raise ValueError("init should be hot or cold")

        self.initial_state = init
        self.s[:, :] = -1
        if linit == 'cold': return

        nlen = self.nlen
        naccept = 0

        while naccept < nlen*nlen/2:
            ic = numpy.random.randint(nlen)
            jc = numpy.random.randint(nlen)
            if self.s[ic, jc] == -1:
                naccept += 1
                self.s[ic, jc] = +1


    def run(self, nsteps, file = None, report_freq = 1, ndiscard = 0,
            random_update = False):

        """
        Run a number of MC steps and produce some information
        """

        n = 0
        self._report_open(file)

        while n < nsteps:

            n += 1
            self._monte_carlo_sweep(random_update)

            if n > ndiscard and numpy.mod(n, report_freq) == 0:
                self._report_update(n)

        self._report_close()


    def observables(self):
        """
        Return current observable values S, M. These are normalised
        by the system size (ie., average per site). Recall tat the
        total energy is - J S - mu H M (with mu = 1).
        """

        nlen = self.nlen
        s = 0.0
        m = 0.0

        for ic in range(nlen):
            im = numpy.mod(ic - 1 + nlen, nlen)
            ip = numpy.mod(ic + 1, nlen)
            assert im >= 0 and im < nlen
            assert ip >= 0 and ip < nlen

            for jc in range(nlen):
                jm = numpy.mod(jc - 1 + nlen, nlen)
                jp = numpy.mod(jc + 1, nlen)
                assert jm >= 0 and jm < nlen
                assert jp >= 0 and jp < nlen

                s0 = self.s[ic, jc]

                s += 0.5*s0*(self.s[im, jc] + self.s[ip, jc] + self.s[ic, jm] + self.s[ic, jp])
                m += s0

        # Compute averge per site
        s = s/(nlen*nlen)
        m = m/(nlen*nlen)

        return [s, m]

    def _monte_carlo_sweep(self, random_update = False):

        """
        Single MC sweep across lattice. Returns the number of
        accepted moves for this step.

        If you want detailed balance, use a random update; for
        general use, ordered update gives shorter correlation
        times and so is quicker.
        """

        nlen = self.nlen
        j = self.j
        h = self.h
        de = 0.0
        naccept = 0

        for i0 in range(nlen):
            ic = i0
            if random_update: ic = numpy.random.randint(nlen)
            im = numpy.mod(ic - 1 + nlen, nlen)
            ip = numpy.mod(ic + 1, nlen)

            for j0 in range(nlen):
                jc = j0
                if random_update: jc = numpy.random.randint(nlen)
                jm = numpy.mod(jc - 1 + nlen, nlen)
                jp = numpy.mod(jc + 1, nlen)

                # Existing state and new trial state (just flip
                # with probability one)

                s0 = self.s[ic, jc]
                s1 = -s0

                # Compute difference in energy

                ds = s1 - s0
                delta = -j*ds*(self.s[im, jc] + self.s[ip, jc] + self.s[ic, jm] + self.s[ic, jp]) - h*ds

                # Metropolis

                if delta < 0.0:
                    self.s[ic, jc] = s1
                    naccept += 1
                    de += delta
                elif numpy.random.uniform() < numpy.exp(-delta/self.kT):
                    self.s[ic, jc] = s1
                    naccept += 1
                    de += delta

        return naccept, de


    def _report_open(self, filename):

        self.av = {'s' : 0.0, 'm' : 0.0}
        self.sq = {'s' : 0.0, 'm' : 0.0}
        self.ncount = 0

        if filename is None: return

        f = open(filename, "w")
        f.write("Ising Model Report\n")
        f.write("System size: {}\n".format(self.nlen))
        f.write("J          : {}\n".format(self.j))
        f.write("mu H       : {}\n".format(self.h))
        f.write("kT         : {}\n".format(self.kT))
        f.write("Random seed: {}\n".format(self.seed))
        f.write("State (t=0): {}\n".format(self.initial_state))
        f.write("Observable : S interaction energy/J \n")
        f.write("Observable : M magnetization\n")
        self._file = f

    def _report_update(self, nt):

        [s, m] = self.observables()

        self.ncount += 1
        self.av['s'] += s
        self.av['m'] += m

        self.sq['s'] += s*s
        self.sq['m'] += m*m

        if self._file is None: return

        self._file.write("{:7d} {:14.7e} {:14.7e}\n".format(nt, s, m))

    def _report_close(self):

        r = 0.0
        if self.ncount > 0:
            r = 1.0/self.ncount

        for k in self.av:
            self.av[k] *= r
        for k in self.sq:
            self.sq[k] *= r

        if self._file is None: return

        self._file.write("# Summary\n")
        self._file.write("# Samples, mean observables:\n")
        self._file.write("# {:7d} {:14.7e} {:14.7e}\n".format(self.ncount, self.av['s'], self.av['m']))
        self._file.write("# Samples, mean square observables\n")
        self._file.write("# {:7d} {:14.7e} {:14.7e}\n".format(self.ncount, self.sq['s'], self.sq['m']))
        self._file.close()

        if self.verbose:
            sys.stdout.write("Wrote results to {:s}\n".format(file))
