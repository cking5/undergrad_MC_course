--------------------------------
Tutorial 4: Semi-grand Ensembles
--------------------------------

Authors: Chris King, James Grant - r.j.grant@bath.ac.uk

Introduction:
=============

So far in this course, we have looked at three types of ensemble that can be used in Monte Carlo simulations: NVT (Canonical), NpT (Gibbs), and :math:`\mu`VT (Grand Canonical).  In the first two ensembles, the total number of objects in the simulated system, *N*, remains constant ans can also be used in other deterministic simulation techniques, such as Molecular Dynamics.  These are the simplest ensemble conditions to model computationally, and so are widely used in the study of many different phenomena.  However, there are limitations to using these ensembles, which you have seen in previous sessions.  In these cases, if one is using Monte Carlo techniques, one can apply the :math:`\mu`VT ensemble, where now the chemical potential, :math:`\mu`, is kept constant, instead of *N*, allowing objects to be inserted/removed from the system as a simulation progresses, as you have seen in the previous session.  In this session, we will be looking at another, slightly different ensemble, known as the semi-grand canonical ensemble (or simply 'semi-grand' ensemble, SG).  SG ensembles are used exclusively in Monte Carlo simulations of mixtures.  In SG ensembles, the volume or pressure, the temperature and the total number of objects are kept constant, but the composition of the system is allowed to change, *i.e.* if:

.. math::

   N = \sum_i N_i
   
where :math:`N_i` is the total number of objects of type *i* in the system.  So, in SG ensembles, *N* must remain constant, but :math:`N_i` need not.  This allows so-called 'exchange' moves, where an object of one type is exchanged with that of another type.  E.g. If we have a system with two types of particles, A and B, with :math:`N_A` particles of type A and :math:`N_B` particles of type B, the total number of particles in the system is :math:`N = N_A + N_B`.  If we propose an exchange move where one particle of type A is changed to type B, if the move is accepted, then the total number of particles after the move is :math:`N = (N_A - 1) + (N_B + 1) = N_{Afin} + N_{Bfin}`.

The 2D Square Ising Model that you were introduced to in session 1 used such an ensemble, unbeknownst to you.  In that example, our system was a square lattice of two types of particles, A and B, representing the two possible orientations of the spin degree of freedom.  The moves proposed there were exchanges between A and B particles, representing a spin flip in the lattice.

SG ensembles are employed in any simulations of mixtures, these include: alloys and solid solutions, liquid mixtures and solutions, colloids, and isomeric and polymorphic transitions.  SG is particularly useful in exploring how the system behaves as its composition changes.

In this session, you will explore the use of the SG ensemble by investigating how the composition of a Lennard-Jones system comprised of two types of Lennard-Jones particles with different interaction potentials and different sizes affects the behaviour of the system.  Principally, we will look at how the particles move and segragate by changing the chemical potential of the two particles and then by changing the temperature.  The DL

