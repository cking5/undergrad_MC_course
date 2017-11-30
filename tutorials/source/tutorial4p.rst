--------------------------------
Tutorial 4: Semi-grand Ensembles
--------------------------------

Authors: Chris King, James Grant - r.j.grant@bath.ac.uk

Introduction:
=============

So far in this course, we have looked at three types of ensemble that can be used in Monte Carlo simulations: NVT (Canonical), NpT (Gibbs), and :math:`\mu`VT (Grand Canonical).  In the first two ensembles, the total number of objects in the simulated system, *N*, remains constant and is used in other deterministic simulation techniques, such as Molecular Dynamics.  These are the simplest ensemble conditions to model computationally, and so are widely used, despite their limitations, which you encountered in previous sessions.  In these cases, if one is using Monte Carlo techniques, one can apply the :math:`\mu`VT ensemble, where now the chemical potential, :math:`\mu`, is kept constant, instead of *N*, allowing objects to be inserted/removed from the system.  In this session, we will be looking at another, slightly different ensemble, known as the semi-grand canonical ensemble (or simply 'semi-grand' ensemble, SG).  SG ensembles are used exclusively in Monte Carlo simulations describing systems with two or more objects, such as mixtures.  In SG ensembles, the volume or pressure, the temperature and the total number of objects,*N*, are kept constant, but the composition of the system is allowed to change, *i.e.* if:

.. math::

   N = \sum_i N_i
   
where :math:`N_i` is the total number of objects of type *i* in the system.  So, in SG ensembles, *N* must remain constant, but :math:`N_i` can vary.  This allows so-called 'exchange' moves, where an object of one type is exchanged with that of another type.  E.g. If we have a system with two types of particles, A and B, with :math:`N_A` particles of type A and :math:`N_B` particles of type B, the total number of particles in the system is :math:`N = N_A + N_B`.  If we propose an exchange move where one particle of type A is changed to type B, if the move is accepted, then the total number of particles after the move is :math:`N = (N_A - 1) + (N_B + 1) = N_{Afin} + N_{Bfin}`.

The total system energy in the SG ensemble of a system of pure A particles, :math:`E_{N_{A/B}}` is given by:

.. math::
 
   E_{N_{A/B}} = \DeltaV_{AA/AB} - \Delta\mu\N_{A/B}
   
where :math:`\Delta\mu` is the difference in chemical potential between the two particle types, which in our case is :math:`\Delta\mu = \mu_A - \mu_B` and :math:`\DeltaV_{AA/AB} = V_{AA} - V_{AB}` is the interaction potential between either two A particles or between A and B particles, respectively.  This can then be inserted into the Metropolis Acceptance criterion for proposed move from an initial configuration with :math:`N_{Ai}` particles of type A and :math:`M_{Bi}` particles of type B to a new configuration containing :math:`N_{Af}` A particles and :math:`M_{Bf}` B particles:

.. math::

         P_{\mathrm{acc}}((N_{Ai}, N_{Bi}) \rightarrow (N_{Af}, N_{Bf})) = \min(1, \exp \ \Bigl(- \frac{E(N_{Ai}, N_{Bi}) - E(N_{Af}, N_{Bf})}{kT}\Bigr) \ )

where *k* is the Boltzmann constant and *T* is temperature.

The 2D Square Ising Model used in session 1 made use of this ensemble, unbeknownst to you.  In that example, our system was a square lattice of two types of particles representing the two possible orientations of the spin degree of freedom.  The moves proposed there were exchanges between spin 'up' and spin 'down' particles, representing a spin flip in the lattice.

SG ensembles are employed in any simulations of mixtures, these include: alloys and solid solutions, liquid mixtures and solutions, colloids, and isomeric and polymorphic transitions.  SG is particularly useful in exploring how the system behaves as its composition changes.

In this session, you will explore the use of the SG ensemble by investigating how the composition of a Lennard-Jones system comprised of two types of Lennard-Jones particles, A and B, with different interaction potentials and different sizes affects the behaviour of the system.  Principally, we will look at how the particles behave when changing both the chemical potential of the particles and the temperature.  The DL_MONTE input files for the system we will be simulating in this session are shown below:

CONFIG
------

The first part of the DL_MONTE CONFIG file is shown below, you may notice that it is identical to the CONFIG from Session 2.  So our initial configuration is a face-centred cubic (fcc) structure of pure A.

.. code-block::
   :linenos:
   
   NCONFIG 1 
            0         0
          12.6400000000        0.0000000000        0.0000000000
           0.0000000000       12.6400000000        0.0000000000
           0.0000000000        0.0000000000       12.6400000000
   NUMMOL 1 1 
   molecule lj     2048     2048
   A    CORE
           0.5312500000        0.5312500000        0.5312500000
   A    CORE
           0.5312500000        0.5937500000        0.5937500000
   ...

You will notice that no particles of B are present in the initial configuration, this is becuse we will initially be attempting to introduce particles of B into the system.

CONTROL
-------

The CONTROL file that will be used in this session is shown below:

.. code-block::
   :linenos:
   
   Semi-grand simulation of Lennard-Jones system 
   finish
   ranseed               # Seed RNG seeds explicitly to the default                          
   nbrlist auto                    # Use a neighbour list to speed up energy calculations
   maxnonbondnbrs 80                # Maximum number of neighbours in neighbour list  
   verlet 0.0
   temperature     2269 # T* = \frac{2}{\ln(1+sqrt(2))} 
   steps           100000            # Number of moves to perform in simulation 
   equilibration    0              # Equilibration period: statistics are gathered after this period
   print           1000         # Print statistics every 'print' moves
   stats           1000
   stack           10              # Size of blocks for block averaging to obtain statistics
   sample coords   100000         # Sample coordinates to ARCHIV every <N> steps
   revconformat dlmonte            # REVCON file is in DL_MONTE CONFIG format
   archiveformat dlmonte           # ARCHIV file is in DL_MONTE CONFIG format
   #yamldata 1000                   # Print yamldata every <N> steps
   move semigrandatoms 1 1 100.0     # Move type semigrandatoms <n> <f> <deltamu>
   A core B core                   # Semi-grand exchanges of atom A and B 
   #move atom 1 1
   #A core
   #move atom 1 1
   #B core
   move volume cubic linear 1
   start 
   
As you can see, the basic format of the CONTROL file is unchanged, except for one new line, line 17, which instructs DL_MONTE to attempt moves which swap A to B and vice-versa.  The first two integers represent how many atoms are removed from the system and how frequently swaps are attempted, respectively.  The final number represents :math:`\Delta\mu`.

|think| Referring to the expression for system energy in the semi-grand ensemble, does :math:`\Delta\mu > 0` favour the insertion of B particles?  What about :math:`\Delta\mu < 0`?

.. |think| image:: images/General/think.png
   :height: 100 px
   :scale: 25 %

FIELD
-----

The simulations we will be running will use potentials described by the following FIELD file:

.. code-block::
   :linenos:
   
   semi-grand LJ field file
   CUTOFF 2.5
   UNITS internal
   NCONFIGS 1
   ATOMS 2
   A core 1.0  0.0
   B core 1.0  0.0
   MOLTYPES 1
   lj
   MAXATOM 2048
   FINISH
   VDW 3
   A core  A core lj   1.0 1.0
   B core  B core lj   1.0 1.1
   A core  B core lj   1.0 1.05
   CLOSE
   
The layout is the same as the previously used FIELD files, where we define two neutral particles: A and B, which both have a mass of 1.0.  We have one molecule type, 'lj', that has a maximum number of 2048 atoms.  The three potentials define the interaction between two A particles, two B particles and between A and B particles, all of them are described by the Lennard-Jones potential:

.. math::

   \phi(r_{ij}) = 4\epsilon\biggl[\Bigl(\frac{\sigma}{r_{ij}}\Bigr)^{12}-\Bigl(\frac{\sigma}{r_{ij}}\Bigr)^{6}\biggr]
   
where :math:`\epsilon` = 1.0 eV for each interaction but :math:`\sigma` is larger for B particles.  |think| Based on this, how does the size of B particles compare to that of A particles?

Exercise 1)
===========

The first thing you will be doing is looking at how many B particles are inserted into the system with changing :math:`\Delta\mu`.  Let's begin by making :math:`\Delta\mu` increasingly positive.  

|action| Run DL_MONTE on the initial input files.  The calculation should take about a minute to complete.

.. |action| image:: images/General/action.png
   :scale: 5 %
   
|action| Once your calculation is complete, identify and extract the total number of B present in the system at the end of the calculation from your output data. (HINT: look to the OUTPUT.000 file). |action| Record this value and the value of :math:`\Delta\mu`.

|action| Repeat the calculation with increasingly positive values of :math:`\Delta\mu` by changing the corresponding value in the CONTROL file.  Make sure you copy the input files into a new directory for each separate calculation that you undertake. |action| Record the number of B present in the system at each value of :math:`\Delta\mu`.

|think| How does the value of :math:`\Delta\mu` affect the ease at which A is swapped for B? Is this what you expect?

|action| Now make the value of :math:`\Delta\mu` negative and run the calculation.  |think| How does this affect the number of successful swap moves?

|action| Repeat the calculations with increasingly negative values of :math:\Delta\mu` and |think| identify the value of :math:\Delta\mu` at which the system remains pure A, *i.e.* no B is present in the system at the end of the calculation.

Exercise 2)
===========

We have seen how :math:`\Delta\mu` affects the amount of B that is added to the system, but what is the final distribution of B in the system? Are B particles spread homogeneously through the system or do they aggregate into larger 3D structures? To answer this question, one can either look at the final system configuration files directly (*i.e.* view the REVCON.000 file) or one can calculate the *radial distribution function* (rdf) of B particles in the system.  

The rdf describes the density variation of one or more particles as a function of separation measuring from one reference particle.  In essence, it tells you how many particles can be found at a given distance from a chosen reference particle (see Figure 1).  If one plots the rdf against distance, as in Figure 2, you will observe a curve with several peaks.  Each corresponds to an increase in the density of particles, correlating to an increase in the number of particles at that distance from the reference.  To summarise, the rdf can be used to find the number of particles type is present and how far apart they are from each other by looking at how the density changes with distance from an arbitrary reference.  In DL_MONTE, the option to collect data for an rdf is specified by the 'sample rdfs *n* :math:`R_c` *x*' directive in the CONTROL file.  The entire system is divided into *n* equally-sized bins and the density variation around each particle in the system is calculated up to a maximum separation of :math:`R_c` and this data is collected every *x* number of Monte Carlo steps.  The rdf data is stored in the RDFDAT.000 file.

.. figure:: images/Tut_4p_images/RDF-illustrations-med.png
   :align: center

   **Figure 1**: An example of how the rdf is constructed from a system in a given configuration, using an arbitrary reference particle.
   
.. figure:: images/Tut_4p_images/RDF-illustrations-med.png
   :align: center
   
   **Figure 2**: Example rdf plot for a Lennard-Jones system.  Each peak corresponds in an increase in the density of particles, and thus the number of particles at a given distance, can be found and the configuration of the system can thence be inferred.
   
|action| For each of your previous calculations, plot the RDFDAT.000 data, showing how density varies with particle separation.

|think| From your plots, do B particles tend to cluster together or remain spread out relatively homogeously thoughout the system? Explain your answer.

|think| Do these observations reflect the behaviour of real solid-state binary systems?

Exercise 3)
===========

We have seen how the value of :math:`\Delta\mu` affects the ease at which B is added to our initially pure A system, but you may recall from the acceptance criterion in the Metropolis Algorithm, the temperature of the system, *T*, will also have an effect on the probability of accepting a move, and hence the ease of which B is added to the system.

|think| By examining the acceptance criterion, how do you think the temperature will affect the amount of B added to the system under a constant positive :math:`\Delta\mu`? Explain your answer.

|action| Select one of your positive :math:`\Delta\mu` calculations.  |action| Create a new directory and copy the input files for this calculation into it.  |action| Now change the temperature value in the CONTROL file and re-run the calculation.  |action| Once the calculation completes, record the total number of B particles present in the system.

|action| Repeat the calculation at different temperature values and record the final number of B particles present.  |think| From your observations, how does *T* affect the ease at which B is added to the system? Is this what you expect? Explain your answer.

Extension:
----------

|action| Plot the rdfs for each of your temperature calculations.  

|think| Compare you new rdf plot.  Does the temperature have an affect on the shape of the rdf? |think| Is this what you expect?

Conclusions:
============

In this session, you have been introduced to the semi-grand ensemble and demonstrated its use in simulations of mixtures.  You have applied the semi-grand ensemble to a solid mixture of two types of Lennard-Jones particles, A and B, and investigated how both temperature and the difference in chemical potential between the two particles affects the ideal composition of the system.  You have been introduced to the radial distribution function and used it to examine the distribution of added B particles in the system.  In the next session, you will be given an open-ended practical exercise to investigate using the knowledge and experience you have gathered throughout this course.

Extensions:
===========

1. This session has had you look at a system which is initially pure solid A, but one can also start from the pure B solid and investigate the addition of A to the system in a manner analagous to that conducted in this session.  

|action| To do this, go to the 'solid B' directory, where you will find input files for the solid B system.  You will notice that the FIELD and CONTROL files are identical to those used in the pure A system.  However, the CONFIG file is slightly different, you will notice that our system is larger and the B particles have different positional coordinates than the CONFIG file for the pure A solid.  This is due to the fact that B is a different size to A, and so the lattice parameters of the pure solid will be different.

|action| Investigate the composition of the system and the distribution of added A particles in a manner analagous with that done in this session.  |action| Determine how :math:`\Delta\mu` and *T* affect the final composition of the system and |action| determine whether clusters of A are preferable over a roughly homogenous distribution.

|think| Compare your results of adding A to solid B and vice versa, are there any similarities between your results? Rationalise any observed differences in terms of the relative sizes of the two particles.