.. _tutorial_1:

TUTORIAL 1 : NVT Lennard-Jones fluid
====================================

Authors James Grant, John Purton, r.j.grant@bath.ac.uk

Introduction to DLMONTE
-----------------------

We will start our DLMONTE tutorial with the simplest Monte Carlo simulation, a canonical, NVT simulation of a  Lennard-Jones fluid. In this example only the atoms are allowed to move and a random walk is constructed following the approach of Metropolis *et al.* [#F1]_ . There are three steps to this approach.

1. Select an atom at random (**r_1**) and calculate its energy :math:`(U(\mathbf{r}_1))`.
2. Displace this atom by a random amount to a new position (**r_2**) and calculate its new energy :math:`(U\mathbf{r}_2)`.
3. Accept/reject the displacement according to the Metropolis algorithm

.. math::

         P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2) = \min(1, \exp \{- \beta [U(\mathbf{r}_2) - U(\mathbf{r}_1)] \} )

DLMONTE Input Files
-------------------

The purpose of this test example or tutorial is to introduce some of the basic commands and file formats that are used by DLMONTE. Like DLPOLY, DLMONTE requires (at least) 3 standard input files:

FIELD, CONFIG and CONTROL 

We will start by introducing these files and will discuss each of these in turn.  
Full details of options available are given the DLMONTE manual.

FIELD
-----

The FIELD file contains the interatomic potentials and a detailed discussion of this topic can be found in the DLMONTE-2 manual. The file is similar to that used by DLPOLY, but contains some important differences.
The FIELD file for the first tutorial is shown below

.. code-block:: html
   :linenos:

   CUTOFF 2.5
   UNITS internal
   NCONFIGS 1
   ATOMS 1
   LJ core 1.0  0.0
   MOLTYPES 1
   lj
   MAXATOM 512
   FINISH
   VDW 1
   LJ core  LJ core lj   1.0 1.0
   CLOSE
   
The important points to note are the *cutoff = 2.5* for the short-range interactions in Angstroms and the energy units expected by DLMONTE are in *internal* units, the CONFIG file contains 1 configuration *NCONFIGS = 1*.  
There is one *ATOM*/particle type called called *LJ* having *mass = 1.0* and *charge = 0.0*. 
These are contained within a single molecule called *lj*. The size of this molecule is restricted to 512 particles, which in this case is the same number of particles that will be read in from the CONFIG file.
The particles interact through a *lj == Lennard-Jones* interaction having depth *\epsilon = 1* and range *\sigma = 1*.  
Unlike DL_POLY there is no one-to-one correspondence between the FIELD and CONFIG file thus we only need to specify the particle type under the keyword ATOMS.

CONFIG
------

As with DLPOLY the CONFIG file contains the positions of the particles and a sample of the configuration is provided below.

.. code-block:: html
   :linenos:
   
   Example 1: LJ NVT
   0         0
   11.7452  0.00000  0.00000
   0.00000  11.7452  0.00000
   0.00000  0.00000  11.7452
   NUMMOL 1 1
   MOLECULE lj 512 512
   LJ core
   0 0 0
   LJ core
   0 0 0.125
   ......
   
The first line is simply a title. 
In line 2 the first integer is similar to levcfg in DLPOLY and controls the amount of input and one can use this to skip the forces and velocities in output from DLPOLY configurations, note however that these will not work in DLMONTE without some modification.  
The second integer reflects the style of coordinates *0 == fractional* (but also for Cartesian: 1 == cubic lattice, 2 == orthorhombic, 3 == palleliped). 
These are followed by the lattice vectors. 

NUMMOL specifies the number of types of molecules followed by the number of each type.  
In this case there is one MOLECULE type lj, and there is one of them in the CONFIG file.  
The molecule has 512 atoms/particles to initially and is limited to a maximum number of 512.
Finally the atoms/particles are read in. 
The atom data follows and for each particle, its name (usually chemical symbol) and type (core, semi and metal that can be abbreviate to c, s, and m) are specified.  
On the next line these are followed by the coordinates in a format consistent with that given on line 2.  
 
CONTROL
-------
 
The CONTROL provides directives to DLMONTE how to undertake the calculations and switches on or off functionality. The CONTROL file in this example is:

.. code-block:: html
   :linenos:
  
   NVT simulation of Lennard-Jones fluid
   finish
   seeds 12 34 56 78               # Seed RNG seeds explicitly to the default
   nbrlist auto                    # Use a neighbour list to speed up energy calculations
   maxnonbondnbrs 512              # Maximum number of neighbours in neighbour list
   temperature     1.4283461511745 # Corresponds to T*=1.1876; T(in K) = T* / BOLTZMAN (see constants_module.f90)
   steps          10000            # Number of moves to perform in simulation
   equilibration    0              # Equilibration period: statistics are gathered after this period
   print           1000            # Print statistics every 'print' moves
   stack           1000            # Size of blocks for block averaging to obtain statistics
   sample coord   10000            # How often to print configurations to ARCHIVE.000
   revconformat dlmonte            # REVCON file is in DLMONTE CONFIG format
   archiveformat dlpoly4           # ARCHIVE.000/HISTORY.000/TRAJECTORY.000 format 
                                   # In this case: HISTORY.000 in DLPOLY4 style
   move atom 1 100                 # Move atoms with a weight of 100
   LJ core
   start
 
The first line is the title and the second contains the keyword *finish*. 
We will see later that there a number of directives in DLMONTE that *must* be placed before this directive which must be present in the CONTROL file.
*Seeds* followed by a series of 4 integers provides a reproducible seed, otheriwse *ranseed* generates a random seed from the system clock at initialisation.
The diretives *nbrlist* and *maxnonbondnbrs* control the size and administration of the neighbourlist used by DLMONTE to optimise performance, and are explained further in one of the exercises.
*Temperature* is self explanatory while *steps* is the length of the simulation in attempted  moves.
*Equilibration* etc control the detail of how data is output from DLMONTE.
The *revconformat dlmonte* instruction describes the format of the output file REVCON.000 which contains the final configuration of the simulation and can be used for continuing a simulation, *dlpoly2*, *dlpoly4* are other options.
*archiveformat dlpoly4* describes the format of the trajectory file here it will be HISTORY.000, equivalent of the HISTORY in DLPOLY4.  The *dlpoly{2,4}* formats are readable by common visualisation packages such as vmd.  Full options are detailed in the manual.

The directive *move atoms* is where we begin to touch on the fundamental control of the simulation.
The key feature here is that DLMONTE will not do anything unless told to do so (N.B. While this gives DLMONTE great flexibility it means also means that it may be possible to ask DLMONTE to perform ill-defined calculations). 
In this simple *NVT* ensemble only the particles move, thus *move atoms* tells DLMONTE to move *1* atom type, with a weight of  *100*. 
The line(s) following this detail the atom and type.
Finally the *start* directive ends the *CONTROL* file and instructs DLMONTE to start the simulation.

.. You are now ready to run DLMONTE

Running and Output
------------------

A successful DLMONTE calculation will produce a number of output files::

* OUTPUT.000 contains details of the simulation, statistics, running time, or errors if the calculation failed.
* REVCON.000 contains the final configuration in the format specified
* PTFILE.000 contains statistics though will eventually be deprecated in favour of
* YAMLDAT.000 which contains statistics in the yaml format
* ARCHIVE.000/HISTORY.000/TRAJECTORY.000 contains the trajectory in the specified format
* note that GCMC is not yet supported by standard packages.

For analysis we will typically process the YAMLDAT.000 and visualise the trajectory files.  
However for understanding how the simulation proceeds it is useful to have some familiarity with the OUTPUT file.
The file begins with a header detailing the version, authors and suggested citations, followed by the brief summary of details of the simulation as specified in the input files.
As section headed *simulation parameters* then specifies all parameter values that will be used within the simulation.

The final step before starting the calculation is to determine the initial energy of the system and the details of this are printed in a block::

 --------------------------------------------------
                  initial energies
 --------------------------------------------------

 break down of energies for box:   1

 total energy                       -0.7037212307E+03
 reciprocal space coulomb            0.0000000000E+00
 real space coulomb                  0.0000000000E+00
 external mfa coulomb                0.0000000000E+00
 nonbonded two body (vdw)           -0.7037212307E+03
 bonded two body (pair)              0.0000000000E+00
 nonbonded three body                0.0000000000E+00
 bonded three body (angle)           0.0000000000E+00
 bonded four body (angle)            0.0000000000E+00
 many body energy                    0.0000000000E+00
 external potential energy           0.0000000000E+00
 total virial                        0.0000000000E+00
 volume                              0.1620247087E+04

This is followed by a partial breakdown per molecule type in the system and the time taken to initialise the calculation.  There after every *print* steps as specified in the CONTROL file a block is printed to the output::

 Iteration       1000 - elapsed time (seconds)     0.0730

      step      en-total            h-total             coul-rcp            coul
 -real
      step      en-vdw              en-three            en-pair             en-a
 ngle
      step      en-four             en-many                                 volu
 me
      step      cell-a              cell-b              cell-c
      step      alpha               beta                gamma
      r-av      en-total            h-total             coul-rcp            coul
 -real
      r-av      en-vdw              en-three            en-pair             en-a
 ngle
      r-av      en-many             vir-tot             volume              pres
 sure
      r-av      cell-a              cell-b              cell-c
      r-av      alpha               beta                gamma
 ----------------------------------------------------------------------------------------------------
      1000    -0.7650295178E+03   -0.6130828709E+02    0.0000000000E+00    0.0000000000E+00
              -0.7650295178E+03    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
               0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
               0.1174520000E+02    0.1174520000E+02    0.1174520000E+02
               0.9000000000E+02    0.9000000000E+02    0.9000000000E+02


              -0.7351584460E+03   -0.3143721526E+02    0.0000000000E+00    0.0000000000E+00
              -0.7351584460E+03    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
               0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
               0.1174520000E+02    0.1174520000E+02    0.1174520000E+02
               0.9000000000E+02    0.9000000000E+02    0.9000000000E+02

 LJ       c        512.0000        512.0000

 lj                1.0000          1.0000
 ----------------------------------------------------------------------------------------------------

This specifies the same breakdown in a tabulated form, including the iteration number and time taken.  By using the command::

  grep time OUTPUT.000

you can quickly gauge the progress of your calculation.
DLMONTE runs a *check* of the system at regular intervals to verify that the system is behaving correctly.
This typically involves calculating the energy from scratch and comaring it with the running total, the result is printed to the OUTPUT.000 with a line looking like::

 Workgroup    0, box    1 check: U_recalc - U_accum =  0.32742E-10  0.27570E-10  0.53847E-13 -0.32246E-13 (internal, kT, kT/atom, dU/U)

The final value is the relative energy difference which should be of the order of the working precision, typically ending E-13 or E-14.

Finally at the end of the simulation a summary block detailing the average values during the simulation (excluding the first *equilbration* steps) and their fluctations, the final energies, and 'processing data' detailing time, move data and final parameters::

 ----------------------------------------------------------------------------
                          processing data


 ----------------------------------------------------------------------------

 total no of atom moves         :     10000
 successful no of atom moves    :      5279      0.52790000

 displacement (Angstroms) for LJ       c :   0.5346

 total elapsed time (seconds)          0.7730
 normal exit

The OUTPUT.000 of a succssfully completed job will end with *normal exit*.

.. Here begins Balena specific submission

Submitting your job
-------------------

We will be using the University of Bath's HPC, Balena for the workshop.
You should have received a crib sheet on accessing Balena.
To begin with submit your job using the command::

[username@balena-01 tutorial1]$ sbatch single.sub

Most jobs in the workshop can be run using this script.  
You can monitor the job using::

[user@balena-01 tutorial1]$ squeue -u $USER

In order to analyse output and visualise trajectories you will need to log into a copmute node using the command::

[user@balena-01 tutorial1]$ sint

To visualise your calculation run vmd::

[user0@node-sw-039 tutorial1]$ vmd

To visualise the trajectory::

   File
      New molecule

   Determine file type
      -> Select DLPOLY V3 History
      -> Browse and Select HISTORY.000 

Now add the line::

  yamldata 1000

to the control file and rerun the simulation.

To analyse the trajectory, e.g. the energy evolution during the calculation, first run the script:: 

[user0@node-sw-039 tutorial1]$ strip_yaml.sh energy

and then plot the time sequence::

   [user0@node-sw-039 tutorial1]$ gnuplot
   gnuplot> plot './energy.dat' u 1:2

.. Here ends Balena specific submission

.. Links to extension exercises

Now try the extension exercises to learn more about funcitonality within DLMONTE and to optimise your calculation:

  :ref:`tut1_ex1`

  :ref:`tut1_ex2`

  :ref:`tut1_ex3`

With each exercise we recommend you create a copy of the inputs in a sub-directory::

[user0@node-sw-039 tutorial1]$ mkdir ex1
[user0@node-sw-039 tutorial1]$ cp CONFIG CONTROL FIELD ex1
[user0@node-sw-039 tutorial1]$ cd ex1

.. Link to next tutorial

Or move on to  :ref:`tutorial_2` and NPT simulations.

.. rubric:: Footnotes

.. [#f1] N. Metropolis, A.W. Rosenbluth, M.N. Rosenbluth, A.N. Teller, and E. Teller. Equation of state calculations by fast computing machines. J. Chem. Phys. , 21:1087 1092, 1953.
