------------------------------------------------------
*TUTORIAL 3: Introduction to Ensembles and Move Types*
------------------------------------------------------

Authors: Chris King, James Grant - r.j.grant@bath.ac.uk

Introduction
============

In the previous session, you were introduced to the ideas underlying the Monte Carlo simulation technique by modelling a simple system and extracting useful output data which (hopefully) concurred with theory.  In this session, you will be formally introduced to the general Monte Carlo program, DLMONTE, which you unknowingly used in the previous session.  This session will begin to touch upon the concepts of ensembles and how the choice of ensemble affects the outcome of a simulation and explore the possible types of moves available in each ensemble.

Ensembles are very important in thermodynamics and statistical mechanics, where one can readily observe the macroscopic properties of a system without being able to know or control its microscopic properties.  In simple terms, an ensemble can be thought of as a collection of a large number of replicas of a system, akin to running repeating an experiment many times with the same initial conditions.  By averaging the outcomes of each of these replicas, one can create a probability distribution of possible outcomes for an experiment in a given ensemble.  In this course and statistical mechanics you will come across several different ensembles which are described in the form of three letters, describing which parameters of the system are kept constant for the ensemble.  Ensembles are defined in this manner to ensure that statistical equilibrium is maintained, *i.e.* ensembles not to evolve over time, these include:

- *N*: total number of individual objects (atoms, molecules, particles etc.) in the system kept constant.
- *V*: system volume kept constant
- *P*: system pressure kept constant
- *T*: system temperature kept constant
- *E*: the total energy of the system kept constant
- :math:`\mu`: the chemical potential of the system kept constant (you'll see this later in the course)

For example, 'NVT' means that the total number of particles, total volume and temperature are kept constant (known as the Canonical Ensemble) and 'NPT' means that the total number of particles, total pressure and temperature are kept constant (a.k.a. the Gibbs ensemble).  Figure 1 displays physical representations of the different ensembles that are used in thermodynamics.

.. figure:: images/Tut_3_images/Statistical_Ensembles.png
   :align: center
   
   **Figure 1**: Visual representation of different ensembles [#f1]_. 

As you may recall, a key part of Monte Carlo simulations is the sampling over all possible states of the system by 'moving' through configuration space, where a move is a change from an initial configuration to a new one.  So far in this course we have looked at moves which consist of physically moving an object from one set of coordinates to another (within reason), which are known as translational moves.  But this is not the only move that can be made in Monte Carlo simulations, Other moves are possible, depending on the type of ensemble used.  For instance, one can define a volume move, where the proposed move is changing the total system volume by altering one or more length scales of the system and as we will see later in the course one can also perform insert/delete and swap moves of objects. 

The sequence of moves defines the 'trajectory' of the system in configurational space.  It is important to understand that this is not the same as the physical motion of individual objects in the system, or the system itself, which is the traditional definition of a trajectory that you are probably used to. There is no 'time' in Monte Carlo simulations; nothing *in* the system evolves with time.  This limits Monte Carlo to simulation of static systems only, so it cannot determine any dynamic properties of the system, like diffusion coefficients or rates of reaction..  Any reference to time in this context refers to the computational time required to complete the calculation. 

Part 1: NVT Lennard-Jones fluid
===============================

To start with, you will be given an overview of the program that you will use (or have used) throughout this course: DLMONTE.  DLMONTE is the name of a program which provides a general Monte Carlo framework for use in scientific research.  This program is a direct Monte Carlo analogue to the Molecular Dynamics program DL_POLY (which you used in Tutorial 1). 

The best way to learn how to use DLMONTE is by applying it to a simple system, in this case determining the thermal behaviour of a Lennard-Jones solid.  This part of the tutorial is directly analagous to Tutorial 1, except now we shall be using Monte Carlo to simulate our Lennard-Jones material instead of Molecular Dynamics.  This will hopefully help highlight the differences in their respective approaches to this problem.

The purpose of this part of the tutorial is to introduce some of the basic commands and file formats that are used by DLMONTE. DLMONTE requires at least 3 input files: FIELD, CONFIG and CONTROL files and we will discuss each of these in turn.

CONFIG
------

The CONFIG file contains the positions of the particles and a sample of the configuration is provided below.

.. code-block:: html
   :linenos:

   Example 1: LJ NVT               # Title                      
   0         0                     # Integers describing how the input is read in and the style of coordinates, respectively
   11.7452  0.00000  0.00000       # These lines describe the dimensions of the system in terms of basis lattice vectors, with 'x y z' components, respectively
   0.00000  11.7452  0.00000       # Since our system is 3D, we need three basis vectors to fully describe it
   0.00000  0.00000  11.7452       # In this case, the system is a cube with sides of length 11.7452 Angstroms
   NUMMOL 1 1                      # NUMMOL specifies the number of types of molecules followed by the number of each type. In this case there is one MOLECULE type lj, and there is one of them in the CONFIG file. 
   MOLECULE lj 512 512             # The molecule has 512 atoms/particles to initially and is limited to a maximum number of 512
   LJ core                         # Now each particle is read in to the file, in the form:
    0 0 0                          # NAME core
   LJ core                         # x y z position in the system, the first particle is located at the origin
    0 0 0.125                      # and the second is located at x = 0, y = 0, z = 0.125 
   ......

The rest of the file reads in the rest of the particles in our molecule in the same format.  If there were more than type of molecule, it would read in all the atoms/particles in each molecule sequentially.

CONTROL
-------
 
The CONTROL file provides parameters and conditions to DLMONTE about how to undertake the calculations. The CONTROL file in this example is:

.. code-block:: html
   :linenos:
  
    NVT simulation of Lennard-Jones fluid  # Title
    finish                                 # Needs to be here, some conditions must be placed before this word, but there aren't any in this case 
    seeds 12 34 56 78                      # Sets the initial configuration
    nbrlist auto                           # Use a neighbour list to speed up energy calculations (don't worry about this)
    maxnonbondnbrs 512                     # Maximum number of neighbours in neighbour list
    temperature     1.4283461511745        # Temperature of the system in Kelvin
    steps          10000                   # Number of moves to perform over the course of the simulation
    equilibration    0                     # Equilibration period: statistics are gathered after this period
    print           1000                   # Print statistics every 'print' moves to the output file
    stack           1000                   # Number of moves over which average values are calculated
    sample coord   10000                   # How often to print configurations to ARCHIVE.000
    revconformat dlmonte                   # REVCON file (the final configuration) is in DLMONTE CONFIG format
    archiveformat dlpoly4                  # Sets the format for the ARCHIVE.000/HISTORY.000/TRAJECTORY.000 files, in this case: HISTORY.000 in DLPOLY4 style
    move atom 1 100                        # Move atoms called 'LJ' with a weight of 100
    LJ core
    start                                  # Tells DLMONTE to begin calculation

The directives *nbrlist* and *maxnonbondnbrs* control the size and administration of the neighbourlist used by DLMONTE to optimise performance, and where a detailed explanation is given in one of the extensions.  The key feature here is that DLMONTE will not do anything unless told to do so (*N.B.* While this gives DLMONTE great flexibility it means also means that it may be possible to ask DLMONTE to perform ill-defined calculations). 

FIELD
-----

The FIELD file contains a full description of the interatomic potentials present in the system.  The FIELD file used in this instance is:

.. code-block:: html
   :linenos:
   
    Lennard-Jones                    # Title
    CUTOFF 2.5                       # The maximum distance between two particles at which the potential energy is calculated
    UNITS internal                   # Set the units of energies, internal = 10 J mol^-1
    NCONFIGS 1                       # Number of configurations described in the CONFIG file
    ATOMS 1                          # Number of atom types in the system
    LJ core 1.0  0.0                 # In this case there is one atom type called 'LJ' with mass = 1.0 and charge = 0.0
    MOLTYPES 1                       # Number of molecule types in the system...
    lj                               # ...called 'lj'...
    MAXATOM 512                      # ...with a maximum number of 512 atoms
    FINISH                           # Completes the list of atom and molecule types in the system
    VDW 1                            # The number of potentials present in the system, this must be the same number as the number of interaction types presented on the following lines
    LJ core  LJ core lj   1.0 1.0    # Defines the interaction between two LJ atoms as a Lennard-Jones (lj) potential with epsilon = 1.0 and sigma = 1.0
    CLOSE                            # This ends the FIELD file once all interaction are described
    
The CUTOFF keyword is defined as 2.5:math:`\sigma` by convention.  The UNITS can also be electron volts (eV), kJmol\ :sup:-1 \, kJ or kcal.  The *NCONFIGS* keyword refers to the number of configurations and this is usually set to 1. You will see in later tutorials that more than one types of interactions can be defined at the end of the file.

For more information on these files, refer to the DLMONTE manual in 'this directory'.

You are now ready to run DLMONTE

Exercise 1
----------

The aim of this exercise is to mirror the exercises from the first session and will hopefully illustrate another way to model a Lennard-Jones solid.  In this case, we will simulate the system under the NVT ensemble at various temperatures in order to estimate the melting point of the solid.

What types of moves are possible in this ensemble?

[instructions to get onto Balena]We will be using the University of Bath's HPC, Balena for the workshop.  You should have received a crib sheet on accessing Balena.

Once you have successfully logged onto Balena, type the following command into the command line::

   cd tutorial NVT inputs filepath

and press 'Enter'.  This will navigate you to the appropriate directory where you will run your calculations and analyse your data.  

In general you can use::

   cd thefilepathintowhateverfolderyouwant

to navigate from one folder to another on Balena.  *N.B.* to go back one directory, type '../' into the cd command, write this for each folder that you want to go back from, *i.e.* if you wanted to go back three folders, you would type::

   cd ../../../

You can then add the filepath of the folder to the end of this to navigate in one fluid line.  You can auto-complete by pressing the 'Tab' key once (this can be a good way to check whether a folder name that you type is present in that directory).  If you can't remember what files are present or lose where you are while typing your command, you can double-tap 'Tab' to list all files in your current directory.

Open the folder called '++++', you will see your DLMONTE input files: CONFIG, CONTROL and FIELD, as well as some scripts that you will use to analyse your output data.

Now you will submit a calculation to Balena from this folder using the current input files.

Submitting your job
^^^^^^^^^^^^^^^^^^^

To begin with submit your job using the command::

   [username@balena-01 tutorial1]$ sbatch single.sub

Most jobs in the workshop can be run using this script.  
You can monitor the job using::

   [user@balena-01 tutorial1]$ squeue -u $USER

Outputs
^^^^^^^

A successful DLMONTE calculation will produce a number of output files:

* OUTPUT.000 contains details of the simulation, statistics, running time, or errors if the calculation failed.
* REVCON.000 contains the final configuration in the format specified
* PTFILE.000 contains statistics though will eventually be deprecated in favour of...
* YAMLDAT.000 which contains statistics in the yaml format
* ARCHIVE.000/HISTORY.000/TRAJECTORY.000 contains the trajectory in the specified format

In this exercise we will analyse the YAMLDAT.000 and visualise the trajectory files.  
However for understanding how the simulation proceeds it is useful to have some familiarity with the OUTPUT file.

*N.B.* The OUTPUT.000 of a successfully completed job will end with 'normal exit'.

To visualise the trajectory of your system, you can use a program called 'vmd' by running the following commands::

   [user0@balena-01 tutorial1]$ module load vmd
   [user0@balena-01 tutorial1]$ vmd

Once it has loaded, you should see a 'VMD Main' window, in this window, go to:

   File :math:`\rightarrow` New molecule :math:`\rightarrow` Determine file type :math:`\rightarrow` Select DLPOLY V3 History :math:`\rightarrow` Browse and Select HISTORY.000 

Repeat the calculation at different temperatures.  Create a new folder for each new temperature and copy the CONFIG, CONTROL and FIELD files from one of your other calculations into it.  Then change the temperature value in the CONTROL file to a value of your choosing (HINT: you won't need to go above ++).  Run this calculation in the same manner as described above.  Do this for a range of temperatures and identify the melting temperature by visualising it's evolution using VMD.  

Exercise 2:
-----------

So far in this course, we have assumed that the system has reached equilibrium with its surroundings, *i.e.* that the system has reached its most thermodynamically stable state with minimal net exchange of energy with its surroundings.  Equilibration is incredibly important to Monte Carlo (and many other computational modelling techniques) as it ensures reproducibility of results.  If we start from an arbitrary initial state with a given set of parameters, the first stage of the calculation will be establishing equilibrium, with the output during this period being of little use and should be omitted from any statistical analysis of the output.  In DLMONTE (and DLPOLY) we account for this period o time using the 'equilibration' parameter in the CONTROL file.  This states the point at which output data is included in any statistical analysis.  This 'equilibration time' will be different for every system with a given set of initial parameters and is usually estimated during preliminary analysis of the data.

One way of determining when a system has reached equilibrium is by plotting the time evolution of total energy over the course of the simulation, which is what you will now do:

Navigate to one of your completed calculations and run the following command:: 

   [user0@node-sw-039 tutorial1]$ strip_yaml.sh energy

This will give you a yaml file containing the total energy after a given number of steps has completed.  Plot this data, either by using software packages like Excel, or by running the following commands in your directory::

   [user0@node-sw-039 tutorial1]$ gnuplot
   gnuplot> plot './energy.dat' u 1:2
   gnuplot> set term png
   gnuplot> set output energyvst.png

This opens a data plotting program called 'gnuplot' and tells it to plot the file called 'energy.dat' and save it as a png file in your directory.  While in gnuplot, you can navigate to different directories with the 'cd' command, but with the filepath in 'quotes'.  To exit gnuplot, type::
   
   exit 

From these energy plots, how can you tell whether the system has equilibrated? Estimate the equilibration time for your system.

How do you think the equilibration time will change with temperature? Explain your answer.

Exercise 3:
-----------

In this part of the tutorial, we will again be looking at the phase transition of a Lennard-Jones solid, but under the NPT ensemble.  This allows not only translational moves of individual particles, but also volume moves (system expansion/contraction).

Navigate to the folder 'NPT', you will find the same input files as in the NVT ensemble and a directory for each temperature that you will run a calculation.  The CONFIG and FIELD files are unchanged but the CONTROL has a few modifications:

.. code-block:: html
   :linenos:
  
   NPT simulation of Lennard-Jones material # Title
   finish                                   # Needs to be here, some conditions must be placed before this word, but there aren't any in this case 
   seeds 12 34 56 78                        # Sets the initial configuration
   nbrlist auto                             # Use a neighbour list to speed up energy calculations (don't worry about this)
   maxnonbondnbrs 512                       # Maximum number of neighbours in neighbour list
   temperature     1.4283461511745          # Temperature of the system in Kelvin
   pressure     0.0179123655568             # pressure of the system
   steps          10000                     # Number of moves to perform over the course of the simulation
   equilibration    0                       # Equilibration period: statistics are gathered after this period
   print           1000                     # Print statistics every 'print' moves to the output file
   stack           1000                     # Number of moves over which average values are calculated
   sample coord   10000                     # How often to print configurations to ARCHIVE.000
   revconformat dlmonte                     # REVCON file (the final configuration) is in DLMONTE CONFIG format
   archiveformat dlpoly4                    # Sets the format for the ARCHIVE.000/HISTORY.000/TRAJECTORY.000 files, in this case: HISTORY.000 in DLPOLY4 style
   move atom 1 100                          # Move atoms called 'LJ' with a weight of 100
   LJ core
   move volume cubic linear 1               # Move volume, the system is cubic, linearly scaled with a weight of 1
   start

Specifically, with the additional lines::

   move volume cubic linear 1     

which is the instruction to introduce volume moves as *move volume*, *linear* refers to how volume is sampled, and the inclusion of pressure::

   pressure     0.0179123655568

In these calculations, volume moves are attempted less frequently than translational moves, this is because typically volume moves are more computationally intensive than single atom moves, why do you think that this is the case?

Run and analyse the output data in the same manner as for the previous exercise. Ensure that the system has equilibrated at each instance.  Remember to create a new directory for each temperature you attempt. 

Examine the evolution of your system using VMD, rationalise any observed differences between the behaviours of the system under the NVT and NPT ensembles?

Estimate the melting point of the Lennard-Jones solid under the NPT ensemble.  How does it compare with the value you obtained from the NVT calculations?

Plot the total energy of the system as a function of temperature under both NPT and NVT ensembles on the same graph.  How do they compare with each other? (HINT: think about the different types of energy transfer that could be taking place in each case.)

Additionally, by using the command::

   strip_yaml.sh volume

you can extract the time evolution of the system volume from YAMLDATA.000.  Plot this data for each temperature on the same graph.  What trends do you observe as you change the temperature? Is this what you expect from a material? 

For at least one of your calculations, plot the volume and energy time evolutions on the same graph, are there any similarities between the shape of the two plots?

Conclusions:
============

After this session, you should now be familiar with the input/output files of DLMONTE as well as running calculation with the program.  You have demonstrated its use by running simulations on the simple Lennard-Jones solid system and confirmed that it shows thermodynamic behaviour consistent with real materials.  You have been introduced to the concept of ensembles in thermodynamics, in particular the NVT and NPT ensembles.  You should also have an appreciation for the possible types of Monte Carlo moves that can be proposed within NVT and NPT ensembles and the differences between them.  In the next session, we will move onto simulations under the :math:`\mu` VT (Grand Canonical) ensemble and troduce the concept of detailed balance in the Monte Carlo technique.

Extensions (optional):
======================

1. Lennard-Jones phase diagram
------------------------------ 

You have found the melting point of the system by varying the temperature while keeping the pressure constant.  Now you will examine how the melting point of the system changes with both temperature and pressure.  You can readily change the pressure of the system by altering the associated value in the CONTROL file.  Create a new directory for each pressure value you make, and for each pressure, run the calculation at various temperatures and estimate the melting point of the system.  How does pressure and temperature affect the melting point of the system?

.. figure:: images/Tut_3_images/LJ_phase_diagram.png
   :align: center

   **Figure 1**: Phase diagram for the Lennard-Jones system [#f2]_

Figure 1 shows the phase diagram for the Lennard-Jones system, compare your results with the phase diagram.  Why do you not see the coexistence of solid-liquid phases in your system?

2. Neighbour lists
------------------

DLMONTE uses neighbour list to improve the performance of the energy calculation, particles only have to check interactions with their neighbours, not every particle in the simulation.  This is particularly beneficial when particles retain the same numbers for the whole simulation or for any attempted moves. In DLMONTE, this functionality is described by the following lines in the CONTROL file::

   nbrlist auto

This rebuilds a particle's neighbourlist whenever necessary.  The size of the neighbour list is determined by::

   maxnonbondnbrs <int>

This determines the memory allocated for each particles neighbourlist.  The size will be determined by the size of your system, its density and the interaction cut-off as specified in the FIELD file.  

Go to the directory named 'helloneighbour' and view the CONTROL file, you will notice that it looks identical to the CONTROL files that you have already seen, with one extra line called::

   verlet <float>  

Run calculations using different values for this parameter and see how they affect the time taken to complete the calculation.  You can extract this for a given calculation by using the command::

   grep "total elapsed" OUTPUT.000

or alternatively the script::

   time.sh

How does tuning the parameter affect the duration of this calculation?  Why might this be the case?

*N.B.* For short simulations system time can dominate the apparent performance of the calculation.  To see this affect try running DLMONTE consecutively and checking the duration.

3. Different Sampling Schemes
-----------------------------

It is also possible to conduct a random walk in the logarithm of the volume.  Create a new directory and copy the CONFIG, CONTROL and FIELD files from one of your completed calculations into it.  Open the CONTROL file and change the keyword *linear* to  *log* in the volume move command.  This changes the way in which the the move is generated; logarithmic rather than linear, and is often claimed to be more efficient.  The acceptance criterion for the move in the Metropolis algorithm is now:

.. math::

  P_{\mathrm{acc}}([\mathbf{r}_{1},V_1] \rightarrow [\mathbf{r}_2,V_2]) = \min(1, \exp \{- \beta [U(\mathbf{r}_2) - U(\mathbf{r}_1) + P_{ext}(V_{2}-V_{1}) - ( N + 1 ) \beta^{-1} \ln(V_{2} / V_{1}) ] \} )
         
Extract the time evolution of the volume for this calculation and compare it with the volume evolution from the equivalent linear calculation.  Rationalise the observed differences.

.. rubric:: Footnotes

.. [#f1] Statistical_Ensembles.png - Wikipedia Commons: from - https://commons.wikimedia.org/wiki/File:Statistical_Ensembles.png Author: NZjacobmartin

.. [#f2] B. L. Holian, "Shear viscosities away from the melting line: A comparison of equilibrium and nonequilibrium molecular dynamics", *J. Chem. Phys.*, **78**, 11, pp. 5147-5150, 1983.