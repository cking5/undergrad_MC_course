------------------------------------
TUTORIAL 2: Introduction to DL_MONTE
------------------------------------

Authors: Chris King, James Grant - r.j.grant@bath.ac.uk

In the previous session, you were introduced to the ideas underlying Monte Carlo simulations by modelling a simple system and extracting useful output data which (hopefully) concurred with theory.  In this session, you will be formally introduced to the general Monte Carlo program, DL_MONTE, which you unknowingly used in the previous session to plot the thermodynamic behaviour of a simple Lennard-Jones material using Monte Carlo methods.  This is a direct analogy with the system studied in the first session using Molecular Dynamics and we will be able to compare the results from both methods in this session.  

Part 1: NVT Lennard-Jones material
==================================

To start with, you will be given an overview of the program that you will use (or have used) for all Monte Carlo calculations throughout this course: DL_MONTE.  DL_MONTE is the name of a program which provides a general Monte Carlo framework for use in scientific research.  This program is a direct Monte Carlo analogue to the Molecular Dynamics program DLPOLY (which you used in Tutorial 1). 

The best way to learn how to use DL_MONTE is by applying it to a simple system, in this case determining the thermal behaviour of a Lennard-Jones solid.  This part of the tutorial is directly analagous to Tutorial 1, except now we shall be using Monte Carlo to simulate our Lennard-Jones material instead of Molecular Dynamics.  This will hopefully help highlight the differences in their respective approaches to this problem.  In this part of the tutorial, each of our simulations will contain a constant number of particles, *N*, have a constant volume, *V*, and be at a fixed temperature, *T*, so these are called 'NVT' simulations. 

The purpose of this part of the tutorial is to introduce some of the basic commands and file formats that are used by DL_MONTE. DL_MONTE requires at least 3 input files: FIELD, CONFIG and CONTROL files and we will discuss each of these in turn.

CONFIG
------

The CONFIG file describes the initial configuration of our system before the simulation begins.  A CONFIG file used in this tutorial is shown below:

.. code-block:: html
   :linenos:

   Example 1: LJ NVT               # Title                      
   0         0                     # Integers describing how the input is read in and the style of coordinates, respectively
   11.7452  0.00000  0.00000       # These lines describe the dimensions of the system in terms of lattice vectors, with 'x y z' components, respectively
   0.00000  11.7452  0.00000       # Since our system is 3D, we need three basis vectors to fully describe it
   0.00000  0.00000  11.7452       # In this case, the system is a cube with sides of length 11.7452 Angstroms
   NUMMOL 1 1                      # NUMMOL specifies the number of types of molecules followed by the number of each type. 
   MOLECULE lj 512 512             # The molecule has 512 atoms/particles to initially and is limited to a maximum number of 512
   LJ core                         # Now each particle is read in to the file, in the form: NAME core
    0 0 0                          # x y z position in the system, the first particle is located at the origin
   LJ core                         # And the second is located at x = 0, y = 0, z = 0.125
    0 0 0.125                      # And so on until all atoms are defined 
   ......

The first line of the file gives a short title to describe what system the file is showing, it is not read by DL_MONTE so could be anything the author wishes.  The second line shows two numbers, the first relates to how any inputs for the file are read by DL_MONTE and the second describes the format of the coordinates of each object in the system.  Lines 3-5 describe the three dimensions of our system: x, y and z.  In this case, our system is a cube with sides of length 11.7452 Angstroms.  The *NUMMOL* line specifies how many different types of molecule we have in our system followed by the number of each type present.  In this case, we have one type of molecule and there is only one of it in our system.  The molecule's name and the number of atoms it contains are stated on the following line.  In this case, our molecule is called 'lj' and contains 512 Lennard-Jones particles and can only have up to 512 particles in it. The rest of the file describes the name and position of every particle present in the system in the form:

.. code-block:: html
   
   NAME core
   x y z


CONTROL
-------
 
The CONTROL file is where you set all the conditions that you want your calculation to obey in DL_MONTE. The CONTROL file in this tutorial is given below:

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
    revconformat DL_MONTE                   # REVCON file (the final configuration) is in DL_MONTE CONFIG format
    archiveformat dlpoly4                  # Sets the format for the ARCHIVE.000/HISTORY.000/TRAJECTORY.000 files
    move atom 1 100                        # Move atoms called 'LJ' with a weight of 100
    LJ core
    start                                  # Tells DL_MONTE to begin the calculation

The first line, as in our CONFIG file, is simply a title to describe the system that this CONTROL file belongs to.  *finish* tells DL_MONTE that there are no conditions that need to be immediately stated.  Line 3 describes a combination of four 'seeds' from which the initial configuration is based ('grown').  The lines beginning with *nbrlist* and *maxnonbondnbrs* are used by DL_MONTE to optimise its performance.  The line below these specifies what the temperature of our system is going to be, in Kelvin.  *steps* allows you to set how long the simulation will last by specifying how many steps you want it to take, where a step is one Monte Carlo move proposal and outcome.  Equilibration is how many steps at the beginning of the simulation that are needed for the system to reach equilibrium.  Lines 9-13 allow you to adjust the format of the outputs of DL_MONTE.  The *move atom* line states that we want a Monte Carlo move to represent moving a particle in our system from one position to another.  In our case, our particles are called 'LJ', corresponding to the name of the atoms in the CONFIG file.  The key feature of the CONTROL file is that DL_MONTE will not do anything unless told to do so (*N.B.* While this gives DL_MONTE great flexibility it means also means that it may be possible to ask DL_MONTE to perform ill-defined calculations). 

FIELD
-----

The FIELD file contains a full description of the interatomic potentials present in the system.  An example of a FIELD file used in this tutorial is shown below:

.. code-block:: html
   :linenos:
   
    Lennard-Jones                    # Title
    CUTOFF 2.5                       # The maximum distance between two particles for which the interaction energy is calculated
    UNITS internal                   # Set the units of energies, internal = 10 J mol^-1
    NCONFIGS 1                       # Number of configurations described in the CONFIG file
    ATOMS 1                          # Number of atom types in the system
    LJ core 1.0  0.0                 # In this case there is one atom type called 'LJ' with mass = 1.0 and charge = 0.0
    MOLTYPES 1                       # Number of molecule types in the system...
    lj                               # ...called 'lj'...
    MAXATOM 512                      # ...with a maximum number of 512 atoms
    FINISH                           # Completes the list of atom and molecule types in the system
    VDW 1                            # The number of potentials present in the system
    LJ core  LJ core lj   1.0 1.0    # Defines the interaction between two LJ atoms 
    CLOSE                            # This ends the FIELD file once all interaction are described
    
Like with the other input files, the first line is a title which describes the system being described by the file.  *CUTOFF* states the maximum distance between two particles that DL_MONTE will calculate the interaction energy between them, if the distance between two particles is greater than this, then the interaction energy is neglible and ignored by DL_MONTE.  The cutoff is defined as 2.5 :math:`\sigma` by convention.  The *UNITS* line tells DL_MONTE what units of energy you want to use, in DL_MONTE, you can choose between electron volts (eV), kJmol\ :sup:`-1` \, kJ, kcal or 'internal'.  Internal units are the simplest units for the program to use and for DL_MONTE, this is 10 J mol\:sup:`-1`.  Sometimes, a CONFIG file can specify several possible configurations, so *NCONFIGS* tells DL_MONTE how many configurations are present in the associated CONFIG file.  In our case (and all future cases), there is only one configuration in the CONFIG file. The line beginning with ATOMS simply states how many different atoms are present in the CONFIG file.  In our case, there is only one (the 'LJ' in the CONFIG).  The number stated on this line must match the number of atom types present in the CONFIG file.  The next line describes all the different atom types in the form:

.. code-block:: html

   NAME core mass charge

In our case, we have only one atom type, called 'LJ' of mass = 1.0 and no charge.  *FINISH* tells DL_MONTE that all the atom types have been defined and it can move onto the next lines.  Lines 11 and 12 describe all the interactions present in the system, *VDW* tells DL_MONTE how many different interactions are present in the system.  The line(s) after this list each interaction between all relevant atom types. In our case we only have one interaction that is between two of our 'LJ' atoms that is described by a Lennard-Jones ('lj') potential.  For the Lennard-Jones potential, one must give a value for both :math:`\epsilon` and :math:`\sigma`, in this case, :math:`\epsilon = 1.0` eV and :math:`\sigma = 1.0` Angstroms.  Other potentials may require a different number of constants to be defined, these are tabulated in the DL_MONTE manual.

For more information on these files, refer to the DL_MONTE manual in 'this directory'.

Exercise 1)
-----------

The aim of this exercise is to mirror some of the exercises from the first session and will hopefully illustrate another way to model a Lennard-Jones solid.  In this case, we will simulate this system under 'NVT' conditions, with simulations at various temperatures in order to estimate the melting point of the solid.

As you may recall, a key part of Monte Carlo simulations is sampling over all possible states of the system by 'moving' through configuration space, where a move is a change from an initial configuration to a new one.  The most intuitive move is translational moves which consist of physically moving an object from one set of coordinates to another (within reason).  But other types of moves are possible, depending on the type of ensemble used.  For instance, one can define a volume move, where the proposed move is changing the total system volume by altering the length scales of one or more dimensions of the system and (as we will see later in the course) one can also perform insert/delete and swap moves of objects.

|think| What types of moves are possible in a system under these conditions?

.. |think| image:: images/General/think.png
   :height: 100 px
   :scale: 25 %

**instructions for running a calculation here**

Outputs
^^^^^^^

A successful DL_MONTE calculation will produce a number of output files:

* OUTPUT.000 contains details of the simulation, statistics, running time, or errors if the calculation failed.
* REVCON.000 contains the final configuration in the format specified
* PTFILE.000 contains statistics though will eventually be deprecated in favour of...
* YAMLDAT.000 which contains statistics in the yaml format
* ARCHIVE.000/HISTORY.000/TRAJECTORY.000 contains the trajectory in the specified format

In this exercise we will analyse the YAMLDAT.000 and visualise the trajectory files.  
For your understanding of how the simulation proceeds it may nonetheless be useful to have some familiarity with the OUTPUT file.

*N.B.* The OUTPUT.000 of a successfully completed job will end with 'normal exit'.

The sequence of moves defines the 'trajectory' of the system in configurational space.  It is important to understand that this is not the same as the physical motion of individual objects in the system, or the system itself, which is the traditional definition of a trajectory. There is no 'time' in Monte Carlo simulations; nothing *in* the system evolves with time.  This limits Monte Carlo to simulation of static systems only, so it cannot determine any dynamic properties of the system, like diffusion coefficients or rates of reaction.  Any reference to time in this context refers to the computational time required to complete the calculation.
 
|action| Repeat the calculation at different temperatures.  Create a new folder for each new temperature and copy the CONFIG, CONTROL and FIELD files from one of your other calculations into it.  Change the temperature value in the CONTROL file to a value of your choosing (HINT: you won't need to go above ++).  Run this calculation in the same manner as described above.  Do this for a range of temperatures.

.. |action| image:: images/General/action.png
   :scale: 5 % 

|think| Identify the melting temperature of the system both by visualising the trajectories of your system at each temperature and by plotting energy vs temperature in the same way that you did in Tutorial 1.  Are the melting temperatures identified in your simulations consistent?

Remember, you can use the 'collate temperature and plot' script to plot the system energy against temperature, once you have run calculations at a variety of temperatures.

|think| How does your estimate(s) of the melting point compare with that based on your Molecular Dynamics calculation?

Exercise 2)
-----------

So far in this course, we have assumed that the system has reached equilibrium with its surroundings, *i.e.* that the system has reached its most thermodynamically stable state with minimal net exchange of energy with its surroundings.  This concept of 'equilibration' is incredibly important to Monte Carlo (and many other computational modelling techniques) as it ensures reproducibility of results.  If we start from an arbitrary initial state with a given set of parameters, the first stage of the calculation will be establishing equilibrium, with the output during this period being of little use and should be omitted from any statistical analysis of the output.  In DL_MONTE (and DLPOLY) we account for this period of time using the 'equilibration' parameter in the CONTROL file.  This states the point at which output data is included in any statistical analysis.  This 'equilibration time' will be different for every system with a given set of initial parameters and is usually estimated during preliminary analysis of the data.

One way of determining when a system has reached equilibrium is by plotting the time evolution of total energy over the course of the simulation, which is what you will now do.

|action| Navigate to one of your completed calculations and run the following command:: 

   [user0@node-sw-039 tutorial1]$ strip_yaml.sh energy

|think| From these energy plots, how can you tell whether the system has equilibrated? Estimate the equilibration time for your system.

|think| How do you think the equilibration time will change with temperature? Explain your answer.

As you may recall, detailed balance is a sufficient condition for ensuring that our simulation reflects the laws of thermodynamics.  It is generally stated as:

.. math::

   W(\mathbf{r}_1 \rightarrow \mathbf{r}_2)P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2) = W(\mathbf{r}_2 \rightarrow \mathbf{r}_1)P_{\mathrm{acc}}(\mathbf{r}_2 \rightarrow \mathbf{r}_1)

where :math:`W(\mathbf{r}_1 \rightarrow \mathbf{r}_2)` is the statistical weight of moving from an initial configuration, :math:`\mathbf{r}_1` to a final configuration, :math:`\mathbf{r}_2` (and vice-versa for :math:`W(\mathbf{r}_2 \rightarrow \mathbf{r}_1)`) and :math:`P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2)` is the probabilty of accepting the move from :math:`\mathbf{r}_1` to :math:`\mathbf{r}_2` (and similarly for :math:`P_{\mathrm{acc}}(\mathbf{r}_2 \rightarrow \mathbf{r}_1)`).

|think| Describe the condition for detailed balance in this series of simulations, where only translational moves are permitted.

Part 2: NPT Lennard-Jones Material:
===================================

Exercise 3)
-----------

In this part of the tutorial, we will again be looking at the phase transition of a Lennard-Jones solid, but this time, the each simulation will be run such that they contain a fixed number of particles, *N*, have a constant pressure, *P*, and are at constant temperature, *T*, these are known as 'NPT' simulations.  This allows not only translational moves of individual particles, but also volume moves (system expansion/contraction).

|action| Navigate to the folder 'NPT', you will find the same input files that you used in the previous exercises, if you have a look at these, you will see that the CONFIG and FIELD files are unchanged but the CONTROL has a few modifications:

.. code-block:: html
   :linenos:
  
   NPT simulation of Lennard-Jones material 
   finish                                    
   seeds 12 34 56 78                        
   nbrlist auto                             
   maxnonbondnbrs 512                       
   temperature     1.4283461511745          
   pressure     0.0179123655568             # pressure of the system
   steps          10000                     
   equilibration    0                       
   print           1000                     
   stack           1000                     
   sample coord   10000                     
   revconformat DL_MONTE                     
   archiveformat dlpoly4                    
   move atom 1 100                          
   LJ core
   move volume cubic linear 1               # Move volume, the system is cubic, linearly scaled with a weight of 1
   start

Specifically, with the additional lines::

   move volume cubic linear 1     

which is the instruction to introduce volume moves as *move volume*, *linear* refers to how volume is sampled, and the inclusion of pressure::

   pressure     0.0179123655568

In these calculations, volume moves are attempted less frequently than translational moves, this is because typically volume moves are more computationally intensive than single atom moves. |think| Why do you think that this is the case?

|action| Run the calculation at the same temperature values you used in exercise 1.  Ensure that the system has equilibrated at each instance.  Remember to create a new directory for each new temperature. 

|action| Examine the evolution of your system at each temperature and compare them with their equivalent NVT calculation. |think| Rationalise any observed differences between the behaviours of the NVT and NPT systems.

|think| Estimate the melting point of the NPT simulations of the Lennard-Jones material.  |think| How does it compare with the value you obtained from the NVT calculations?

|action| Plot the total energy of the system as a function of temperature under for both NPT and NVT calculations on the same graph.  |think| How do they compare with each other? (HINT: think about the different types of energy transfer that could be taking place in each case.)

|action| Additionally, by using the command::

   strip_yaml.sh volume

you can extract the time evolution of the system volume from YAMLDATA.000.  |think| Plot this data for each temperature on the same graph.  |think| What trends do you observe as you change the temperature? Is this what you expect from a real material? 

|action| For at least one of your calculations, plot the volume and energy time evolutions on the same graph. |think| Are there any similarities between the shape of the two plots?

Exercise 4)
-----------

You have found the melting point of the system by varying the temperature while keeping the pressure constant.  Now you will examine how the melting point of the system changes with both temperature and pressure.  You can readily change the pressure of the system by altering the associated value in the CONTROL file.  |action| Navigate to the 'changep' directory, open the CONTROL file, and change the pressure (*N.B.* values anywhere between x and y should be sufficient).  |action| Create a new directory for each pressure, and within each of these directories, run the calculation at various temperatures and |think| estimate the melting point of the system at each pressure.  

|think| How does the melting point vary with both temperature and pressure? Is this consistent with the behaviour of real materials?

*N.B.* make sure to copy the strip_yaml script into each new directory you make.

.. figure:: images/Tut_2_images/LJ_phase_diagram.png
   :align: center

   **Figure 2**: Phase diagram of the Lennard-Jones system, plotting (reduced) temperature against (reduced) density [#f1]_.

Figure 2 shows the phase diagram for the Lennard-Jones system.  *N.B.* Don't be put off by the fact that density is shown instead of pressure, they are equivalent in our system.

|action| Compare your results with Figure 2.  

|think| Why do you not see the coexistence of solid and liquid phases in your system?

*N.B.* If you want to know about reduced units, try the following sources: (1)_, (2)_

.. _(1): http://www4.ncsu.edu/~franzen/public_html/CH795N/modules/ar_mod/comp_output.html,  

.. _(2): http://cbio.bmt.tue.nl/pumma/index.php/Manual/ReducedUnits

Conclusions:
============

After this session, you should now be familiar with the input/output files of DL_MONTE as well as running calculation with the program.  You have demonstrated its use by running simulations on a simple Lennard-Jones solid system and confirmed that it shows thermodynamic behaviour consistent with real materials.  You have compared the results from Monte Carlo and Molecular Dynamics techniques and understood the differences between them.  You should also have an appreciation for the possible types of Monte Carlo moves that can be proposed, depending on the constraints of our simulation and the differences between them.  In the next session, we will move onto simulations under the :math:`\mu`\VT (Grand Canonical) ensemble, where the total number of particles in the system is not constant.

Extensions (optional):
======================

1.  Move size update
--------------------

DL_MONTE is able to automatically tune the size of attempted moves to optimise performance. By altering the maximum proposed move size during the simulation DL_MONTE is able to optimise for the particular problem.

|think| If the proposed moves are very small, how does this affect the acceptance probability? |think| How would this affect the evolution of the system?

|think| Similarly, what happens when proposed moves are very big?

|action| Navigate to 'inputs' :math:`\rightarrow` 'Tut_4' :math:`\rightarrow` 'extensions' :math:`\rightarrow` 'movesize' to find your standard DL_MONTE input files for this part of the tutorial.  If you open the CONTROL file, you will notice three new lines::

    maxatmdist   0.1               # Maximum atom displacement for a proposed move is 0.1 Angstroms
    acceptatmmoveupdate      100   # Adjust the maximum atom displacement every 100 steps
    acceptatmmoveratio    0.37     # The desired ratio of successful translational moves to all attempted translational moves

There are two key parts of code that are needed for this performance optimisation-generating the move:

.. code-block:: html
   :linenos:
   
   atm = random_number * natoms + 1                           # Randomly select a particle in the system
   delta_pos = (random_number - 0.5) * max_atm_displacement   # Change its position by moving it a random distance that is less than the maximum
   pos_new = pos_old + delta_pos                              # Define the new position of the particle 

and updating the move size: 

.. code-block:: html
   :linenos:
   
    do iter = 1 to max_iterations                                   # Start a 'do' loop that represents the entire simulation
    
        DO MONTE CARLO STUFF                                        # Use DL_MONTE to run the simulation as normal
    
        if mod(iter / accept_atm_move_update) == 0                  # Execute the following lines if the step number is divisible by acceptatmmoveupdate
        
            ratio = accepted_moves / attempted_moves                # The acceptance ratio at this point in the simulation
            
            if ratio > accept_atm_move_ratio                        # Execute the following line if the ratio is greater than acceptatmmoveratio
            
                max_atm_displacement = max_atm_displacement * 1.05  # Increase the maxatmdist value by a factor of 1.05
                
            else                                                    # Execute the following line if the ratio is less than acceptatmmoveratio
            
                max_atm_displacement = max_atm_displacement * 0.95  # Decrease the maxatmdist value by a factor of 0.95
                
            endif                                                   # End the if statement at line number = 9
            
        endif                                                       # End the if statement at line number = 5
        
    enddo                                                           # End the 'do' loop, i.e. end the calculation

The maximum displacement of an atom is controlled by the variable *max_atm_displacement*. The *max_atm_displacement* can not be known prior to the start of the simulation and the most suitable valuable changes as the simulation progresses. The acceptance ratio (ratio of accepted moves to all proposed moves) can determine the rate of equilibration and the efficiency of the sampling. For these reasons DL_MONTE provides a mechanism for adjusting the value of *max_atm_displacement* as the simulation proceeds.

The initial values in the CONTROL file are the default values for DL_MONTE but by altering these values you can improve the efficiency of sampling and minimise the equilbration time.

|action| Vary each of these values and investigate how the energy equilibrates during the course of the simulation. Try and determine the set of values that give the most efficient equilibration.

You can use::

   grep displacement OUTPUT.000

or the script::

   disp.sh

to print the initial values and the final value of the maximum displacement(s).

*N.B.*  This functionality should be used to identify the optimum move size for sampling a given system.  Beware using this functionality in a calculation as it can break detailed balance.

|think| How could the condition of detailed balance be broken by using this functionality?

2. Detailed balance for volume moves
------------------------------------

Establishing the condition for detailed balance in a simulation where volume moves are enabled is more complicated than for translational moves alone.  To maintain detailed balance with volume moves, the acceptance probability for a move from an initial configuration of particles in positions :math:`\mathbf{r}_1` in a volume, :math:`V_1` to a new configuration, :math:`\mathbf{r}_2` with a volume :math:`V_2`, changes in the Metropolis algorithm to:

.. math::

   P_{\mathrm{acc}}([\mathbf{r}_{1},V_1] \rightarrow [\mathbf{r}_2,V_2]) = \min(1, \exp \{- \beta [U(\mathbf{r}_2) - U(\mathbf{r}_1) + P_{ext}(V_{2}-V_{1}) - N \beta^{-1} \ln(V_{2} / V_{1}) ] \} )

where :math:`P_{ext}` is the external pressure acting on the system and :math:`\beta = \frac{1}{kT}`.  In most simulations, the positions of every object in a system is expressed as dimensionless, scalable position coordinates, which scale with the size of the system such that when the volume changes, the *relative* positions of the objects in the new size remains the same, but the distance between objects in the system goes up or down depending on whether the volume has increased or decreased.  However, the number of possible configurations of a system is determined in part by its total volume, such that a larger system will have more possible configurations.  This must be accounted with the :math:`N \beta^{-1} \ln(V_{2} / V_{1})` term.  The other terms in the exponent come from the probability distributions of isothermal-isobaric systems, where the :math:`P_{ext}(V_{2}-V_{1})` represents the work done *on* the system by an external pressure.  For more information on detailed balance for volume moves, see [#f2]_.

|think| When performing volume moves on molecular systems, the position of the centre of mass of the molecule is scaled, as opposed to the positions of all of its constituent atoms, rationalise this caveat. (Hint: what would happen to all the chemical bonds if the atom positions were scaled instead? How would this affect the likelihood of accepting the move?) 

3. Neighbour lists
------------------

DL_MONTE uses neighbour list to improve the performance of the energy calculation, particles only have to check interactions with their neighbours, not every particle in the simulation.  This is particularly beneficial when particles retain the same numbers for the whole simulation or for any attempted moves. In DL_MONTE, this functionality is described by the following lines in the CONTROL file::

   nbrlist auto

This rebuilds a particle's neighbourlist whenever necessary.  The size of the neighbour list is determined by::

   maxnonbondnbrs <int>

This determines the memory allocated for each particles neighbourlist.  The size will be determined by the size of your system, its density and the interaction cut-off as specified in the FIELD file.  

|action| Go to the directory named 'helloneighbour' and view the CONTROL file, you will notice that it looks identical to the CONTROL files that you have already seen, with one extra line called::

   verlet <float>  

|action| Run calculations using different values for this parameter and see how they affect the time taken to complete the calculation.  You can extract this for a given calculation by using the command::

   grep "total elapsed" OUTPUT.000

or alternatively the script::

   time.sh

|think| How does tuning the *verlet* parameter affect the duration of this calculation? |think| Why might this be the case?

4. Logarithmic Volume Moves
---------------------------

The 'linear' keyword in the 'volume move' line of the NPT CONTROL file represents how the volume will change, in this case, on a linear scale.  However, one can also set the volume change to a logarithmic scale, this can be more efficient in simulations where large volume changes are required to representatively sample configuration space. 

|action| Create a new directory and copy the CONFIG, CONTROL and FIELD files and the strip_yaml script from one of your completed calculations into it. Open the CONTROL file and change the keyword *linear* to  *log* in the volume move command.  This changes the way in which the the move is generated; logarithmic rather than linear, and is often claimed to be more efficient.  The acceptance criterion for the move in the Metropolis algorithm is now:

.. math::

  P_{\mathrm{acc}}([\mathbf{r}_{1},V_1] \rightarrow [\mathbf{r}_2,V_2]) = \min(0,  ( N + 1 ) \ln(\frac{V_{2}}{V_{1}}) - \beta [U(\mathbf{r}_2) - U(\mathbf{r}_1) + P_{ext}(V_{2}-V_{1})])
         
|action| Extract the time evolution of the volume for this calculation and compare it with the volume evolution from the equivalent linear calculation.  |think| Rationalise the observed differences.

.. rubric:: Footnotes

.. [#f1] B. L. Holian, "Shear viscosities away from the melting line: A comparison of equilibrium and nonequilibrium molecular dynamics", *J. Chem. Phys.*, **78**, 11, pp. 5147-5150, 1983.

.. [#f2] M. S. Shell, "Monte Carlo simulations in other ensembles"[online], University of California at Santa Barbara: Engineering, 2012.  Available from: https://engineering.ucsb.edu/~shell/che210d/Monte_Carlo_other_ensembles.pdf