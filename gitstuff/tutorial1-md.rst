.. _tutorial_1:

TUTORIAL 1 : NVT Lennard-Jones solid
====================================

Authors James Grant, John Purton, r.j.grant@bath.ac.uk

Introduction to Molecular Dynamics
-----------------------

Molecular Dynamics (MD) is a classical potential modelling simulation technique.  Potential models employ several approximations compared to quantum mechanical modelling, such as: integrating electrons out of calculations, all charges are point charges and describing interatomic (intermolecular) interactions as potential functions.  MD uses these potential functions aand their associated force fields to plot trajectories of mobile objects by integrating Newton's Equations of Motion (EOMs).  This allows the velocities of any object at any given instant in the simulation to be calculated.  The time interval at which velocities are recalculated depends on the timescale over which interactions occur.  Molecular dynamics is a particularly useful way to model liquids and gases as well as predict transport properties of chemical species in solid and/or liquid solutions, such as diffusion, electrical conductivity and viscocity.  In this session, you will use molecular dynamics to simulate the phase behaviour of a model system.

The Lennard-Jones Solid
-----------------------

Now that you’ve become acquainted with the nature of MD simulations, we shall now apply some of the concepts and methods to explore the properties of a model system: a 2D Lennard-Jones material.
A Lennard-Jones material is a hypothetical material (LennardJonesium) which is made up of discrete chargeless particles with mass = 1 au that only interact with each other via the Lennard-Jones potential.  The Lennard-Jones potential has the following form:

.. math::

         \phi(r_{ij}) = 4\epsilon(((\frac{\sigma}{r_{ij})^{12})-((\frac{\sigma}{r_{ij}})^{6}))
         
where *\phi(r_{ij})* is the potential energy of the interaction between particles i and j separated by a distance *r_{ij}*, *\epsilon* is parameter relating to the depth of the potential well and *\sigma* is finite distance at which the interaction is zero.  The positive term represents the repulsive component of the interactions which stems from phenomena like the Pauli Exclusion Principle and is dominant at short distances.  The negative term represents the attractive component of the interaction and stems from phenomena like London Dispersion and dominates at longer distances.  The figure below illustrates the shape of *\phi_{ij}* and the corresponding force acting between two objects, :math:`F_{ij} = - \frac{d\phi_{ij}}{dr_{ij}}`, as a function of *r_{ij}*.  When *F_{ij}* is positive, the resultant force acting between the two objects is repulsive and when *F_{ij}* is negative, the resultant force are attractive.  As you can see, at :math: `r_{ij} <  r_{min}`, *F_{ij}* is repulsive and the two objects repel each other, and when :math: `r_{ij} <  r_{min}`, *F_{ij}* is attractive and the two objects move towards each other.

.. figure:: Tut_0_images/LJ_potential.png

This model system has limited use in real world applications as it neglects all other possible intermolecular interactions, such as the Coulombic, permanent dipole-dipole and hydrogen bonding.  However, LennardJonesium has been used to model the phase behaviour of noble gases, like argon, with a reasonable degree of accuracy [#F1], [#F2], [#F3].

Exercise 1)
-----------

In this exercise, we will be running molecular dynamics calculations on a 3D model of solid LennardJonesium in an face-centred cubic (fcc) arrangement to predict its behaviour when the temperature is raised.  

Navigate to ________________.  You will see input files: CONFIG, CONTROL, and FIELD.  The CONFIG file displays the size of the system and the position and velocity components in each direction (x, y, z) of every particle at the beginning of the simulation, the FIELD file defines the interactions between particles, and the CONTROL file sets parameters and conditions for running the simulation. Create a new folder (directory) in your area and copy the files into it.  You should also open these files to familiarise yourself with their contents.
Now you will attempt to run the simulation using these input files.  Open the Command Prompt on your Windows machine and type the following command::

	___________

and press 'Enter' on your keyboard.  This calls a molecular dynamics program called DL_POLY to run the simulation according to the parameters in the CONFIG, CONTROL and FIELD files.  The calculation should take around __ minutes to complete, you will know when it is complete when your working directory is displayed next to the keyboard prompt, you can also look at the latest modification times on your output files.  If the calculation finishes almost immediately after submitting it, an error has occurred.  In this instance, ask a demonstrator to help you find the problem (HINT: start by looking at the bottom of the OUTPUT file).  

As the calculation runs and completes, you will notice several new files appear in your directory.  These are standard output files for DLPOLY.  The files that you will need for the purposes of this exercise are the HISTORY, OUTPUT and REVCON files.  The HISTORY file contains the configuration of the system at different times within the simulation (the spacing of these intervals is specified in the 'print' line of the CONTROL file), OUTPUT records various properties like energies, pressure and temperature over the course of the simulation (at intervals specified by the CONTROL file) and REVCON displays the final configuration of the system.  

N.B. you may notice that the REVCON is bigger than the initial CONFIG file, this is because the initial CONFIG doesn't contain any velocity data (it assumes that the initial velocities of all particles are zero).

So you've successfully run a simulation and obtained some outputs, now what? Well, there are many possible actions one could take, depending on what you are looking for.  For this exercise, we will start with viewing both the initial and final configurations of your system and seeing what differences, if any, are present.  To do this, click the windows button and search for a program called 'Vesta' and open it.  Vesta is a program which allows you to visualise systems of atoms/molecules in three dimensions. To view your initial and final system configurations, go to file -> New Structure; select import and browse to find your CONFIG and REVCON, respectively.  Do this separately for each file.  You should see both structures in separate tabs. 

Describe the initial and final configurations.  Rationalise any observed differences.

You can visualise the time evolution of the system using a program called VMD, find and open this program.  You should see several windows appear, an example of what you should see is shown below:

.. figure:: Tut_0_images/VMD_initialise.png

In the 'VMD Main' window, click file -> New Molecule, then click 'Browse' and navigate to your chosen HISTORY file and open it.  From the 'Determine file type' drop-down menu, select either 'DL_POLY_4 HISTORY' or 'DL_POLY_C HISTORY' and then press 'Load'.  You should see an animation appear in the display window, each frame is a system configuration from the HISTORY file.  It will initially display the particles as lines and won't appear very informative.  You can change this by navigating to: Graphics -> Representations and choosing from the 'Drawing Method' drop-down menu ('VDW' is probably the most intuitive way to visualise the system).  You can reduce the frame rate of the animation by adjusting the 'speed' scale in the 'VMD Main' window.  

Repeat the calculation at increasingly high temperatures, following the instructions above, but changing the temperature value in the CONTROL file.  You will not need to go above 10 K.  You may wish to create a new directory for each temperature and copy the CONFIG, CONTROL and FIELD files into each.  

Note, you can add '&' to the end of the run command to make the calculation run in the background, allowing you to use the Command Prompt to run the calculations from the other directories.  The calculations will then run simultaneously in the background, though running a lot calculations may result in some performance issues for your machine while they all run.  To avoid this, try to avoid running any more than __ calculations at one time.  

View the REVCON from each calculation in Vesta (the CONFIG file will be the same for each one) and view the evolution of the system in VMD.  What do you notice about the final configuration of the system as the temperature increases? What happens to the solid as the temperature is increased? Qualitatively determine and record the temperature(s) at which any significant transitions occur.  

N.B. You will only be able to reliably view one animation at a time in VMD, so you will either need to quit VMD (by closing the 'VMD Main' window) or by deleting your 'molecule' from the 'VMD Main' window by selecting the entry in the window, then selecting: 'Molecule' -> 'Delete Molecule'. 

Part 2: Energy in Molecular Dynamics Simulations
------------------------------------------------

This part of the tutorial aims to help solidify your understanding of how kinetic energy and potential energy are treated and used to control and monitor a molecular dynamics simulation.  The total energy of any thermodynamic system, *E*, can be broken down into the contributions from both kinetic, *KE*, and potential energy, *U*, such that:
.. math::

	E = U + KE
	
The conservation of total energy (*E* = constant) is critical to maintaining physicality of the system.  So if *KE* decreases, *U* must increase to keep *E* constant and vice-versa.  According to Kinetic Theory, the kinetic energy is directly proportional to the mean square speed of our particles, which in turn defines the temperature of the system:
.. math::

	KE = frac{1}{2}m\langle c^{2} \rangle = frac{3}{2}RT

where *m* is the total mass of all the particles, *R* is the molar gas constant, and *c* is the speed of the particle (in an arbitrary direction), the <…> represent taking the average value of the variable inside them.  In this case, the average is conducted over all particles.
For our model (and many other classical models), the total potential energy of the system is the sum of the potential energies of each particle with the rest of the system: 
..math::

	U = \sum_{i} \psi_{i}

where:
..math::
	\psi_i = \sum_{j=1,j\neq}^{N-1} \phi(r_{ij})

where *N* is the total number of particles in a system and *\psi_i* is the total interaction energy of particle i with all other particles in the system (excluding itself).  
The Lennard-Jones potential represents a short-range interaction (:math: `(r_{ij})^-6` and :math: `(r_{ij})^-12)`, the contributions from interactions between particles become infinitesimal the further away they are from each other.  Also, the calculation time increases considerably if we explicitly calculate the interaction energy for each particle pair, so it is common to often invoke a cut-off distance.  By convention, this is taken as 2.5*\sigma* and is stated in the CONTROL file of the simulation. For a given particle, only particles within the cut-off are assumed to significantly contribute to the interaction energy.  This introduces a small but easily-correctable error in our calculated values. 


Exercise 2)
-----------

In this part of the tutorial, we will extract total, potential and kinetic energies of the system from the OUTPUT file and plot them as a function of temperature.  To do this, start by navigating to one of your directories in the Command Prompt and run the following command::

	_______________

This will activate a script which will extract *T*, *E*, and *U*, from the OUTPUT file and place them into a new file called _____.  It also calculates the average kinetic energy as :math: `E - U' appends it to ______.  Run this command on each of your simulations so that you have a data file in each of your repositories.  Now plot *E*, *U* and *KE* against *T* on the same graph, using whichever program you're most comfortable with (Excel, MATLab, gnuplot etc.).  It may also be helpful to run more simulations around the transition temperature to improve the accuracy of your plotted data at the transition.  Comment on the shape of the plots.  Do these indicate the presence of a phase transition?

Part 3: Cooling in Molecular dynamics simulations
-------------------------------------------------

As you have seen from the tutorial so far, potential modelling of physical systems can reliably and accurately simulate the thermodynamic behaviour when increasing the temperature.  However, for reasons that we will discuss, it can be a lot harder to cool a system back down in a way that reflects observed physical behaviour.

Exercise 3)
-----------

In this final exercise, you will observe what happens when you cool your LennardJonesium liquid.  To do this, go to a directory where the simulation has *just* melted (*i.e.* at a temperature just above the estimated melting point) and copy the REVCON, CONTROL and FIELD files into a new directory.  Rename REVCON to CONFIG and change the temperature in the CONTROL file to a value *just* below your system's melting point.  Now you should have everything ready to simulate the cooling of your liquid back into a solid.  We take the REVCON and not the CONFIG as we want the final melted configuration from the 'hot' simulation to be the starting configuration in the 'cool' simulation.  Now run the simulation and view the results in both Vesta and VMD.  Record your observations.  Is this what you expect? Is this behaviour supported by thermodynamic theory?

It is far more difficult to accurately model a systems thermodynamic behaviour when reducing the temperature using MD (or any potential modelling technique) primarily because of entropy, *S*, and the Third Law of Thermodynamics.  The Third Law of Thermodynamics can be stated as :math:`S \geq 0`.  When you cool a system, its entropy decreases, but this corresponds to an increase in entropy of its surroundings such that the Third Law of Thermodynamics is obeyed.  In a computational simulation, it is difficult to define 'entropy' and 'the surroundings' in this way, so when you cool a system from a temperature where it is liquid to one where it is solid, the observed 'disorder' of a system will not change, and the system will still appear to be liquid (or it may become a glass, if you run for long enough times).  Also, a system crystallises when the atoms within the system to enter into a fixed orientation relative to one another, if all the atoms are freely moving, this outcome is **highly** unlikely.

Conclusions:
------------

Congratulations, you have applied molecular dynamics to a model system of LennardJonesium to observe its thermodynamic behaviour as you change its temperature and related it back to the behaviour of real-life systems.  You have determined a phase transition, both qualitatively from the time-evolution of the system and more quantitatively from plots of system energies.  You have seen how potential modelling techniques deal with thermodynamic quantities like energy, entropy and particle trajectories and the limitations of such techniques in recovering the full range of observed thermal behaviour of real-life systems.

Extensions (optional):
----------------------

In your studies you may have come across the idea of latent heat of phase transitions.  Latent heat, *L*, can be described as the energy required for all particles in a material to overcome thermal activation barriers and become more mobile in a less condensed phase (solid-liquid, liquid-gas).  This is observed as a plateau at the transition temperatures of heating curves, where no change in temperature is seen despite heat flowing into the system, or as a step-change in the potential energy at the phase transition as a function of temperature.  From your plot of *U* vs *T*, estimate the latent heat for the solid-liquid phase transition of LennardJonesium.

A widely-used classification of phase transitions is the Ehrenfest classification, which describes phase transitions as n^{th} order, where n is the n^{th} order temperature derivative of internal energy where a discontinuity occurs.  For instance, the liquid-gas phase transition is described as a 1^{st} order phase transition as there is a discontinuity in :math:`C_{v} = \frac{\partial U}}{\partial T}_{V}.  While a solid-solid phase transition is a 2^{nd} order phase transition as there is a discontinuity in :math:\fract{\partial C_{v}}{\partial T} = fract{(\partial)^{2} U}{\partial T^{2}}.  The figure below illustrates the behaviour of 0^{th}, 1^{st} and 2^{nd} order phase transitions in terms of: Gibbs Free Energy, *G*, volume, *V*, enthalpy, *H*, entropy, *S*, and heat capacity at constant pressure, *C_{p}*.

.. figure::Tut_0_images/Ehrenfest.png

With this in mind, what type of phase transition is your LennardJonesium system undergoing and why?

.. rubric:: Footnotes

.. [#f1] W. T. Ashurst and W. G. Hoover, "Argon Shear Viscosity via a Lennard-Jones Potential with Equilibrium and Nonequilibrium Molecular Dynamics", *Phys. Rev. Lett.*, 31, 4, 206-208, July 1973.
.. [#F2] B. W. Davies, "Radial Distribution Function for Argon: Calculations from Thermodynamic Properties and the Lennard-Jones 6:12 Potential", *J. Chem. Phys.*, 54, 11, pp.4616-4625, June 1971. 
.. [#F3] R. O. Watts, "Percus-Yevick Approximation for the Truncated Lennard-Jones (12, 6) Potential Applied to Argon", *J. Chem. Phys.*, 50, 2, pp. 984-988, January 1969.  