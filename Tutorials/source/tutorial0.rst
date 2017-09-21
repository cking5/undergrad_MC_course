.. _tutorial_0:

TUTORIAL 0: Introduction to Molecular Dynamics
==============================================

Authors: Chris King, James Grant - r.j.grant@bath.ac.uk

Molecular Dynamics (MD) is a classical potential modelling simulation technique.  Potential models employ several approximations compared to quantum mechanical modelling, such as: integrating electrons out of calculations, all charges are point charges and describing interatomic (intermolecular) interactions as potential functions.  MD uses these potential functions and their associated force fields to plot trajectories of mobile objects by integrating Newton's Equations of Motion (EOMs).  This allows the velocities of any object at any given instant in the simulation to be calculated.  The time interval at which velocities are recalculated depends on the timescale over which interactions occur.  Molecular dynamics is a particularly useful way to model liquids and gases as well as predict transport properties of chemical species in solid and/or liquid solutions, such as diffusion, electrical conductivity and viscocity.  In this session, you will use molecular dynamics to simulate the phase behaviour of a model system.

Part 1: The Lennard-Jones Solid
===============================

Now that you've become acquainted with the nature of MD simulations, we shall now apply some of the concepts and methods to explore the properties of a model system: a 3D Lennard-Jones material.  A Lennard-Jones material is a hypothetical material which is made up of chargeless particles that only interact with each other via the Lennard-Jones potential.  The Lennard-Jones potential has the following form:

.. math::

  \phi(r_{ij}) = 4\epsilon\biggl[\Bigl(\frac{\sigma}{r_{ij}}\Bigr)^{12}-\Bigl(\frac{\sigma}{r_{ij}}\Bigr)^{6}\biggr]
         
where :math:`\phi(r_{ij})` is the potential energy of the interaction between particles i and j separated by a distance :math:`r_{ij}`, :math:`\epsilon` is the parameter relating to the depth of the potential well and :math:`\sigma` is the distance at which there is zero interaction.  The positive term represents the repulsive component of the interactions which stems from phenomena like the Pauli Exclusion Principle and is dominant at short distances.  The negative term represents the attractive component of the interaction and stems from phenomena like London Dispersion and dominates at longer distances.  The figure below illustrates the shape of :math:`\phi(r_{ij})` and the corresponding force acting between two objects, :math:`F(r_{ij}) = -\frac{\mathrm{d}\phi(r_{ij})}{\mathrm{d}r_{ij}}`, as a function of :math:`r_{ij}`.  When :math:`F(r_{ij})` is positive, the resultant force acting between the two objects is repulsive and when :math:`F(r_{ij})` is negative, the resultant force are attractive.  As you can see, at :math:`r_{ij}<r_{min}`, :math:`F(r_{ij})` is repulsive and the two objects repel each other, and when :math:`r_{ij} <  r_{min}`, :math:`F(r_{ij})` is attractive and the two objects move towards each other.

.. figure:: images/Tut_0_images/LJ_potential.png
   :align: center

   **Figure 1:** Shapes of (a) the Lennard-Jones potential, :math:`\phi(r_{ij})`, between two particles and (b) the corresponding force, :math:`F(r_{ij})`, acting between the two particles.

This model system has limited use in real world applications as it neglects all other possible intermolecular interactions, such as the Coulombic, permanent dipole-dipole and hydrogen bonding.  However, Lennard-Jones solids have been used to model the phase behaviour of noble gases, like Argon, with a reasonable degree of accuracy [#f1]_, [#f2]_, [#f3]_.

Exercise 1)
-----------

In this exercise, we will be running molecular dynamics calculations on a 3D model of a Lennard-Jones solid in an face-centred cubic (fcc) arrangement to predict its behaviour when the temperature is raised.  

|action| Navigate to 'inputs' :math:`\rightarrow` 'Tut_1' :math:`\rightarrow` 'main' :math:`\rightarrow` 'Init'.  You will see input files: CONFIG, CONTROL, and FIELD.  The CONFIG file displays the size of the system and the position and velocity components in each direction (x, y, z) of every particle at the beginning of the simulation, the FIELD file defines the interactions between particles, and the CONTROL file sets the parameters and conditions for running the simulation. 

.. |action| image:: images/General/action.png
   :scale: 5 %

**(instructions on how to run the simulation)**  

As the calculation runs and completes, you will notice several new files appear in your directory.  These are standard output files for DLPOLY.  The files that you will be using in this exercise are the HISTORY, OUTPUT and REVCON files.  The HISTORY file contains the configuration of the system at different times within the simulation (the spacing of these intervals is specified in the 'print' line of the CONTROL file), OUTPUT records various properties like energies, pressure and temperature over the course of the simulation (at intervals specified by the CONTROL file) and REVCON displays the final configuration of the system.  

*N.B.* you may notice that the REVCON is bigger than the initial CONFIG file, this is because the initial CONFIG doesn't contain any velocity data (it assumes that the initial velocities of all particles are zero).

Once the calculation has successfully completed, it is now time to begin the analysis of the output data.  For this exercise, we begin by examining both the initial and final configurations of your system and seeing what differences, if any, are present.  

**(alternate instructions for visualising configurations here)**

|action| Click the windows button and search for a program called \'Vesta\' and open it.  Vesta is a program which allows you to visualise systems of atoms/molecules in three dimensions. 

|action| To view your initial and final system configurations, go to file :math:`\rightarrow` New Structure; select import and browse to find your CONFIG and REVCON, respectively.  Do this separately for each file.  You should see both structures appear in separate tabs. 

|think| Describe the initial and final configurations.  Rationalise any observed differences.

.. |think| image:: images/General/think.png
   :height: 100 px
   :scale: 25 %

**alternate instructions for visualising system evolution here** 

**use script to change temperature here**

|action| Repeat the calculation at increasingly high temperatures, following the instructions above, but changing the temperature value in the CONTROL file.  You will not need to go above 10 K.  You may wish to create a new directory for each temperature and copy the CONFIG, CONTROL and FIELD files into each.  

|action| View the REVCON from each calculation in Vesta (the CONFIG file will be the same for each one) and view the evolution of the system in VMD.  

|think| What do you notice about the final configuration of the system as the temperature increases? 

|think| Qualitatively determine and record the temperature(s) at which any significant transitions occur.  

Part 2: Energy in Molecular Dynamics
====================================

This part of the tutorial aims to help solidify your understanding of how kinetic and potential energy are treated in molecular dynamics and help to monitor the simulation.  The total energy of any thermodynamic system, *E*, can be broken down into the contributions from both kinetic, *KE*, and potential energy, *U*, such that:

.. math::

  E = U + KE
	
The conservation of total energy (*E* = constant) is critical to maintaining physicality of the system.  So if *KE* decreases, *U* must increase to keep *E* constant and vice-versa.  According to Kinetic Theory, the kinetic energy is directly proportional to the mean square speed of our particles, which in turn defines the temperature of the system:

.. math::

  KE = \frac{1}{2}m\langle c^{2} \rangle = \frac{3}{2}RT

where *m* is the total mass of all the particles, *R* is the molar gas constant, and *c* is the speed of the particle (in an arbitrary direction), the <...> represent taking the average value of the variable inside them.  In this case, the average is conducted over all particles in the system.

In our model (and many other classical models), the total potential energy of the system is the sum of the interaction energies of each particle with the rest of the system: 

.. math::

  U = \sum_{i} \psi_{i}

where:

.. math::

  \psi_i = \sum_{j=1,j \neq i}^{N-1} \phi(r_{ij})

where *N* is the total number of particles in a system and :math:`\psi_i` is the total interaction energy of particle i with all other particles in the system (excluding itself). 
 
The Lennard-Jones potential represents a short-range interaction (:math:`r_{ij}^{-6}` and :math:`r_{ij}^{-12}`), the contribution to the total interaction becomes infinitesimal as particles become further apart.  Also, the calculation time increases considerably if we explicitly calculate the interaction energy for each particle pair, so it is common to often invoke a cut-off distance.  By convention, this is taken as 2.5 :math:`\sigma` and is stated in the CONTROL file of the simulation. For a given particle, only particles within the cut-off are assumed to significantly contribute to the interaction energy.  This introduces a small but easily-correctable error in our calculated values. 

Exercise 2)
-----------

In this part of the tutorial, we will extract total, potential and kinetic energies of the system from the OUTPUT file and plot them as a function of temperature.  

|action| start by navigating to one of your directories (in the Command Prompt) and run the following command:

**script like analysis.sh here**

This will activate a script which will extract *T*, *E*, and *U*, from the OUTPUT file and place them into a new file called 'output'.  It also calculates :math:`KE = E - U` and appends it to 'output'.  

|action| Run this command on each of your simulations so that you have a data file in each of your repositories.  

|action| Plot *E*, *U* and *KE* against *T* on the same graph, using whichever program you're most comfortable with (Excel, MATLab, gnuplot etc.).  It may also be helpful to run more simulations around the transition temperature.  

|action| Comment on the shape of the plots.  |think| Do these indicate the presence of a phase transition?

Part 3: Cooling in Molecular Dynamics
=====================================

As you have seen from the tutorial so far, potential modelling of physical systems can reliably and accurately simulate the thermodynamic behaviour when increasing the temperature.  However, it can be a lot harder to cool a system back down in a way that reflects observed behaviour of real materials.

Exercise 3)
-----------

In this final exercise, you will observe what happens when you cool your Lennard-Jones fluid.  

|action| Go to a directory where the simulation has *just* melted (*i.e.* at a temperature just above the estimated melting point) and copy the REVCON, CONTROL and FIELD files into a new directory.  

|action| Rename REVCON to CONFIG and change the temperature in the CONTROL file to a value *just* below your system's melting point.  Now you should have everything ready to simulate the cooling of your liquid back into a solid.  We take the REVCON and not the CONFIG as we want the final melted configuration from the 'hot' simulation to be the starting configuration in the 'cool' simulation.  

|action| Run the simulation and record your observations.  |think| Is this what you expect given your knowledge of thermodynamics?

It is far more difficult to accurately model a system's thermodynamic behaviour when reducing the temperature using MD (or any potential modelling technique) primarily because of entropy, *S*, and the Third Law of Thermodynamics.  The Third Law of Thermodynamics can be stated as :math:`S \geq 0`.  When you cool a system, its entropy decreases, but this corresponds to an increase in entropy of its surroundings such that the Third Law of Thermodynamics is obeyed.  In a computational simulation, it is difficult to define 'entropy' and 'the surroundings' in this way, so when you cool a system from a temperature where it is liquid to one where it is solid, the observed 'disorder' of a system will not change, and the system will still appear to be liquid (or it may become a glass, if you run for long enough times).  Also, a system crystallises when the atoms within the system enter into a fixed orientation relative to one another, if all the atoms are freely moving, this outcome is **highly** unlikely.

In this simulation, we have been modelling the solid-liquid phase transition of our Lennard-Jones material primarily because we have operated under the constraint that the volume of the system is constant and the volume change between the solid and liquid phases is small compared to the solid-gas and liquid-gas volume change.  Our system volume is slightly larger than is required for the solid state to form so that the phase transition to liquid can be readily observed, but this also means that when trying to freeze the liquid back into a more condensed solid is more difficult. Trying to re-create the more condensed solid from the liquid in the expanded volume creates an additional energy barrier that needs to be overcome before freezing can occur. 

Conclusions:
============

In this session, you have been introduced to the potential modelling technique, Molecular Dynamics (MD).  You should now be aware of the approximations employed by potential models and how MD can be used to calculate useful properties in dynamic systems.  You have illustrated the use of MD to simulate a model system of a Lennard-Jones solid to observe its thermodynamic behaviour as you change its temperature and compared it to the behaviour of real systems.  By the end of this session, you should have:

- determined a phase transition, both qualitatively from the time-evolution of the system and more quantitatively from plots of system energies
- seen how potential modelling techniques deal with thermodynamic quantities like energy, entropy and particle trajectories 
- appreciated the limitations of such techniques in recovering the full range of observed thermal behaviour of real-life systems

Now that you have an awareness of MD techniques, we will move onto introducing the general theory and methodology of Monte Carlo simulations.

Extensions (optional)
=====================

1. Latent Heat
--------------

In your studies you may have come across the idea of latent heat of phase transitions.  Latent heat, *L*, can be described as the energy required for all particles in a material to overcome thermal activation barriers and become more mobile in a less condensed phase (solid-liquid, liquid-gas).  This is observed as a plateau at the transition temperature of heating curves, where no change in temperature is seen despite heat flowing into the system, or as a step-change in the potential energy at the phase transition as a function of temperature.  |think| From your plot of *U* vs *T*, estimate the latent heat for the solid-liquid phase transition of the Lennard-Jones material.

2. Ehrenfest classification:
----------------------------

A widely-used classification of phase transitions is the Ehrenfest classification, which describes phase transitions as n\ :sup:`th` \ order, where n is the n\ :sup:`th` \ order temperature derivative of an intrinsic quantity where a discontinuity occurs (see Figure 2).  For instance, the liquid-gas phase transition is described as a 1\ :sup:`st` \ order phase transition as there is a discontinuity in :math:`C_{v} = \frac{\partial U}{\partial T}`.  While a solid-solid phase transition is a 2\ :sup:`nd` \ order phase transition as there is a discontinuity in :math:`\frac{\partial C_{v}}{\partial T} = \frac{\partial^{2} U}{\partial T^{2}}`.

.. figure:: images/Tut_0_images/Ehrenfest.png
   :align: center

   **Figure 2:** Gibbs Free Energy, *G*, volume, *V*, enthalpy, *H*, entropy, *S*, and heat capacity at constant pressure, :math:`C_{p}` graphs against temperature for 0\ :sup:`th`\, 1\ :sup:`st` \ and 2\ :sup:`nd` \ order Ehrenfest phase transitions..

|think| With this in mind, what type of phase transition is your Lennard-Jones system undergoing and why?

.. rubric:: Footnotes

.. [#f1] W. T. Ashurst and W. G. Hoover, "Argon Shear Viscosity via a Lennard-Jones Potential with Equilibrium and Nonequilibrium Molecular Dynamics", *Phys. Rev. Lett.*, **31**, 4, 206-208, July 1973.
.. [#F2] B. W. Davies, "Radial Distribution Function for Argon: Calculations from Thermodynamic Properties and the Lennard-Jones 6:12 Potential", *J. Chem. Phys.*, **54**, 11, pp.4616-4625, June 1971. 
.. [#F3] R. O. Watts, "Percus-Yevick Approximation for the Truncated Lennard-Jones (12, 6) Potential Applied to Argon", *J. Chem. Phys.*, **50**, 2, pp. 984-988, January 1969.  