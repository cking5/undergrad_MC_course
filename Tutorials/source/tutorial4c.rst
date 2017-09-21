.. _tutorial_4c: 

TUTORIAL 4: CH\ :sub:`4` \ adsorption in a zeolite
==================================================

Authors: Chris King, James Grant - r.j.grant@bath.ac.uk

Introduction
============

In the previous session, we introduced the Grand Canonical (GC) ensemble, where the total number of particles, *N*, is allowed to vary but the chemical potential, :math:`\mu`, the system volume, *V* and the temperature, *T* are fixed (:math:`\mu`\VT).  We explored the advantages of operating under the GC ensemble in accurately simulating the thermal behaviour of a material as opposed to using 'fixed *N*' ensembles.  In this session, we will again be using the :math:`\mu`\VT ensemble, but this time we will be applying it to adsorption of gas molecules into a zeolite.

A zeolite is a solid (typically silica or alumina) that has a porous crystal structure.  This creates channels and pores within the solid that are often large enough to allow small molecules like H \:sub:`2`\ O and CO\ :sub:`2` \ to enter the solid.  This gives zeolites a very large surface area for interactions with other molecular species.  This high surface area allows zeolites to function as catalysts and 'molecular sieves', where they can separate components of a mixture by molecule size and/or its affinity to adsorb onto the zeolite surface.  Any molecule that is too large simply won't be able to enter the zeolite and molecules that are too small or have low affinity for the surface will pass through the material with little adsorption.  Zeolites are also very resistant to environmental conditions due to their chemical inertness; they have high melting points, resist high pressures, are insoluble in water and many other solvents and do not readily undergo redox reactions.

There are a wide array of applications for zeolites, such as high-density storage of gases like hydrogen for fuel cells, separating mixtures like crude oil and the products of 'cracking', and 'carbon scrubbers' for greenhouse gas emissions from power stations.  For a more comprehensive overview of zeolite materials and applications, see [#f1]_, [#f2]_ and [#f3]_.

One potential application of zeolites that has been the subject of much research is use as 'carbon scrubbers', where greenhouse gases like CO\ :sub:`2` \ and CH\ :sub:`4` \ are removed from the flue gases from power stations before they are released into the atmosphere.  Carbon scrubbers are just one of the many ways that greenhouse gas emissions can be reduced thereby limiting the effects of climate change. However, more work is needed to find the optimal zeolite structure that can preferentially adsorb, and therefore separate each greenhouse gas (or any other gas of interest) from a gaseous mixture.  This can be done relatively quickly by using computational modelling using different (and often hypothetical) zeolite structures and determining how much gas is adsorbed as the partial pressure and temperature change.  This is where Grand Canonical Monte Carlo (GCMC) techniques are particularly useful. 

In this session, we will be conducting a GCMC simulation of a bulk zeolite structure containing a variable number of CH\ :sub:`4` \ molecules. In this simulation, our Monte Carlo moves will be insert/delete moves which insert and remove CH\ :sub:`4` \ molecules from the system.  Unlike the previous session, we will allow translational moves for CH\ :sub:`4` \ so that each molecule can explore its local environment.  We will also define a new move type known as rotation moves, these are similar to translation moves but instead propose a rotation by a random value between zero and a pre-determined maximum number of degrees.  This move type only applies to molecules or other objects that lack full rotational symmetry, *i.e.* Lennard-Jones particles  used in previous sessions have been individual featureless spheres, which are rotationally symmetrical, so rotation moves would not change the configurational energy of the system.

First, we shall give a breakdown of the DL_MONTE input files for this system, now that we have real molecules and structures to deal with in our system.  An example of each one are detailed below:

CONFIG
------

As always the CONFIG file contains the starting structure:

.. code-block:: html
   :linenos:

   Zeolite (Si and O, with some Xe)
   0    0
        24.4750735      0.0000000      0.0000000
         0.0000000     24.4750735      0.0000000
         0.0000000      0.0000000     24.4750735
   NUMMOL       1       1     200
   MOLECULE zeolite     584    584
    Si       c
         0.4513898     -0.3668470     -0.4576217     0
    Si       c
        -0.1875927     -0.3694940     -0.4576217     0
   .....
   
As you can see, the CONFIG has the same format as the previous CONFIG files.  Our system is contained within a cube with dimensions of 24.4750735 Angstroms.  There is one molecule present: the 'zeolite' molecule, containing 584 atoms (up to a maximum of 584) which are either silicon, 'Si', oxygen, 'O\_', or xenon, 'Xe'.  

CONTROL
-------

The CONTROL file is shown below:

.. code-block:: html
   :linenos:

   GCMC simulation of CO2 in zeolite
   use gaspressure   		# use the partial pressure in GCMC moves (as opposed to chemical potential)
   use orthonormal              
   finish
   temperature       273.0
   acceptmolmoveupdate  200        # Period (in moves) at which the maximum move size is recalculated
   acceptmolrotupdate  200         # Period (in moves) at which the maximum rotation angle is updated
   steps             1000000        # Number of moves to perform in simulation
   equilibration     50000        # Equilibration time before statistics are gathered (in moves)
   print               1000        # Information is output every 'print' moves     
   revconformat dlmonte             # REVCON file is in DL_POLY CONFIG format
   stack              10000        # Size of blocks (in moves) for block averaging
   maxmolrot           0.005       # Initial maximum rotation angle (degrees) 
   move molecule 1 20              # Perform translation moves for 1 molecule type (ch4) 20% of the time
   ch4
   move gcinsertmol 1 60 0.5       # Perform insertion/removal moves for ch4 60% of the time, with a min. distance of 0.5 from atoms for inserts 
   ch4 0.0001                     # Use a partial pressure of 0.0001 (katm) for ch4
   start

The CONTROL looks a little different to what you're used to, this is primarily because we are now trying to move and insert/delete real molecules in our simulation, rather than simple spherical particles as we have used in previous sessions.  The 'use gaspressure' directive specified at the beginning of the CONTROL file means that the partial pressure of the gas, rather than the activity are specified.

.. math::

    a = \gamma \frac{P}{P_0},

where *a* is the activity, :math:`\gamma` is the fraction of the component within the gaseous mixture and is assumed to be 1 in this case, as we are dealing with pure CH4, and *P*, :math:`P_0` the pressure and reference pressure respectively.

The activity relates to chemical potential according to

.. math::

    a = \exp(\frac{\mu - \mu_0}{RT})

where :math:`\mu` and :math:`\mu_0` are the chemical potential and reference chemical potential (usually that of an ideal gas), *R* gas constant and *T* temperature.  

The 'use orthonormal' directive tells DL_MONTE to keep our coordinates in each dimension (x, y and z) 90\ :sup:`o` \ from each other. Lines 6 and 7 state how often to update the maximum move distance for translational moves and maximum rotation angle for rotation moves, respectively.  Lines 8-13 have the same function  as in the previous CONTROL files.  'maxmolrot' states the initial maximum rotation angle for CH\ :sub:`4` \ in the system.  The four lines proceeding this line define the translational and rotation moves for CH\ :sub:`4`, the first number states how many molecule types the move applies to and the second number states the relative weight at which the moves are conducted.  'move gcinsertmol' defines the insert/delete moves for CH\ :sub:`4`, it applies to just the one (CH\ :sub:`4`) molecule type with a weight of 60 like the other move types specified.  The third number defines the minimum distance that you can insert a CH\ :sub:`4` molecule from any other atoms already present in the system, any insertions below this distance are automatically rejected moves.  The final line states the partial pressure of CH\ :sub:`4`.

The CH\ :sub:`4` \ molecules are considered to be rigid during the simulation, this restriction typically has to be in place for standard GCMC in order to satisfy detailed balance.

FIELD
-----

The FIELD file is shown below:

.. code-block:: html
   :linenos:

   Force fields and bond constraints for for CH4 in a zeolite
   CUTOFF 12.0
   UNITS kcal
   NCONFIGS 1
   ATOMS 4
   Si core 28 0.0
   O_ core 16 0.0
   CH core 16 0.0
   Xe core 1 0.0
   MOLTYPES 2                                  # There are two molecules present in this system
   zeolite                                     # The first molecule is the zeolite 
   MAXATOM 584                                 # with (a maximum of) 584 atoms
   ch4                                         # The second molecule is methane
   ATOMS 1 1                                   # 1 atom type with a maximum of 1 atom in the molecule (?!)
   CH core  0.00000000 0.0000000 0.0000000     # 1 CH 'atom' positioned at the origin of the molecule
   FINISH
   VDW       4
   CH core       CH core       lj    0.31494  3.72
   O_ core       CH core       lj    0.224466  3.3765
   CH core      Xe core      12-6   16777216 0.0
   CH core      Si core      12-6   16777216 0.0
   CLOSE
  
The cutoff distance in this system is 12 Angstroms and the units of energy are in kcal.  There are four atom types: silicon, Si, atoms with mass of 28 amu, oxygen, O\_, atoms with mass = 16 amu, CH 'atoms' with mass = 16 amu and xenon, Xe, with mass = 1 amu.  All atoms have no net charge for the sake of simplicity.  As you may have noticed, the mass of the Xe atoms is not the same as its atomic mass because, in this simulation, the actual mass of Xe has no impact on the course of the simulation.

You will have also noticed that the methane molecules only have one CH 'atom', which might be unexpected given that a methane molecule actually contains one carbon and four hydrogen atoms with four C-H single covalent bonds.  This alternative description is used because, in computational simulations, calculations should be as efficient as possible.  One way of doing this is to reduce the system to the simplest representation possible while attempting to retain as much accuracy in the results as possible.  Consider the CH\ :sub:`4` \ molecule: a heteroatomic, tetrahedral, spherically-symmetrical molecule, containing (roughly speaking) non-polar C-H bonds.  This means that it has no net dipole moment and can be adequately described by one CH unit or 'atom' with the molecular mass of CH\ :sub:`4` \; 16.  This approximation of the full CH\ :sub:`4` \ structure and bonding is adequate for the purposes of this simulation.  More intuitive representations of CH\ :sub:`4` \ that more accurately describe CH\ :sub:`4` \ behaviour and properties exist, but these would add unnecessary complexity to our simulation.

At the end of the FIELD file, there are four defined interactions: one between two CH\ :sub:`4` \ molecules, one between CH\ :sub:`4` \ and the oxygen atoms in the zeolite, one between CH\ :sub:`4` \ and xenon and the final one between CH\ :sub:`4` \ and the silicon atoms in the zeolite.  You will see two different interaction types: the familiar 'lj' potential and the '12-6' potential.  12-6 is the name given to an alternative form of the Lennard-Jones potential:

.. math::

   \phi(r_{ij}) = \frac{A}{r_{ij}^{12}} - \frac{B}{r_{ij}^6}

where :math:`\phi(r_{ij})` is the potential energy between two particles, i and j, separated by a distance, :math:`r_{ij}`, *A* and *B* are constants. The first term therefore represent the repulsive part of the Lennard-Jones potential and the second term represents the attractive part of the potential. The two numbers specified in the lines for the '12-6' interactions are *A* and *B*, respectively.  For more information, please refer to the DL_POLY manual. 

|think| By visualising the structure, or otherwise, identify why the zeolite contains Xe atoms. 

.. |think| image:: images/General/think.png
   :height: 100 px
   :scale: 25 %

HINT: The zeolite contains two different-sized pores in its unit cell, and experiments show that only one of these is involved in gas adsorption.  

Exercise 1)
===========

In this exercise, you will be running simulations of the zeolite solid with the potential to add/remove CH\ :sub:`4` \ over the course of the calculations.  Each of these calculations will be run at a constant temperature but with increasing partial pressure of CH\ :sub:`4`.  From the output of these calculations, you will be able to plot an *adsorption isotherm* of CH\ :sub:`4` \ in this zeolite.  An adsorption isotherm is a graph of the amount of gas adsorbed onto a surface plotted against partial pressure of the gas.  These are used to find the partial pressure at which maximum adsorption is obtained (the saturation pressure).

|action| Navigate to 'inputs' :math:`\rightarrow` 'Tut_5' :math:`\rightarrow` 'main' :math:`\rightarrow` 'init'.  You should find everything you need to run the calculation and perform the subsequent data analysis.  You will **not** need to find the equilibration time for this system at this temperature, it has been given in the CONTROL file already.  

.. |action| image:: images/General/action.png
   :scale: 5 %

|action| Run the calculation for the first value of the partial pressure stated in the CONTROL file.  It should take around *x* minutes to complete.  

|action| Once the calculation is complete, open the OUTPUT.000 file and note the average number and fluctuations of CH\ :sub:`4` \ which can be found near the very end of the OUTPUT.000 file: **script to extract these values instead?**

.. code-block:: html

   ...

   zeolite           1.0000          0.0000

   ch4               average number  fluctuations about the average


   ---------------------------------------------------------------------------------------
                            final energies
   ...

|action| Extract a time sequence of molecule numbers over the course of the simulation by running the script::

  strip_adsorb.sh

This script extracts the number of CH\ :sub:`4` \ from the OUTPUT.000 file each time it records a measurement from the simulation along with the number of steps that have elapsed and places these into a new file called 'adsorb.dat'.

|action| Plot the contents of 'adsorb.dat'.  |think| Is this what you might expect given the conditions specified in the CONTROL file?

|action| Re-run the calculations, increasing the partial pressure in the CONTROL file for each calculation until you reach the saturation pressure, :math:`P_s`, for this temperature.

|action| Create an estimate of the adsorption isotherm by plotting the average number of adsorbed CH\ :sub:`4` \ molecules against partial pressure.  |think| From your graph, identify :math:`P_s`.

|action| Also plot the time sequence of the number of CH\ :sub:`4` \ for each of your calculations. |think| How does the shape of these plots change with increasing partial pressure.  Rationalise your observations.

|think| By looking at the time sequences, what do you need to consider to ensure the accuracy of your calcuation? HINT: Remember, the equilibration time in the CONTROL file tells DL_MONTE how much of the output data is used to calculate final averages.

Exercise 2)
===========

Now that you understand the procedure of estimating an adsorption isotherm from these simulations, this exercise will focus on obtaining isotherms for a range of temperatures to see how varying the temperature changes the adsorption behaviour of the zeolite.

|think| Consider how temperature may affect the number of molecules adsorbed onto the zeolite surface.

|action| First, use the following script which changes the temperature in the CONTROL file and copies the CONFIG, CONTROL and FIELD files into a new directory named after the new temperature value::

  **change_temp script**

*N.B.* You will not need to go above *y* K.

|action| Remember, you will now need to ensure that your system has had sufficient time to equilibrate, do this for each new temperature you use.

|action| Now estimate an adsorption isotherm for each of your chosen temperatures by repeating Exercise 1 and identify :math:`P_s` at each temperature.

|action| Plot :math:`P_s` against temperature.  |think| From this graph, what set of conditions gives the maximum adsorption?

|think| Given that this particular zeolite is thermally-stable up to around *z* K, are the conditions for maximum adsorption feasible?

Conclusions:
============

In this session, you will have appreciated the application of GCMC in the wider context of computational chemistry research and used GCMC to model the adsorption properties of methane onto a siliceous zeolite.  You will have considered how to modify the simulation to improve the accuracy of the results.  You will also have compared the results of your model with those from experiments and thence considered ways to improve upon the existing model.  The next and final session of this course will encourage you to apply all that you have learned in this and previous sessions to solve problems.

.. rubric:: Footnotes

.. [#f1] ed. Cejka, J., Herman, H.V., Corma, A., Schuth, F., *Introduction to Zeolite Science and Practice*, Elsevier Science, Burlington, 2007, **168**, 1-1058.
.. [#f2] Breck, D. W., *Zeolite molecular sieves: structure, chemistry, and use*, Wiley, 1973.
.. [#f3] Chester, A. W., *Zeolite Chemistry and Catalysis*, Springer Netherlands, Dordrecht, 2009.


