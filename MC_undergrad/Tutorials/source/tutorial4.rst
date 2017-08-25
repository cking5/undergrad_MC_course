.. _tutorial_4:

-------------------------------------
TUTORIAL 4: Grand Canonical Ensemble
-------------------------------------

Authors: Chris King, James Grant - r.j.grant@bath.ac.uk

Introduction
============

So far in this course, we have been running simulations on systems with a constant number of particles (NVT or NPT ensembles).  This is useful for studying systems where there is no exchange of particles from one defined type to another and where predicting system behaviour by optimising particle number, *N*, is not required.  However, not all systems can be accurately modelled in this way, for instance, phase coexistence, where more than one phase of a material exist in equilibrium.  Other examples include: adsorption onto a surface (as you will see in this and future sessions), effect of solid-state defects on bulk properties and dissolution rates.  In this session, we will explore simulations where the number of particles is not constant and particles can be readily added or removed from our system.

However, we know that ensembles require three parameters to be kept constant, so what replaces *N*? The answer is chemical potential, :math:`\mu`. There are several definitions of chemical potential in the literature, but the most common definition is the effect of number of particles on the total Gibbs Free Energy, *G*, of the system:

.. math::

   \mu = \Bigl|\frac{\partial G}{\partial N}\Bigr|_{V, T} 

Or in a multi-component system, the effect of system composition on *G*, where the i:sup:`th` component of the system has its own chemical potential, :math:`\mu_i`:

.. math::

   \mu_i = \Bigl|\frac{\partial G}{\partial N_i}\Bigr|_{V, T, N_{j, j \neq i}}

:math:`\mu` can be thought of as the 'chemical driving force' for a system to reach equilibrium (see Figure 1) where the equilibrium composition of a system is the composition with minimum Gibbs free energy.  For more information about the chemical potential, you can look in any general physical chemistry textbook or refer to some of these online resources [1] http://www.icsm.fr/Local/icsm/files/286/JFD_Chemical-potential.pdf, [2] http://chem.atmos.colostate.edu/AT620/Sonia_uploads/ATS620_F11_Lecture5/Lecture5_AT620_083111.pdf.

.. figure:: images/Tut_4_images/chemical_potential.png
   :align: center

   **Figure 1**: A graphical representation of the chemical potential of a two-component (A and B particles) system driving it towards its equilibrium composition.

The :math:`\mu`\VT ensemble is known as the Grand Canonical ensemble and is the most widely-used 'variable *N*' ensemble.  Other simulation techniques, like Molecular Dynamics, are very inefficient at accurately modelling systems where the total number of particles isn't constant, which is why Monte Carlo is the preferred technique in these situations.
In fact, the problems that we encountered in the previous session are readily overcome if we operate under the Grand Canonical ensemble.

|think| One could, in principle, also choose to use an :math:`\mu`\PT ensemble, what kind of problems could arise when running simulations under this ensemble?

.. |think| image:: images/General/think.png
   :height: 100 px
   :scale: 25 %

Removing the constraint of having fixed particle number allows several new Monte Carlo move types to be defined, the most fundamental being insert/delete moves, where the proposed move adds or removes a particle from the system, respectively.  Other move types like swaps and replacements can be defined in terms of insert/delete moves.

|think| Define a swap move, where a particle at one position is swapped with another particle at a different position, as a sequence of insert and delete moves.

|think| Define a replacement move, where a particle in one position is changed to a particle of a different type, but remains in the same position, as a sequence of insert and delete moves.

The total energy, :math:`E` of the system performing an insert/delete move is defined as:

.. math::

  E(\mathbf{r}_2,N_2) = U(\mathbf{r}_1,N_1) + \mu \Delta N

where :math:`U(\mathbf{r}_1)` is the energy of the initial configuration, :math:`\mathbf{r}_1`, of the system, :math:`N_1, N_2` are the initial and final numbers of particles, respectively, and :math:`\Delta N = N_2 - N_1`. |think| Is :math:`\Delta N` positive or negative for insertions? |think| What about deletions?

As in the previous sessions, we will be using DLMONTE run Monte Carlo calculations on the phase behaviour of our all-too-familiar Lennard-Jones material.  However, all of our calculations in this tutorial will be conducted under the Grand Canonical ensemble.  You will hopefully see that we can get a more accurate reflection of the phase behaviour of real systems than if we are restricted to either NVT or NPT ensembles.

Exercise 1)
===========

|action| Navigate to 'inputs' :math:`\rightarrow` 'Tut_4' :math:`\rightarrow` 'main' :math:`\rightarrow` 'Equil'.  You will see the standard DLMONTE inputs files: CONFIG, CONTROL and FIELD, as well as some scripts for use later.  If you open the input files, you will see several differences which reflect the variable nature of *N* in the simulation:

.. |action| image:: images/General/action.png
   :scale: 5 %

CONFIG
------

Below shows the general CONFIG file structure used in this tutorial:

.. code-block:: html
   :linenos:

   Lennard-Jones muVT; particles are molecules, not atoms       # Title
         0         1                                            # Integers describing how the input is read in and the style of coordinates, respectively
   10.0000000000000000  0.0000000000000000  0.0000000000000000  # These lines describe the dimensions of the system in terms of basis lattice vectors
   0.0000000000000000  10.0000000000000000  0.0000000000000000  # Since our system is 3D, we need three basis vectors to fully describe it
   0.0000000000000000  0.0000000000000000  10.0000000000000000  # In this case, the system is a cube with sides of length 10 Angstroms
   NUMMOL 8 1000                                                # Specifies the minimum and maximum number of molecules in the system.
   MOLECULE lj 1 1                                              # Molecule 'lj' has 1 atom in it and has a maximum of 1 atom in it
   LJ   core                                                    # Now each particle is read in to the file, in the form: NAME core
   -5.0000000000000000 -5.0000000000000000 -5.0000000000000000  # x y z position
   MOLECULE lj 1 1                                              # continues to define the remaining molecules in the system, in this case, we start with the minimum number: 8
   LJ   core
   0.0000000000000000 -5.0000000000000000 -5.0000000000000000 
   etc

You will notice that the CONFIG is much smaller than its counterpart used in the last session.  This is because the number of particles (or molecules in this case) will vary over the course of the simulation, we need to only specify the initial configuration, which will start with only 8 molecules.  In principle, you can define the locations of any number of molecules in the CONFIG file (as long as that number falls between the minimum and maximum numbers stated in the 'NUMMOL' line), but for the purposes of this tutorial, we start at the minimum number: 8.

CONTROL
-------

The CONTROL file will take the following form in this tutorial:

.. code-block:: html
   :linenos:
  
   GCMC Lennard-Jones              # Title
   finish                          # Needs to be here, some conditions must be placed before this word, but there aren't any in this case
   seeds 12 34 56 78               # Sets the initial configuration
   temperature     1.4283461511745 # Temperature of the system in Kelvin
   # nbrlist auto                  # Use a neighbour list to speed up energy calculations
   # maxnonbondnbrs 512            # Maximum number of neighbours in neighbour list
   steps          10000            # Number of moves to perform over the course of the simulation
   equilibration    0              # Equilibration period: statistics are gathered after this period
   print           1000            # Print statistics every 'print' moves to the output file
   stack           1000            # Number of moves over which average values are calculated
   sample coord   10000            # How often to print configurations to ARCHIVE.000
   revconformat dlmonte            # REVCON file (final configuration) is in DLMONTE CONFIG format
   archiveformat dlpoly4           # ARCHIVE.000/HISTORY.000/TRAJECTORY.000 format 
                                   # Sets the format for the ARCHIVE.000/HISTORY.000/TRAJECTORY.000 files, in this case: HISTORY.000 in DLPOLY4 style
   yamldata 1000                   # Creates a YAMLDATA.000 output file that records data every 1000 steps
   move gcinsertmol 1 100 0.7      # Perform insertion/removal moves for lj, a weight 100 with a min. distance of 0.7 from atoms
   lj  0.06177                     # Use an activity of 0.06177   
   #  move atom 1 512              # Move atoms with a weight of 512
   #  LJ core 
   #  move volume cubic linear 1   # Move volume, box is cubic, linear scaling with a weight of 1
   start                           # Tells DLMONTE to begin calculation

The directives that switch on the neighbour lists: *nbrlist* and *maxnonbondnbrs* have been suspended in this session.  This is because the computational cost of maintaining the list under :math:`\mu`\VT ensembles negate the benefits when calculating the energy.  We have also suspended: atom translation moves for simplicity (though there is nothing in principle wrong with allowing these types of moves), and volume moves since we work under a constant-volume ensemble.

In this calculation DLMONTE is using the activity *a* rather than the chemical potential :math:`\mu`, which are related according to: 

.. math::

  a = \exp \Bigl(\frac{\mu}{RT}\Bigr)

where *R* is the gas constant.

|think| From the activity value given in the above CONTROL file, what value of :math:`\mu` does this correspond to?

FIELD
-----

The FIELD file looks almost identical to the ones from the previous session:

.. code-block:: html
   :linenos:

   Lennard-Jones                  # Title
   CUTOFF 2.5                     # The maximum distance between two particles at which the potential energy is calculated
   UNITS internal                 # Set the units of energies, internal = 10 J mol^-1
   NCONFIGS 1                     # Number of configurations described in the CONFIG file
   ATOMS 1                        # Number of atom types in the system
   LJ core 1.0  0.0               # In this case there is one atom type called 'LJ' with mass = 1.0 and charge = 0.0
   MOLTYPES 1                     # Number of molecule types in the system...
   lj                             # ...called 'lj'...
   ATOMS 1 1                      # ...with a maximum number of 1 atom...
   LJ core 0.0 0.0 0.0            # ...positioned at the origin
   FINISH                         # Completes the list of atom and molecule types in the system
   VDW 1                          # The number of potentials present in the system, must be the same as the number of interactions defined
   LJ core  LJ core lj   1.0 1.0  # Defines the interaction between two LJ atoms as a Lennard-Jones (lj) potential with epsilon = 1.0 eV and sigma = 1.0 Angstroms
   CLOSE                          # This ends the FIELD file once all interaction are described

In the NVT and NPT cases all the particles were declared to be part of the same molecule, now each particle is a molecule in its own right.  This distinction is made to simplify the calculation under :math:`\mu`\VT ensembles.  In principle, atoms can be added or removed from a molecule however, for simplicity, we shall insert or delete whole molecules rather than parts of molecules.  Since we have a single Lennard-Jones particle in each molecule we simply position it at the origin of the molecule.

Remember, there must be correspondence between the CONFIG and FIELD files, *i.e.* the number of molecule and atom types should be the same in both files.

|action| Run the DLMONTE calculations as you have done in the previous session (quick reminder of how to do it). Extract the time-sequence of the number of particles in the system by using the following script::

  [user@node-sw-119 tut_4]  strip_gcmc.sh

|action| By plotting the time-evolution of *N* for each of your simulations, increase the number of steps to determine when the system reaches equilibrium.

*N.B.* You will see that the output files will be mostly unchanged, except the YAMLDATA, which displays the number of molecules present instead of energies.

Exercise 2)
===========

Now that you have determined the ideal simulation length, |action| open the 'GCMC' folder in the 'main' folder.  |action| Replace the number of steps in the CONTROL file  with the value that you obtained from exercise 1.  |action| Run simulations at various different temperatures and activities by varying the appropriate values in the CONTROL file.  |action| Ensure that the system has equilibrated for each of your calculations. 

|action| Plot the time-evolution of *N* for each of your simulations.  |think| What happens to the number of adsorbed particles over the course of the simulation as you vary the temperature and activity? 

|think| From your results and your own knowledge, how does the value of :math:`\mu` change the ease at which particles are:

 a) inserted
 b) deleted 

|think| For what range of temperatures and activities is there:

 a) adsorption
 b) no adsorption
 c) an adsorption-desorption equilibrium

|action| You can also create histograms of the number of particles adsorbed to the surface over the course of the simulation, once you have produced the time sequence, with the script::

  hist.sh nmol.dat j

where *j* is the width of each bin used to generate the histogram.  You must specify the value of *j* in the command.  Though you are free to vary *j*, it is recommended that you set :math:`j = 1`.  Feel free to explore the effect *j* has on the shape of your histogram.

|think| How does the shape of the histogram vary with:

 a) temperature
 b) activity  

Conclusions:
============

In this session, you have been introduced to the Grand Canonical (GC) ensemble, where the total number of particles in the system can vary but the chemical potential of the system remains constnt.  This is useful for systems with interfaces or exchange of one particle type for another and has no equivalent in deterministic simulation techniques like Molecular Dynamics.  You have applied the GC emsemble to the hypothetical problem of Lennard-Jones particles adsorbing onto a surface to determine the ideal conditions for adsorption/desorption.  In the next session, we will apply the GC ensemble to the physical system of methane adsorption onto the surface of a zeolite in order to predict the conditions for ideal adsorption.

Extensions (optional):
======================

1. Detailed balance in the Grand Canonical ensemble
---------------------------------------------------

Like with the inclusion of volume moves in the previous session, the conditions through which detailed balance is maintained when employing insert/delete moves in :math:`\mu'\VT ensemble must be altered, such that, for particle insertions, the acceptance probability in the Metropolis algorithm in moving from an initial configuration, :math:`mathbf{r}_1`, with :math:`N_1 = N` particles, to a final configuration, :math:`mathbf{r}_2', with :math:`N_2 = N + 1` particles is:

.. math::
  
   P_{\mathrm{acc}}([\mathbf{r}_1,N_1] \rightarrow [\mathbf{r}_2,N_2] ) = \min(1,  \frac{V\Lambda^{-3}}{N+1} \exp \{- \beta [E(\mathbf{r}_2,N_2) - E(\mathbf{r}_1,N_1)] \} )

where :math:`V` is the system volume, :math:`\Lambda` represents the characteristic length scale of the system, :math:`E(\mathbf{r}_{1/2},N_{1/2})` are the configurational energies of the initial/final configurations, respectively and :math:`\beta = \frac{1}{kT}`.  The :math:`\frac{V\Lambda^{-3}}{N+1}` coefficient represents the fact that you can insert a particle anywhere in the system (inside a volume, *V*) but the likelihood of deleting that particle is :math:`\frac{1}{\mathrm{total number of particles} = \frac{1}{N + 1}`.  :math:`\Lambda` appears to conserve units and can be readily absorbed into the chemical potential.  Similarly, the acceptance criterion for particle deletions is given by:

.. math::

   P_{\mathrm{acc}}([\mathbf{r}_1,N_1] \rightarrow [\mathbf{r}_2,N_2] ) = \min(1,  \frac{N\Lambda^{3}}{V}\exp \{- \beta [E(\mathbf{r}_2,N_2) - E(\mathbf{r}_1,N_1)] \} )

where :math:`N = N_1` is the initial number of particles (before the deletion) and :math:`N - 1 = N_2` is the final number of particles (after the deletion).  For more information on the treatment of detailed balance in the Grand Canonical ensemble, see _[#f1].

In this session, we have defined our Lennard-Jones particles as 'molecules' made up of one atom.  For larger molecules, there are additional terms which come from the specific orientation of molecules.  Molecular rotations are difficult to model accurately in this way because the molecule can change its orientation between insertion and deletion moves, leading to technically 'different' molecules being inserted and deleted, breaking detailed balance.  |think| Does this apply to both linear and nonlinear molecules?
|think| What are the possible solutions to this problem? 

|think| Can molecular vibrations be modelled in Grand Canonical Monte Carlo simulations in a way that ensures detailed balance?

.. rubric:: Footnotes

.. [#f1] M. S. Shell, "Monte Carlo simulations in other ensembles"[online], University of California at Santa Barbara: Engineering, 2012.  Available from: https://engineering.ucsb.edu/~shell/che210d/Monte_Carlo_other_ensembles.pdf