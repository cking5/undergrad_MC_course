.. _tutorial_3:

-------------------------------------
TUTORIAL 3: Grand Canonical Ensemble
-------------------------------------

Authors: Chris King, James Grant - r.j.grant@bath.ac.uk

Introduction
============

So far in this course, we have been running simulations on systems with a constant number of particles (NVT or NPT ensembles).  This is useful for studying systems where there is no exchange of particles from one defined type to another and where predicting system behaviour by optimising particle number, *N*, is not required.  However, not all systems can be accurately modelled in this way, for instance, phase coexistence, where more than one phase of a material exist in equilibrium.  Other examples include: adsorption onto a surface (as you will see in this and future sessions), effect of solid-state defects on bulk properties and dissolution rates.  In this session, we will explore simulations where the number of particles is not constant and particles can be readily added or removed from our system.  

Ensembles are very important in thermodynamics and statistical mechanics, where one can readily observe the macroscopic properties of a system without being able to know or control its microscopic properties.  In simple terms, an ensemble can be thought of as a collection of a large number of replicas of a system, akin to running repeating an experiment many times with the same initial conditions.  By averaging the outcomes of each of these replicas, one can create a probability distribution of possible outcomes for an experiment in a given ensemble.  In this course and statistical mechanics you will come across several different ensembles which are described in the form of three letters, describing which parameters in the ensemble of the system are kept constant.  Ensembles are defined in this manner to ensure that they demonstrate statistical equilibrium, *i.e.* ensembles do not evolve over time, these include:

- *N*: total number of individual objects (atoms, molecules, particles etc.) in the system kept constant.
- :math:`\mu`: the chemical potential of the system kept constant
- *V*: system volume kept constant
- *P*: system pressure kept constant
- *T*: system temperature kept constant
- *E*: the total energy of the system kept constant

For example, 'NVT' means that the total number of particles, total volume and temperature are kept constant (known as the Canonical Ensemble) and 'NpT' means that the total number of particles, total pressure and temperature are kept constant (a.k.a. the Gibbs ensemble).  Figure 1 displays physical representations of the different ensembles that are used in thermodynamics.

.. figure:: images/Tut_3_images/Statistical_Ensembles.png
   :align: center
   
   **Figure 1**: Visual representation of different ensembles [#f1]_. 

However, we know that ensembles require three parameters to be kept constant, so what replaces *N*? The answer is chemical potential, :math:`\mu`. There are several definitions of chemical potential in the literature, but the most common definition is the effect of number of particles on the total Gibbs Free Energy, *G*, of the system:

.. math::

   \mu = \Bigl|\frac{\partial G}{\partial N}\Bigr|_{V, T} 

where :math:`\partial` means that you take the partial derivative of *G*.  A partial derivative is used for functions that vary with more than one variable, in this case, *G* can vary with: *N*, *V* and *T*, it means that you take the derivative of the function with respect to one of the variables while keeping the others fixed.  In this case, :math:`\mu` is the partial derivative of *G* with respect to *N* while keeping *V* and *T* constant.  For more information about the chemical potential, you can look in any general physical chemistry textbook or refer to online resources such as: (1)_, (2)_.

.. _(1): http://www.icsm.fr/Local/icsm/files/286/JFD_Chemical-potential.pdf

.. _(2): http://chem.atmos.colostate.edu/AT620/Sonia_uploads/ATS620_F11_Lecture5/Lecture5_AT620_083111.pdf

The :math:`\mu`\VT ensemble is known as the Grand Canonical ensemble and is the most widely-used 'variable *N*' ensemble.  Other simulation techniques, like Molecular Dynamics, are very inefficient at accurately modelling systems where the total number of particles isn't constant, which is why Monte Carlo is the preferred technique in these situations.
In fact, the limitations that we encountered in the previous sessions are readily overcome if we operate under the Grand Canonical ensemble.

Removing the constraint of having fixed particle number allows several new Monte Carlo move types to be defined, the most fundamental being insert/delete moves, where the proposed move adds or removes a particle from the system, respectively.  Other move types like swaps and replacements can be defined in terms of insert/delete moves.

The total energy of the system performing an insert/delete move, :math:`E(\mathbf{r}_2,N_2)`, is defined as:

.. math::

  E(\mathbf{r}_2,N_2) = U(\mathbf{r}_1,N_1) + \mu \Delta N

where :math:`U(\mathbf{r}_1)` is the energy of the initial configuration, :math:`\mathbf{r}_1`, of the system, :math:`N_1, N_2` are the initial and final numbers of particles, respectively, and :math:`\Delta N = N_2 - N_1`. |think| Is :math:`\Delta N` positive or negative for insertions? |think| What about deletions?

.. |think| image:: images/General/think.png
   :height: 75 px
   :scale: 25 %

Now that you are more familiar with DL_MONTE and the general Monte Carlo method, it is now a good time to discussed the concept of detailed balance, which ensures that simulations sample from the Boltzmann distribution for thermodynamic systems:

.. math::
  
   \frac{W(\mathbf{r}_1)}{W(\mathbf{r}_2)} = \exp {\Bigl(\frac{E_2 -E_1}{kT}\Bigr)}

where :math:`E_{1/2}` is the energy of configuration :math:`\mathbf{r}_{1/2}`, respectively, *k* is the Boltzmann constant, *T* is the temperature of the system, and :math:`W(\mathbf{r}_{1/2})` is the 'weight' of :math:`\mathbf{r}_{1/2}` of The concept of statistical weight is crucial in thermodynamics and describes how likely a particular configuration is of being observed out of a hypothetically *large* number of replicas of that system.  For instance, consider the possible configurations of the gas molecules in this room, clearly, this system would have a high probability of being in a configuration where the gas molecules are evenly (on average) distributed throughout the volume of the room and so this configuration has a high weighting.  Yet, there is a configuration where every gas molecule sits in one corner of the room, this configuration is highly unlikely to be seen and so its weighting is very low.  The weight of a particular configuration is given by:

.. math::

   W(\mathbf{r}) = \frac{\exp {\Bigl(\frac{- E(\mathbf{r})}{kT}\Bigr)}}{\sum_{i} \exp {\Bigl(\frac{- E(\mathbf{r_{i}})}{kT}\Bigr)} }

where :math:`E(\mathbf{r})` is the energy of a configuration :math:`\mathbf{r}`.  In MC simulations, the statistical weight of moving from a configuration, :math:`\mathbf{r_1}`, to a new configuration, :math:`\mathbf{r_2}`, is:

.. math::

   W(\mathbf{r}_1 \rightarrow \mathbf{r}_2) = \frac{W(\mathbf{r_1})P(\mathbf{r}_1 \rightarrow \mathbf{r}_2)}{N}

where :math:`W(\mathbf{r_1})` is the weight associated with :math:`\mathbf{r}_1`, :math:`P(\mathbf{r}_1 \rightarrow \mathbf{r}_2)` is the probability of moving from configuration :math:`\mathbf{r}_1` to :math:`\mathbf{r}_2` and *N* is the number of possible configurations. Figure 1 demonstrates the concept of statistical weights between moving from two configurations, A and B.  The corresponding weight of going from :math:`\mathbf{r}_2` back to :math:`\mathbf{r}_1` is:

.. math::

   W(\mathbf{r}_2 \rightarrow \mathbf{r}_1) = \frac{W(\mathbf{r_2})P(\mathbf{r}_2 \rightarrow \mathbf{r}_1)}{N}   

.. figure:: images/Tut_3_images/weights.png
   :align: center

   **Figure 1:** The associated statistical weights of moving between two configurations, A and B.

As you may recall, we use the Metropolis algorithm in this course to accept/reject proposed moves according to the following condition:

.. math::

         P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2) = \min(1, \exp \ \Bigl(- \frac{E(\mathbf{r}_2) - E(\mathbf{r}_1)}{kT}\Bigr) \ )

The statistical weight of a configuration amongst a given distribution of configurations and the acceptance probability for a move define the condition of detailed balance:

.. math::

   W(\mathbf{r}_1 \rightarrow \mathbf{r}_2)P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2) = W(\mathbf{r}_2 \rightarrow \mathbf{r}_1)P_{\mathrm{acc}}(\mathbf{r}_2 \rightarrow \mathbf{r}_1)

We can now obtain the required Boltzmann distribution from this condition by rearrangement:

.. math::

   \frac{W(\mathbf{r}_2 \rightarrow \mathbf{r}_1)}{W(\mathbf{r}_1 \rightarrow \mathbf{r}_2)} = \frac{P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2)}{P_{\mathrm{acc}}(\mathbf{r}_2 \rightarrow \mathbf{r}_1)} = exp \ {\Bigl(\frac{E_2 -E_1}{kT}\Bigr)} 

This tells us that so long as we satisfy detailed balance, our system will be sampled according to the Boltzmann distribution and obey the rules of thermodynamics.  Though it is important to note that the condition of detailed balance is *sufficient* but *not necessary* to ensure that are system accurately reflects thermodynamics, *i.e.* there are simpler conditions one could employ that would ensure that our simulation obeys thermodynamics.  For instance, one could ensure that *balance* is achieved from the system which simply states that moving from one state to another state is the same for any initial and final state pairing, *i.e.*:

.. math::
   
   \frac{\mathrm{d}W(\mathbf{r}_1)}{\mathrm{d}t} = 0

However, detailed balance also ensures equilibrium between all states such that the trajectory from one configuration to another via several steps has the same probability as the reverse trajectory (See Figure 3).  This ensures the reliability of the sampling method used without requiring additional corrections in the calculations.

.. figure:: images/Tut_3_images/detailed_balance3.png
   :align: center

   **Figure 3:** A visualisation of detailed balance (right) for a set of different configurations, A-H, in the configurational space of a system.

As in the previous sessions, we will be using DL_MONTE run Monte Carlo calculations on the phase behaviour of our all-too-familiar Lennard-Jones material.  However, all of our calculations in this tutorial will be conducted under the Grand Canonical ensemble.  You will hopefully see that we can get a more accurate reflection of the phase behaviour of real systems than if we are restricted to either NVT or NpT ensembles.

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
   MOLECULE lj 1 1                                              # continues to define the remaining molecules in the system
   LJ   core
   0.0000000000000000 -5.0000000000000000 -5.0000000000000000 
   etc

This takes the same basic structure as the CONFIG files from the previous session.  There are a few minor differences, for instance the '1' in line 2 leads to the slightly different way of presenting the coordinates of each particle.  The dimensions in lines 3-5 describe the system as a cube with sides of length 10 Angstroms.  The *NUMMOL* line tells us that there can be any number of molecules in the system between 8 and 1000.  The rest of the file defines 8 molecules called 'lj` that contain 1 atom called 'LJ' with a maximum of 1 atom in them.  You will notice that the CONFIG is much smaller than its counterpart used in the last session.  This is because the number of particles (or molecules in this case) will vary over the course of the simulation, we need to only specify the initial configuration, which will start with only 8 molecules.  In principle, you can define the locations of any number of molecules in the CONFIG file (as long as that number falls between the minimum and maximum numbers stated in the 'NUMMOL' line), but for the purposes of this tutorial, we start at the minimum number: 8.

CONTROL
-------

The CONTROL file will take the following form in this tutorial:

.. code-block:: html
   :linenos:
  
   GCMC Lennard-Jones              
   finish                          
   seeds 12 34 56 78               
   temperature     1.4283461511745 
   # nbrlist auto                  
   # maxnonbondnbrs 512            
   steps          10000            
   equilibration    0              
   print           1000            
   stack           1000            
   sample coord   10000            
   revconformat DL_MONTE            
   archiveformat dlpoly4            
                                   
   yamldata 1000                   
   move gcinsertmol 1 100 0.7      # Perform insertion/removal moves for lj, a weight 100 with a min. distance of 0.7 from atoms
   lj  0.06177                     # Use an activity of 0.06177   
   #  move atom 1 512              
   #  LJ core 
   #  move volume cubic linear 1   
   start                           

The lines that switch on the neighbour lists: 'nbrlist' and 'maxnonbondnbrs' have been suspended in this session.  This is because the no benefit in maintaining the list under :math:`\mu`\VT ensembles.  We have also suspended atom translation moves for simplicity (though there is nothing in principle wrong with allowing these types of moves), and volume moves since we work under a constant-volume ensemble.  There are two new lines present: the first describes the insert/delete moves for these simulations, with the first number stating how many molecules are inserted/deleted, the second being the weight of the proposed moves and the third being the minimum insertion distance from any other molecules present in the system.

In this calculation DL_MONTE is using the activity *a* rather than the chemical potential :math:`\mu`, which are related according to: 

.. math::

  a = \exp \Bigl(\frac{\mu}{RT}\Bigr)

where *R* is the gas constant.  This means that small changes to :math:`\mu` can have a large impact on the activity (assuming that *T* is constant).  In your inputs folder you will notice a file called 'activity-chempotential.txt', which lists values of :math:`\mu` and *a* in the first and second columns, respectively at *T* = 1.43 K:

|action| Plot the data in 'activity-chempotential.txt'.

.. |action| image:: images/General/action.png
   :scale: 5 %

|think| By using your graph or otherwise, estimate the value of :math:`\mu` from the value of *a* given in the CONTROL file.

FIELD
-----

The FIELD file looks almost identical to the ones from the previous session:

.. code-block:: html
   :linenos:

   Lennard-Jones                  
   CUTOFF 2.5                     
   UNITS internal                 
   NCONFIGS 1                     
   ATOMS 1                        
   LJ core 1.0  0.0               
   MOLTYPES 1                     
   lj                             
   ATOMS 1 1                      
   LJ core 0.0 0.0 0.0            
   FINISH                         
   VDW 1                          
   LJ core  LJ core lj   1.0 1.0  
   CLOSE                          

In the NVT and NpT cases all the particles were declared to be part of the same molecule, now each particle is a molecule in its own right.  This distinction is made to simplify the calculation under :math:`\mu`\VT ensembles.  In principle, atoms can be added or removed from a molecule however, for simplicity, we shall insert or delete whole molecules rather than parts of molecules.  Since we have a single Lennard-Jones particle in each molecule we simply position the particle at the 'origin' of the molecule.

Remember, there must be correspondence between the CONFIG and FIELD files, *i.e.* the number of molecule and atom types should be the same in both files.  Also remember that the number of interactions stated in the 'VDW' line must correspond to the number of interactions defined between it and the 'CLOSE' statement.

Exercise 1)
===========

As in the previous session, we need to ensure that the system has reached its equilibrium state before output data is obtained.  As you may recall, the amount of time that the system needs to equilibrate is stated by the 'equilibration' line in the CONTROL file and is different for every simulated system.  It is standard procedure for the user (*i.e.* you!) to determine what value the equilibration is for their system before obtaining results, so this is what we shall do now.

|action| Navigate to 'inputs' :math:`\rightarrow` 'Tut_4' :math:`\rightarrow` 'main' :math:`\rightarrow` 'Equil'.  You will see the standard DL_MONTE inputs files: CONFIG, CONTROL and FIELD, as well as some scripts for use later.  

|action| Run the DL_MONTE calculations as you have done in the previous session (quick reminder of how to do it). Extract the time-sequence of the number of particles in the system by using the following script::

  [user@node-sw-119 tut_4]  strip_gcmc.sh

When using 'fixed *N*' ensembles, like NVT and NpT, the simplest way to infer the equilibration of a system is to plot the system energy over the course of the simulation and define the equilibration as the number of steps in the simulation needed for the energy to fluctuate around some constant value.  Under the GC ensembles, this does not apply, instead we plot *N* over the course of the simulation and find the number of steps required for *N* to become roughly constant.

|action| By plotting the time-evolution of *N* for each of your simulations, increase the number of steps to determine when the system reaches equilibrium.

*N.B.* You will see that the output files will be mostly unchanged, except the YAMLDATA, which displays the number of molecules present instead of energies.

Exercise 2)
===========

Now that you know how to estimate the equilibration time needed for systems under the GC ensemble, we will now vary both temperature and activity and determine how these parameters affect *N*.

|action| Open the 'GCMC' folder in the 'main' folder.  |action| Replace the number of steps in the CONTROL file  with the value that you obtained from exercise 1.  

|action| Run simulations at various different temperatures and activities by varying the appropriate values in the CONTROL file.  

|action| Ensure that the system has equilibrated for each of your calculations. 

|action| Plot the time-evolution of *N* for each of your simulations.  

|think| What happens to the total number of particles over the course of the simulation as you vary the temperature and activity? 

|think| From your results and your own knowledge, how does the value of :math:`\mu` change the ease at which particles are:

 a) inserted
 b) deleted 

|action| You can also create histograms of the number of particles in the system over the course of the simulation, once you have produced the time sequence, with the script::

  hist.sh nmol.dat j

where *j* is the width of each bin used to generate the histogram.  You must specify the value of *j* in the command.  Though you are free to vary *j*, it is recommended that you set :math:`j = 1`.  Feel free to explore the effect *j* has on the shape of your histogram.

|think| How does the shape of the histogram vary with temperature?

|think| One could, in principle, also choose to use an :math:`\mu`\PT ensemble, what kind of problems could arise when running simulations under this ensemble?

|think| Define a swap move, where a particle at one position is swapped with another particle at a different position, as a sequence of insert and delete moves.

|think| Define a replacement move, where a particle in one position is changed to a particle of a different type, but remains in the same position, as a sequence of insert and delete moves.

Conclusions:
============

In this session, you have been introduced to the Grand Canonical (GC) ensemble, where the total number of particles in the system can vary but the chemical potential of the system remains constant.  You have demonstrated the use of the GC ensemble by investigating the thermal behaviour of a simple Lennard-Jones system and appreciated the advantages of using GC over 'fixed *N*' ensembles.  In the next session, we will apply the GC ensemble to the physical system of methane adsorption onto the surface of a zeolite in order to predict the conditions for ideal adsorption.

Extensions (optional):
======================

1. Detailed balance in the Grand Canonical ensemble
---------------------------------------------------

Like with the inclusion of volume moves in the previous session, the conditions through which detailed balance is maintained when employing insert/delete moves in :math:`\mu`\VT ensemble must be altered, such that, for particle insertions, the acceptance probability in the Metropolis algorithm in moving from an initial configuration, :math:`\mathbf{r}_1`, with :math:`N_1 = N` particles, to a final configuration, :math:`\mathbf{r}_2`, with :math:`N_2 = N + 1` particles is:

.. math::
  
   P_{\mathrm{acc}}([\mathbf{r}_1,N_1] \rightarrow [\mathbf{r}_2,N_2] ) = \min(1,  \frac{V\Lambda^{-3}}{N+1} \exp \{- \beta [E(\mathbf{r}_2,N_2) - E(\mathbf{r}_1,N_1)] \} )

where :math:`V` is the system volume, :math:`\Lambda` represents the characteristic length scale of the system, :math:`E(\mathbf{r}_{1/2},N_{1/2})` are the configurational energies of the initial/final configurations, respectively and :math:`\beta = \frac{1}{kT}`.  The :math:`\frac{V\Lambda^{-3}}{N+1}` coefficient represents the fact that you can insert a particle anywhere in the system (inside a volume, *V*) but the likelihood of deleting that particle is :math:`\frac{1}{\mathrm{N_{tot}}} = \frac{1}{N + 1}`.  :math:`\Lambda` appears to conserve units and can be readily absorbed into the chemical potential.  Similarly, the acceptance criterion for particle deletions is given by:

.. math::

   P_{\mathrm{acc}}([\mathbf{r}_1,N_1] \rightarrow [\mathbf{r}_2,N_2] ) = \min(1,  \frac{N\Lambda^{3}}{V}\exp \{- \beta [E(\mathbf{r}_2,N_2) - E(\mathbf{r}_1,N_1)] \} )

where :math:`N = N_1` is the initial number of particles (before the deletion) and :math:`N - 1 = N_2` is the final number of particles (after the deletion).  For more information on the treatment of detailed balance in the Grand Canonical ensemble, see [#f1]_.

In this session, we have defined our Lennard-Jones particles as 'molecules' made up of one atom.  For larger molecules, there are additional terms which come from the specific orientation of molecules.  Molecular rotations are difficult to model accurately in this way because the molecule can change its orientation between insertion and deletion moves, leading to technically 'different' molecules being inserted and deleted, breaking detailed balance.  |think| Does this apply to both linear and nonlinear molecules?
|think| What are the possible solutions to this problem? 

|think| Can molecular vibrations be modelled in Grand Canonical Monte Carlo simulations in a way that ensures detailed balance?

.. rubric:: Footnotes

.. [#f1] M. S. Shell, "Monte Carlo simulations in other ensembles"[online], University of California at Santa Barbara: Engineering, 2012.  Available from: https://engineering.ucsb.edu/~shell/che210d/Monte_Carlo_other_ensembles.pdf