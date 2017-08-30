.. _tutorial_2:

---------------------------------------------
TUTORIAL 2: 2D Ising Model using Monte Carlo
---------------------------------------------

Authors: Chris King, James Grant - r.j.grant@bath.ac.uk

Introduction to Monte Carlo Methods
===================================

In the previous session, you were introduced to classical potential modelling simulation techniques, specifically Molecular Dynamics.  This provides useful background knowledge for the entirely different simulation technique that we will be discussing in this session: Monte Carlo. 

Monte Carlo (MC) is the name given to the simulation technique that attempts to solve a problem by randomly sampling out of all of its possible outcomes ('configurational space')and obtaining a result based on numerical analysis of the sampling.  MC is a stochastic method, which means that the final state of the system cannot be predicted precisely based on the initial state and parameters, but through numerical analysis, reproducible results can be obtained.  This contrasts with other techniques like molecular dynamics, which are deterministic, where if you know the initial state and the inputs for the calculations, you can predict what the configuration of the system will be at any and all times thereafter.  This distinction allows MC to be used in a variety of applications across the scientific community where deterministic techniques are ineffective or impossible to use, such as phase co-existence and criticality, adsorption, and development of solid-state defects [#f1]_.

Results from MC simulations are generally accurate and reliable, assuming that the technique has representatively sampled the distribution of possible configurations in the system ('configurational space').  In other words, if our sampling method returns the probability distribution we expect, then we know that are sampling method is reliable.  In thermodynamic systems, the probability distribution of available states is given by the Boltzmann distribution:

.. math::

   W(\mathbf{r}) = \exp {\Bigl(\frac{E}{kT}\Bigr)} 

where :math:`W(\mathbf{r})` is the probability of being in a state of energy, also known as the statistical weight, *E*, at temperature, *T*, and *k* is the Boltzmann constant.  The ratio of Boltzmann distributions at two different energies, :math:`E_2` and :math:`E_1`, is known as the Boltzmann factor:

.. math::
  
   \frac{W(\mathbf{r}_1)}{W(\mathbf{r}_2)} = \exp {\Bigl(\frac{E_2 -E_1}{kT}\Bigr)}

So if our sampling method yields the Boltzmann distribution, we know that our simulation accurately reflects real systems.  There are many possible ways one can sample the configurational space of a simulated system, the intuitive case is simple random sampling in that we move randomly from one configuration to another.  However, this process is only reliable in systems with a constant probability distribution of states as it does not take into account the respective weighting of a given configuration.  For example, it can under-represent a small number of configurations who contribute significantly to the overall state of the system.

The concept of statistical weight is crucial in thermodynamics and describes how likely a particular configuration is of being observed out of a hypothetically *large* number of replicas of that system.  For instance, consider the possible configurations of the gas molecules in this room, clearly, this system would have a high probability of being in a configuration where the gas molecules are evenly (on average) distributed throughout the volume of the room and so this configuration has a high weighting.  Yet, there is a configuration where every gas molecule sits in one corner of the room, this configuration is highly unlikely to be seen and so its weighting is very low.  The weight of a particular configuration is given by:

.. math::

   W(\mathbf{r}) = \frac{\exp {\Bigl(\frac{- E(\mathbf{r})}{kT}\Bigr)}}{\sum_{i} \exp {\Bigl(\frac{- E(\mathbf{r_{i}})}{kT}\Bigr)} }

where :math:`E(\mathbf{r})` is the energy of a configuration :math:`\mathbf{r}`.  In MC simulations, the statistical weight of moving from a configuration, :math:`\mathbf{r_1}`, to a new configuration, :math:`\mathbf{r_2}`, is:

.. math::

   W(\mathbf{r}_1 \rightarrow \mathbf{r}_2) = \frac{W(\mathbf{r_1})P(\mathbf{r}_1 \rightarrow \mathbf{r}_2)}{N}

where :math:`W(\mathbf{r_1})` is the weight associated with :math:`\mathbf{r}_1`, :math:`P(\mathbf{r}_1 \rightarrow \mathbf{r}_2)` is the probability of moving from configuration :math:`\mathbf{r}_1` to :math:`\mathbf{r}_2` and *N* is the number of possible configurations. Figure 1 demonstrates the concept of statistical weights between moving from two configurations, A and B.  The corresponding weight of going from :math:`\mathbf{r}_2` back to :math:`\mathbf{r}_1` is:

.. math::

   W(\mathbf{r}_2 \rightarrow \mathbf{r}_1) = \frac{W(\mathbf{r_2})P(\mathbf{r}_2 \rightarrow \mathbf{r}_1)}{N}   

.. figure:: images/Tut_2_images/weights.png
   :align: center

   **Figure 1:** The associated statistical weights of moving between two configurations, A and B.

There are more sophisticated ways of sampling configurational space, such as the Metropolis Algorithm, which is one of the most widely used sampling schemes in MC simulations (including this one).  The function of the Metropolis algorithm is summarised in Figure 2.  First, it randomly selects a particle in the system and proposes a move to another configuration.  It then calculates the new energy of the configuration and compares it with the energy of the previous configuration before the move was proposed.  It then applies the following condition:

.. math::

         P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2) = \min(1, \exp \ \Bigl(- \frac{E(\mathbf{r}_2) - E(\mathbf{r}_1)}{kT}\Bigr) \ )

where :math:`P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2)` is the probability of accepting the move from the initial configuration, :math:`\mathbf{r}_1`, with an energy, :math:`E(\mathbf{r}_1)`, to the new configuration, :math:`\mathbf{r}_2`, with an energy, :math:`E(\mathbf{r}_2)`.  The function min() means that the smallest value in the brackets is chosen.  If the energy of the new configuration is less than that of the original, *i.e.* :math:`E(\mathbf{r}_2) < E(\mathbf{r}_1)`, then :math:`E(\mathbf{r}_2)-E(\mathbf{r}_1) < 0` and so :math:`\exp \ \Bigl(- \frac{E(\mathbf{r}_2) - E(\mathbf{r}_1)}{kT}\Bigr) \  > 1` and so the move is accepted with :math:`P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2) = 1`.  If the new energy is greater than the energy of the original configuration, *i.e.* :math:`E(\mathbf{r}_2) > E(\mathbf{r}_1)`, then :math:`E(\mathbf{r}_2)-E(\mathbf{r}_1) > 0` and so :math:`\exp \ \Bigl(- \frac{E(\mathbf{r}_2) - E(\mathbf{r}_1)}{kT}\Bigr) \  > 1` and the move is accepted with probability :math:`P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2) = \exp \ \Bigl(- \frac{E(\mathbf{r}_2) - E(\mathbf{r}_1)}{kT}\Bigr) \ < 1`.  

.. figure:: images/Tut_2_images/Metropolis_algorithm.png
   :align: center

   **Figure 2:** Visual representation of the function of the Metropolis algorithm.  Once one move outcome is complete, the algorithm repeats on the final configuration. 

|think| Even if the proposed move leads to a higher-energy configuration, there is still a non-zero probability of it being accepted! Why should this be the case?

|think| What happens to the total number of accepted moves in a given simulation as we change the temperature?

.. |think| image:: images/General/think.png
   :height: 100 px
   :scale: 25 %

This defines the concept of detailed balance:

.. math::

   W(\mathbf{r}_1 \rightarrow \mathbf{r}_2)P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2) = W(\mathbf{r}_2 \rightarrow \mathbf{r}_1)P_{\mathrm{acc}}(\mathbf{r}_2 \rightarrow \mathbf{r}_1)

We can now obtain the required Boltzmann distribution from this condition by rearrangement:

.. math::

   \frac{W(\mathbf{r}_2 \rightarrow \mathbf{r}_1)}{W(\mathbf{r}_1 \rightarrow \mathbf{r}_2)} = \frac{P_{\mathrm{acc}}(\mathbf{r}_1 \rightarrow \mathbf{r}_2)}{P_{\mathrm{acc}}(\mathbf{r}_2 \rightarrow \mathbf{r}_1)} = exp \ {\Bigl(\frac{E_2 -E_1}{kT}\Bigr)} 

This tells us that so long as we satisfy detailed balance, our system will be sampled according to the Boltzmann distribution and obey the rules of thermodynamics.  Though it is important to note that the condition of detailed balance is *sufficient* but *not necessary* to ensure that are system accurately reflects thermodynamics, *i.e.* there are simpler conditions one could employ that would ensure that our simulation obeys thermodynamics.  For instance, one could ensure that *balance* is achieved from the system which simply states that moving from one state to another state is the same for any initial and final state pairing, *i.e.*:

.. math::
   
   \frac{\mathrm{d}W(\mathbf{r}_1)}{\mathrm{d}t} = 0

However, detailed balance also ensures equilibrium between all states such that the trajectory from one configuration to another via several steps has the same probability as the reverse trajectory (See Figure 3).  This ensures the reliability of the sampling method used without requiring additional corrections in the calculations.

.. figure:: images/Tut_2_images/detailed_balance2.png
   :align: center

   **Figure 3:** A visualisation of the difference between the condition of balance (left) and detailed balance (right) for a set of different configurations, A-H, in the configurational space of a system.

Having discussed the concepts behind MC simulation methods, it is time to demonstrate how to apply them to a physical system.  This tutorial will be centred on a MC simulation of the magnetic properties of solid materials.

Ising Model of Magnetism
========================

An application where MC is more effective than deterministic methods is simulating the magnetic behaviour of solid state materials.  

Our simulation will be based on a 2D Ising model, which describes the macroscopic magnetic behaviour of a solid material as a result of the relative orientation of electron spins within the crystal lattice of a material.  As you may recall, each electron has an intrinsic 'spin'.  In simple terms, the spin of an electron can be thought of as a magnetic moment, with two possible orientations: 'up' and 'down'.  This idea helps define two classes of magnetic materials: diamagnetic and paramagnetic.  Diamagnetic materials are made up of atoms/molecules without unpaired electrons, do not interact with external magnetic fields, making them non-magnetic.  Paramagnetic materials contain unpaired electrons, exhibiting a net magnetic moment that can interact with external magnetic fields and give the material its magnetic properties.  Figure 4 shows an example of a paramagnetic material as a 2D lattice of colour-coded spins.

.. figure:: images/Tut_2_images/paramagnet_config.png
   :align: center

   **Figure 4:** A 2D schematic of a paramagnetic material under an external magnetic field.  Yellow indicates the spins that are aligned with the field and purple are spins that are anti-aligned.

There is another type of magnetism observed known as ferromagnetism, where instead of a uniform alignment of spins as in paramagnetic materials, 'domains' of aligned spins form, bound by domains of oppositely aligned spins (see Figure 5).  Ferromagnetic materials can show unique properties, such as being able to generate their own magnetic field (magnetisation) in the absence of an external magnetic field.  These form the common magnets seen in real-world applications.

.. figure:: images/Tut_2_images/ferromagnet_cand2.png
   :align: center

   **Figure 5:** A 2D schematic of a ferromagnetic material at :math:`T < T_{c}`.  Yellow and purple represent the two different spin orientations, 'up' and 'down', respectively.

The main factor influencing whether a given atom's spin is aligned with its neighbours in a crystal, and hence what type of magnetism the material displays, is its exchange energy, *E*, which in the Ising model is given by:

.. math::

	E = -J \sum_{<i,j>} s_{i}s_{j}

where *J* is the coupling constant between adjacent atoms in a given material and :math:`s_{i/j}` is the spin of the particle in position i/j in the lattice, respectively.  The <...> here mean the sum goes over the nearest neighbours of the atom in position (i,j), *i.e.* over the atoms at positions  (i-1, j), (i+1, j), (i, j-1) and (i, j+1) only.  The sign of *J* determines whether spin alignment (ferromagnetism) or anti-alignment (antiferromagnetism) is favourable.

The exchange energy can be thought of as an activation barrier for an atom to change its spin depending on the spins of its neighbours.  This means that, like with any physical system with an energy barrier, spontaneous thermal fluctuations can overcome the barrier and cause some atoms/domains to flip their spin, with the likelihood of flipping a spin increasing as temperature increases.  Therefore, ferromagnetic materials only show domains at temperatures under a specific critical, or Curie, temperature, :math:`T_{c}`.  Above this point, ferromagnetic materials lose their ability to retain magnetisation because the thermal fluctuations are much larger than the energy required to switch a domain's alignment with respect to other domains.  This results in a loss of the domain structure, and hence loss of magnetisation without an external field.  It is for this reason that paramagnetism can be thought of as high-temperature ferromagnetism.

For more information on the Ising model, consult either [#f2]_ or [#f3]_.

|think| The Metropolis algorithm is employed in these simulations, describe what constitutes a 'move' in the context of this system.

|think| Write an expression for the energy difference between the initial and final configurations, :math:`E(\mathbf{r}_2) - E(\mathbf{r}_1)`, for the 2D Ising model.

Exercise 1)
-----------

The aim of this exercise is to familiarise yourself with running MC calculations on a simple 2D Ising model of a ferromagnetic material. The material is represented by a 64x64 2D lattice of points, each representing an atom with its own net spin.  In this exercise, all atoms are spin-aligned.  We will be running a MC simulation to look at how the overall spin alignment (magnetisation) and energy of the system evolves with both time and temperature.

|action| Go to 'inputs' :math:`\rightarrow` 'Tut_2' :math:`\rightarrow` 'main' :math:`\rightarrow` 'Init' and copy the contents into a new directory in your domain.  The CONFIG file displays the initial configuration of your system, the CONTROL file allows you to set the parameters and constraints for your simulation, and the FIELD file describes all interactions in the system (although they may look slightly different to the ones presented in the last session, they perform the same roles).  Though we will be going through the function of these in detail in the next session, it may be helpful to have a look and familiarise yourself with their contents.  

.. |action| image:: images/General/action.png
   :scale: 5 %

**instructions for running a simulation**

As the calculation runs and completes, you will notice several new files appear in your directory.  These have similar roles to their counterparts from the previous session and will be explained in detail in the next tutorial.  The files you will be using throughout this tutorial will be the OUTPUT.000 and the PTFILE.000.  

Now that you have all the output data you could possibly need from this calculation, we shall proceed with extracting the following data from the OUTPUT.000 and PTFILE.000: the time evolution of magnetisation and the distribution of the magnetisations over the course of the simulation.  

|action| You will need to employ the 'analysis.sh' script by running the following command in the directory containing your output files::

	analysis.sh
	
The command should complete almost instantly and you should see several new files: M_seq.dat, M_hist.dat, M_hist.png, and M.dat.  These files contain: time-evolution of magnetisation, a normalised magnetisation frequency distribution (in both data and plotted forms), and the average magnetisation at the temperature of the simulation, respectively.

We shall now proceed to run the calculation at higher temperatures to obtain the temperature-dependence of the magnetisation.  

|action| Create a new directory for each temperature and copy the CONFIG, CONTROL and FIELD files from your first calculation to them.
  
|action| Open the CONTROL file in each and increase the temperature to a value of your choosing (HINT: you will not need to go above 5.0 K!) and run the calculations.    

|action| Once each calculation is complete, run the analysis script in the same manner as above to obtain the relevant data.

|action| From your calculations, plot magnetisation vs temperature for the system.  |think| Comment on the shape of your graph and estimate the critical temperature, :math:`T_{c}`, from it. *N.B.* it may be wise to run calculations at several temperatures around the perceived critical point.  

For any general 2D lattice where coupling along rows and along columns are equal, :math:`T_{c}` is given by:

.. math::

	T_{c} = \frac{2}{\ln(1+\sqrt{2})} \approx 2.269

|think| Does your estimation of :math:`T_{c}` agree with that predicted by the above equation? Account for any observed discrepancies.

|action| Plot the time-evolution of magnetisation (on the same graph) for:

	a) :math:`T < T_{c}`
	b) :math:`T \approx T_c`
	c) :math:`T > T_{c}`

|think| Comment on any differences between in these plots and rationalise them using your knowledge of ferromagnetism.  Do the results correspond to the predictions of the Ising model?

|action| Also, have a look at the magnetisation histogram for some of your temperatures and describe how the distribution of magnetisations appears to change with temperature.  |think| Does this behaviour support the rest of your output data?

Extension:
----------

You have seen what happens as the system is heated, but you can also look at the magnetisation upon cooling the system from a state above the critical temperature to a state below the critical temperature. 

|action| Take the REVCON from one of your simulations where :math:`T>T_{c}`, copy it into a new directory and rename it 'CONFIG'.  Also copy the CONTROL and FIELD files into this directory and change the temperature to :math:`\sim 10^{-3} K`.  Then run the simulation.  

|action| Once the simulation is complete, use the analysis.sh script to extract the output data and plot the time evolution of magnetisation.  Record your observations.  
|think| Does this agree with magnetic behaviour predicted by the Ising model? 

|think| How does this compare with the time evolution at :math:`T > T_{c}`?

Exercise 2)
-----------

This exercise will demonstrate the stochastic nature of MC simulations as well as how the Metropolis algorithm produces reliable and accurate results for this simple 2D Ising model.

We have seen what happens when we start the simulations from a fixed starting configuration (all spins aligned), but what will happen when we start from a random configuration? 

|action| Create a new directory and copy the CONFIG, CONTROL and FIELD files from one of your previous calculations into it. 

|action| Replace the line starting with 'seeds' to just 'ranseed'.  

|action| Make a note of the temperature and run the calculation and use analysis.sh on the output data as you have done in the previous exercise. 

|action| Run this calculation on these input files several times (WARNING: remember to copy the output files into separate directories each time before running the calculation again!) and plot the time-evolution of the magnetisation for each calculation.  Each of these calculations represent running the simulation on a different, randomly-generated initial configuration at the same temperature.  

|think| How does the final magnetisation of each random initial configuration compare with each other, *i.e.* does the initial configuration have an effect on the outcome of the simulation? 

Extension:
----------

|action| For one of your calculations, find out the initial configuration by typing the following into the command line::

	grep seeds OUTPUT.000

Running this command should return a line containing four integer numbers.  

|action| Create a new directory and copy the CONFIG, CONTROL and FIELD files into it.  

|action| Go to your CONTROL file and replace 'ranseed' with 'seeds int1 int2 int3 int4' where 'int' are the numbers from the command line.

|action| Re-run the calculation with this CONTROL file and plot the magnetisation vs time.  |think| Compare this with the equivalent \'ranseed\' calculation data.  

|think| What do you notice about the magnetisation evolution in the two calculations? Does this confirm that the stochastic nature of Monte Carlo methods can produce reliable results?

Conclusions:
============

Now that you have reached the end of this tutorial, you will hopefully have a better understanding of the Monte Carlo method and the motivation for its use. You have simulated the magnetic properties of a 2D material based on the Ising model and obtained:

- the temperature-dependence of magnetisation
- the evolution of magnetisation with time
- validation of the stochastic nature of Monte Carlo methods

In the next tutorial, you will be introduced to a general Monte Carlo program called DLMONTE and use it to model the thermal properties of a Lennard-Jones material.

Extensions (optional):
======================

1. Antiferromagnetism:
----------------------

So far, you have looked at how the magnetic behaviour of a ferromagnetic system changes over time and temperature, but there is another possible type of magnetism called antiferromagnetism, where the sign of the coupling constant, *J*, changes sign.  This means that it is now favourable for the spin of one atom to be opposed to the spin of its neighbours, resulting in a preferred 'checkerboard' pattern of magnetisation on the 2D lattice (see Figure 6).  You can investigate the magnetic behaviour in this case using the 2D Ising model.

.. figure:: images/Tut_2_images/antiferromagnet.png
   :align: center

   **Figure 6:** The most stable magnetic configuration of an antiferromagnetic material at :math:`T < T_{c}`.

**script for changing the signs of the coupling constants in the FIELD file and create a new directory for the initial calculation**

|action| Now investigate the magnetic properties of this material in a manner similar to what you have done in this tutorial.

|think| Compare your results of the antiferromagnet with the ferromagnet.  Rationalise any observed differences in terms of exchange energy and alignment of spins.

.. Link to next tutorial

.. rubric:: Footnotes

.. [#f1] S. Mordechai (Editor), *Applications of Monte Carlo Method in Science and Engineering* [Online]. Available: https://www.intechopen.com/books/applications-of-monte-carlo-method-in-science-and-engineering 
.. [#f2] J. V. Selinger, "Ising Model for Ferromagnetism" in *Introduction to the Theory of Soft Matter: From Ideal Gases to Liquid Crystals*.  Cham: Springer International Publishing, 2016, pp. 7-24.
.. [#f3] N. J. Giordano, *Computational Physics*.  Upper Saddle River, N.J.: Prentice Hall, 1997. 