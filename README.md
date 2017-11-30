# README - undergrad_MC_course#
(C) Chris King, James Grant, Steve Parker 

undergrad_MC_course is licensed under a Creative Commons Attribution 4.0 International License which you should have received with this work.  If not see <http://creativecommons.org/licenses/by/4.0/>.

AUTHORS
-------
Chris King, James Grant, Steve Parker
Department of Chemistry, University of Bath

Purpose
-------

This repository contains material to support the teaching of Monte Carlo methodologies in molecular simulation for Undergraduate courses e.g. Chemistry, Chemical Engineering, Physics.
The course makes use of the [CCP5](https://www.ccp5.ac.uk)  code [DL_MONTE](https://ccpforge.cse.rl.ac.uk/gf/project/dlmonte2/) for simulations which is availble for free to academic users.
The initial development of this material was supported by an Undergraduate Bursary from CCP5.
The material is a natural extension of the general tutorial material that was developed as part of the DL_MONTE2 development project, which is a valuable resource for accessing the advance methods in DL_MONTE2.
In contrast this material is aimed to be introductory and assumes little or no knowledge of Monte Carlo methods or molecular simulation.

The authors welcome feedback or suggestions raised as issues, or by emailing the authors: r.j.grant@bath.ac.uk.  If you use this material within your course in some or any form we would also be grateful if you could let us know.

Contents
--------

The course contains tutorial documentation and input files for DL_MONTE2 arranged as follows:
 
* <tutorials>
    * <jupyter>
        * <images>
        * <inputs>
    * <source>
        * <images>
    * <teaching_plans>
* <inputs>
    * <tut_{0..4>

The tutorials come are provided in rst which can be used to make html or latex->pdf versions of the worksheets.  We are also experimenting with developing a version of the material in the form of Jupyter notebooks to reduce the barrier to uptake of the material in line with developments in many courses at Bath and across academia.  We also aim to provide a teaching plan or learning outcomes for each of the tutorials which are currently available as Word docs.  DL_MONTE input files for each of the tutorials and in due course optional extenstion projects are also provided.  We are also developing slides to support the teaching of the course whihc will be included in due course.

Overview
--------

The use of Monte Carlo methods in molecular simulation is common in research, being applied to condensed and soft matter, in e.g. chemistry, physics.  
However, until recently a general purpose simulation code has been missing with research groups typically developing in house codes for specific methods/systems.
The CCP5 code DL_MONTE attempts to address this by providing Monte Carlo methods which make use of a format similar to and as compatible as possible with its sister code DL_POLY for Molecular Dynamics.

This material makes use of DL_MONTE to provide a framework for courses for use in Undergraduate degrees to introduce students to Monte Carlo methods.
The material is designed to be delivered as 4 or 5 half day guided practical workshops to demonstrate alternatives to more commonly taught Molecular Dynamics.
The tutorial are arranged as follows:
*   Session 0: an optional session that gives a brief introduction to Molecular Dynamics for those who have no previous experience with computational techniques.  This session illustrates the core concepts of Molecular Dynamics by applying them to the phase behaviour of a small system of Lennard-Jones particles.
*   Session 1: introduces the fundamental concepts of the Monte Carlo method for simulation and demonstrates this using a square 2D Ising model. 
*   Session 2: describes the operation and capabilities of DL_MONTE by using it to run simulations of the Lennard-Jones system from Session 0 under both NVT and NpT ensembles.  This also allows direct comparison between Molecular Dynamics and Monte Carlo techniques.  
*   Session 3: continues to use the Lennard-Jones system, but now explores simulations under the Grand Canonical ensemble.  This session illustrates the advantages of using ‘variable N’ ensembles in simulating certain physical phenomena over ‘fixed N’ ensembles.  
*   Session 4c: (Chemistry/Chem. Eng.) applies the Grand Canonical ensemble and Monte Carlo techniques on a more realistic chemical system; adsorption of CH4 onto a silaceous zeolite.  This session will help illustrate the additional considerations for dealing with real molecules in DL_MONTE and should provide context for the use of Monte Carlo techniques in academic research.
*   Session 4b: (Physics) provides an alternative to 4a by investigating the semi-Grand Ensemble of A and B Lennard-Jones particles.

We are also developing material for exercises that would allow the course to include extended mini-projects based on  sessions 3 and 4c, 4p.

### Who do I talk to? ###

James Grant, Research Software Engineer, University of Bath
Email: r.j.grant@bath.ac.uk
