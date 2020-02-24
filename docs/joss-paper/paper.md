---
title: 'Chaste: Cancer, Heart and Soft Tissue Environment'
tags:
  - C++
  - cardiac electrophysiology
  - monodomain
  - bidomain
  - simulation
  - cancer
  - heart
  - lung
  - cell-based
authors:
 - name: Fergus R Cooper
   orcid: 0000-0002-6221-4288
   affiliation: 1
 - name: Ruth E Baker
   orcid: 0000-0002-6304-9333
   affiliation: 1
 - name: Miguel O Bernabeu
   orcid: 0000-0002-6456-3756
   affiliation: 2
 - name: Rafel Bordas
   orcid: 0000-0002-0856-544X
   affiliation: 3
 - name: Louise Bowler
   orcid: 0000-0002-4910-9205
   affiliation: 3
 - name: Alfonso Bueno-Orovio
   orcid: 0000-0002-1634-3601
   affiliation: 3
 - name: Helen M Byrne
   orcid: 0000-0003-1771-5910
   affiliation: 1
 - name: Valentina Carapella
   orcid: 0000-0001-9750-9425
   affiliation: 4
 - name: Louie Cardone-Noott
   affiliation: 3
 - name: Jonathan Cooper
   orcid: 0000-0001-6009-3542
   affiliation: 5
 - name: Sara Dutta
   affiliation: 3
 - name: Benjamin D Evans
   orcid: 0000-0002-1734-6070
   affiliation: "6, 7"
 - name: Alexander G Fletcher
   orcid: 0000-0003-0525-4336
   affiliation: "8, 9"
 - name: James A Grogan
   affiliation: 1
 - name: Wenxian Guo
   affiliation: 10
   orcid: 0000-0001-8130-3326
 - name: Daniel G Harvey
   affiliation: 3
 - name: Maurice Hendrix
   orcid: 0000-0002-6621-7996
   affiliation: "11,12"
 - name: David Kay
   orcid: 0000-0003-4984-3669
   affiliation: 3
 - name: Jochen Kursawe
   orcid: 0000-0002-0314-9623
   affiliation: 13
 - name: Philip K Maini
   orcid: 0000-0002-0146-9164
   affiliation: 1
 - name: Beth McMillan
   orcid: 0000-0002-2605-1733
   affiliation: 3
 - name: Gary R Mirams
   orcid: 0000-0002-4569-4312
   affiliation: 11
 - name: James M Osborne
   orcid: 0000-0002-5622-0104
   affiliation: 14
 - name: Pras Pathmanathan
   orcid: 0000-0003-2111-6689
   affiliation: 15
 - name: Joe M Pitt-Francis
   orcid: 0000-0002-5094-5403
   affiliation: 3
 - name: Martin Robinson
   orcid: 0000-0002-1572-6782
   affiliation: 3
 - name: Blanca Rodriguez
   orcid: 0000-0002-3435-7388
   affiliation: 3
 - name: Raymond J Spiteri
   affiliation: 8
 - name: David J Gavaghan
   orcid: 0000-0001-8311-3200
   affiliation: 3
affiliations:
 - name: Wolfson Centre for Mathematical Biology, Mathematical Institute, University of Oxford, Oxford, UK
   index: 1
 - name: Centre for Medical Informatics, Usher Institute, University of Edinburgh, Edinburgh, United Kingdom
   index: 2
 - name: Department of Computer Science, University of Oxford, Oxford, UK
   index: 3
 - name: Division of Cardiovascular Medicine, Radcliffe Department of Medicine, University of Oxford, Oxford, UK
   index: 4
 - name: Research IT Services, University College London, London, UK
   index: 5
 - name: Centre for Biomedical Modelling and Analysis, Living Systems Institute, University of Exeter, Exeter, UK
   index: 6
 - name: School of Psychological Science, University of Bristol, Bristol, UK
   index: 7
 - name: School of Mathematics & Statistics, University of Sheffield, Sheffield, UK
   index: 8
 - name: Bateson Centre, University of Sheffield, Sheffield, UK
   index: 9
 - name: Department of Computer Science, University of Saskatchewan, Canada
   index: 10
 - name: Centre for Mathematical Medicine & Biology, School of Mathematical Sciences, University of Nottingham, Nottingham, UK
   index: 11
 - name: Digital Research Service, University of Nottingham, Nottingham, UK
   index: 12
 - name: Faculty of Biology, Medicine and Health, University of Manchester, Manchester, UK
   index: 13
 - name: School of Mathematics and Statistics, University of Melbourne, Victoria, Australia
   index: 14
 - name: Office of Science and Engineering Laboratories (OSEL), Center for Devices and Radiological Health (CDRH), U.S. Food and Drug Administration (FDA), Silver Spring, MD 20993, USA
   index: 15

date: 12th August 2019
bibliography: paper.bib
---

# Summary

Chaste (Cancer, Heart And Soft Tissue Environment) is an open source simulation package for the numerical solution of mathematical models arising in physiology and biology.

To date, Chaste development has been driven primarily by applications that include continuum modelling of cardiac electrophysiology ('Cardiac Chaste'), discrete cell-based modelling of soft tissues ('Cell-based Chaste'), and modelling of ventilation in lungs ('Lung Chaste').
Cardiac Chaste addresses the need for a high-performance, generic, and verified simulation framework for cardiac electrophysiology that is freely available to the scientific community.
Cardiac chaste provides a software package capable of realistic heart simulations that is efficient, rigorously tested, and runs on HPC platforms.
Cell-based Chaste addresses the need for efficient and verified implementations of cell-based modelling frameworks, providing a set of extensible tools for simulating biological tissues.
Computational modelling, along with live imaging techniques, plays an important role in understanding the processes of tissue growth and repair.
A wide range of cell-based modelling frameworks have been developed that have each been successfully applied in a range of biological applications.
Cell-based Chaste includes implementations of the cellular automaton model, the cellular Potts model, cell-centre models with cell representations as overlapping spheres or Voronoi tessellations, and the vertex model.
Lung Chaste addresses the need for a novel, generic and efficient lung modelling software package that is both tested and verified.
It aims to couple biophysically-detailed models of airway mechanics with organ-scale ventilation models in a package that is freely available to the scientific community.

Chaste is designed to be modular and extensible, providing libraries for common scientific computing infrastructure such as linear algebra operations, finite element meshes, and ordinary and partial differential equation solvers.
This infrastructure is used by libraries for specific applications, such as continuum mechanics, cardiac models, and cell-based models.
The software engineering techniques used to develop Chaste are intended to ensure code quality, re-usability and reliability.
Primary applications of the software include cardiac and respiratory physiology, cancer and developmental biology.

## The software

Chaste is available on GitHub [https://github.com/Chaste/Chaste](https://github.com/Chaste/Chaste), and the current stable release is version 2019.1.
Please see the Readme.md file on the Github repository for links to the Chaste wiki and install guides.

Previous publications about Chaste have detailed the rationale for, and design principles behind, Chaste [@Pitt-Francis2009Chaste], as well as the main application areas of Chaste up to 2013 [@Mirams2013Chaste].

Chaste places an emphasis on reproducibility and verification and, as such, extensive automated testing is used to ensure software quality and reliability.
A series of test suites must all pass before any commit is considered a release-candidate.
Most testing is performed on Long Term Support (LTS) versions of Ubuntu Linux, with unit tests additionally being run on macOS.

Testing includes compilation of all libraries with GCC, Clang and Intel C++ compilers; extensive unit testing; performance profiling to identify any slowdowns over time; memory testing with valgrind; verification of code coverage; and running unit tests with different combinations of dependencies to ensure portability.
The output of these tests is available at [https://chaste.cs.ox.ac.uk/buildbot/](https://chaste.cs.ox.ac.uk/buildbot/).

Since 2013, Chaste has substantially changed to modernise its infrastructure and to enable new science.
In terms of infrastructure, Chaste now uses a modern CMake build system, the C++14 language standard, and makes extensive use of BuildBot for continuous integration.
In terms of science, Lung Chaste is entirely new and allows the use of Chaste in a new scientific domain.
In Cardiac Chaste, we can now create algebraic Jacobians for CellML ODE systems, which can improve speed of simulation for cardiac action potential and tissue simulations [@Cooper2015Cellular], and metadata annotations of CellML files have replaced manual specification of variables in configuration files.
Cell-based Chaste has been overhauled to improve flexibility.
Changes include hierarchies of simulation modifiers, information writers, cell-cycle models, subcellular reaction network models, and numerical methods that allow new customisation points in almost every area of all cell-based simulations.
In addition, simulation output has been standardised to use VTK, a standard and powerful visualisation framework, and some cell-centre simulations now run in parallel using MPI.

### Comparison with other software

Chaste provides substantial common infrastructure enabling a wide range of applications across multiple disciplines.
Common elements include meshing, solving differential equations, input/output and continuum mechanics, and these form a platform for Cardiac, Cell-based and Lung Chaste.

A key goal of Chaste is to enable the implementation of many different modelling frameworks.
This not only allows a user to select the most appropriate tool for their research but, importantly, enables the comparison of different modelling frameworks to better understand the benefits and drawbacks of each [@Osborne2017Comparing].
This is an explicit design goal of Chaste, which focusses on the flexibility of implementing multiple models rather than (for example) building a graphical user interface.
See Table 1 for a comparison of alternatives to Chaste in specific domains, with all other software tools implementing a single modelling framework.

| Software    | Open Source |  GUI  |  CA   |  CP   |  PM   |  VT   |  VM   |
| ----------- | :---------: | :---: | :---: | :---: | :---: | :---: | :---: |
| Chaste      |      x      |       |   x   |   x   |   x   |   x   |   x   |
| CompuCell3D |      x      |   x   |       |   x   |       |       |       |
| Morpheus    |      x      |   x   |       |   x   |       |       |       |
| EPISIM      |             |   x   |       |       |   x   |       |       |
| CellSys     |             |   x   |       |       |   x   |       |       |
| PhysiCell   |      x      |       |       |       |   x   |       |       |
| Biocellion  |             |       |       |       |   x   |       |       |
| VirtualLeaf |      x      |   x   |       |       |       |       |   x   |
| EmbryoMaker |      x      |   x   |       |       |   x   |       |       |

**Table 1:**
A comparison of software tools for cell-based modelling.
GUI: graphical user interface.
CA: cellular automata.
CP: cellular Potts.
PM: particle model, a cell-centre model.
VT: Voronoi tessellation, a cell-centre model.
VM: vertex model.
References: CompuCell3D [@Swat2012CompuCell3D], Morpheus [@Starruss2014Morpheus], EPISIM [@Sutterlin2013EPISIM], CellSys [@Hoehme2010CellSys], PhysiCell [@Ghaffarizadeh2018PhysiCell], Biocellion [@Kang2014Biocellion], VirtualLeaf [@Merks2011VirtualLeaf], EmbryoMaker [@Marin-Riera2015EmbryoMaker].

### Installation

Installation of Chaste has been greatly simplified through the development of a Docker image [https://github.com/chaste/chaste-docker](https://github.com/chaste/chaste-docker). Docker is a lightweight, open-source virtualisation technology for running encapsulated applications ('containers') on all major operating systems at near-native speed. This enables Chaste (including all dependencies, environment settings, convenience scripts and the latest precompiled release) to be downloaded and installed with just a single command. Isolating Chaste within a container also means that its dependencies and those installed on the user's host system can coexist without interference or version conflicts.

In addition to simplifying the set-up and execution of Chaste, importantly this also enhances its reproducibility by providing a homogeneous computational environment regardless of the underlying operating system and hardware. Not only is the Chaste source code version-controlled, but so too are the dependencies, configuration settings and environment variables used to build and run it. This means that collaborators and reviewers can easily and consistently reproduce results (to within machine precision) on any platform while developers can seamlessly migrate and scale-up their simulations from a laptop to a workstation or HPC cluster.

## Example usage

Chaste has tutorials to walk users through basic functionality for each application area.
Tutorial examples are bundled for each specific release version, and examples for this release are available at [https://chaste.cs.ox.ac.uk/chaste/tutorials/release_2019.1](https://chaste.cs.ox.ac.uk/chaste/tutorials/release_2019.1).

Tutorials take the form of C++ header files that each define 'tests' in the Chaste testing infrastructure.
These tests must be compiled and run to produce an output, which can be visualised using [ParaView](https://www.paraview.org/).

In the following sections we showcase a specific tutorial for each of cardiac, cell-based, and lung Chaste, with minimal commands necessary to reproduce the output shown.

### Cardiac example

Here we demonstrate how to run and visualise a three-dimensional monodomain cardiac simulation.
This follows the tutorial `TestMonodomain3dRabbitHeartTutorial` which simulates the result of an electrical stimulus being applied to a realistic rabbit heart geometry.
Assuming that

 1. Chaste has been installed on Ubuntu Linux (or is running within a Docker container),
 1. the Chaste source code exists at `$CHASTE_SOURCE_DIR`,
 1. the environment variable `$CHASTE_TEST_OUTPUT` is set to a valid directory,

a minimal set of commands to build and run the tutorial is as follows:

```shell
mkdir build && cd build
cmake $CHASTE_SOURCE_DIR
make TestMonodomain3dRabbitHeartTutorial
ctest -R TestMonodomain3dRabbitHeartTutorial
```

This will produce output in the following directory:

```shell
$CHASTE_TEST_OUTPUT/Monodomain3dRabbitHeart
```

To view the results evolving over time as an animation in ParaView it is necessary to post-process the results with the following command:

```shell
cd $CHASTE_TEST_OUTPUT/Monodomain3dRabbitHeart/vtk_output
python $CHASTE_SOURCE_DIR/python/utils/AddVtuTimeAnnotations.py \
       results.vtu annotated_results.vtu
```

To visualise the output, open the file `annotated_results.vtu` in ParaView, and select to colour by `V` (voltage).

![rabbit heart](figs/RabbitTransparent.png)\

**Figure 1:**
Trans-membrane voltage on the rabbit heart mesh at the end of the simulation. As viewed on the surface from the apex of the heart (left) and on a wireframe showing the ventricular cavities (right).

### Cell-based example

Here we demonstrate how to run and visualise a cell sorting simulation using Chaste's vertex model implementation.
This follows the tutorial `TestRunningDifferentialAdhesionSimulationsTutorial`.
Assuming that

 1. Chaste has been installed on Ubuntu Linux (or is running within a Docker container),
 1. the Chaste source code exists at `$CHASTE_SOURCE_DIR`,
 1. the environment variable `$CHASTE_TEST_OUTPUT` is set to a valid directory,

a minimal set of commands to build and run the tutorial is as follows:

```shell
mkdir build && cd build
cmake $CHASTE_SOURCE_DIR
make TestRunningDifferentialAdhesionSimulationsTutorial
ctest -R TestRunningDifferentialAdhesionSimulationsTutorial
```

This will produce output in the following directory:

```shell
$CHASTE_TEST_OUTPUT/TestVertexBasedDifferentialAdhesionSimulation
```

To visualise the simulation, open the file `results.pvd` in ParaView, choose to colour by 'Cell types', and display 'Surface With Edges'.

![cell sorting simulation](figs/CellSorting.png)\

**Figure 2:**
The initial configuration of cells (left), and the final configuration of cells after sorting has occurred (right).

### Lung example

Here we demonstrate how to run and visualise the lung airway generation tutorial.
This follows the tutorial `TestAirwayGenerationTutorial` which statistically generates lung airways given initial geometry segmented from a CT scan.
Assuming that

 1. Chaste has been installed on Ubuntu Linux (or is running within a Docker container),
 1. the Chaste source code exists at `$CHASTE_SOURCE_DIR`,
 1. the environment variable `$CHASTE_TEST_OUTPUT` is set to a valid directory,

a minimal set of commands to build and run the tutorial is as follows:

```shell
mkdir build && cd build
cmake $CHASTE_SOURCE_DIR
make TestAirwayGenerationTutorial
ctest -R TestAirwayGenerationTutorial
```

This will produce output in the following directory:

```shell
$CHASTE_TEST_OUTPUT/TestAirwayGenerationTutorial
```

To visualise the generated airway geometry, open the file `example_complete_conducting_airway.vtu` in ParaView.
Application of an 'Extract Surface' filter followed by a 'Tube' filter allows the centreline and radius information to be viewed as a series of tubes.

![airway generation](figs/AirwayGen.png)\

**Figure 3:**
The initial geometry of major airways segmented from a CT scan (left), and an example of a complete generated airway tree (right).

## Recent publications enabled by Chaste

Since our last publication on Chaste [@Mirams2013Chaste], over 70 peer-reviewed publications have been enabled, which we mention briefly below.

Publications using Cardiac Chaste have included scientific studies relating to:
basic mechanisms of cardiac electrophysiology in healthy and diseased settings [@Walmsley2013mRNA; @Passini2013Computational; @Sadrieh2014Multiscale; @Bartolucci2014Linking; @Samanta2015Ca; @Pathmanathan2015; @cardone2016human; @Zhou2016; @Mahoney2016Connexin43; @Corsi2017Noninvasive; @Dutta2016; @Reilly2016; @Dutta2017; @Lyon2018;  @Xin2019];
the effects of realistic tissue structure on simulated cardiac electrical activity [@Walmsley2013Estimation; @Lekadir2014Effect; @Lekadir2016Statistically; @Zacur2017];
the sources and consequences of inter-subject electrophysiological variability [@Dutta2013Ionic; @Britton2013Experimentally; @Elkins2013Variability; @Walmsley2015Application; @Britton2017; @Muszkiewicz2018];
predicting the effects of drugs on cardiac activity, including safety assessment [@Beattie2013Evaluation; @Zemzemi2013Computational; @Wallman2014Computational; @Mirams2014Prediction; @Passini2014Late; @Cardone-Noott2014Computational; @Zemzemi2015Effects; @Moreno2015New; @Davies2016Recent; @Hill2016Computational; @Passini2016; @britton2017quantitative; @McMillan2017; @passini2017human; @Lim2018Role];
and the development of associated web-based tools [@Williams2015Web; @Cooper2016Cardiac; @Daly2018].
Other studies enabled by Cardiac Chaste have advanced the
methodologies for parameter identifiability and inference, model
selection and uncertainty quantification for in cardiac
electrophysiology models
[@Daly2015Hodgkin; @Mirams2016White; @Johnstone2016Uncertainty; @Daly2017Comparing];
and for the verification and efficient numerical simulation of such
models
[@Marsh2012Secrets;@Agudelo-Toro2013Computationally; @Pathmanathan2014Verification; @Corrado2016Stability; @Campos2016Lattice; @Spiteri2016Godunov; @Cervi2018HighOrder; @Cardone2018; @Green2019LUT].
The continuum-mechanics solvers in Chaste have been used for studies of dielectric elastomers [@Langham2018Modeling]; and our electro-mechanics code as used in [@Carapella2014Quantitative], has also been used to verify new numerical methods [@gurev2015high].

Work using Cardiac Chaste has also been published on mesh generation and model simulation in the area of gastric electrophysiology, in particular focusing on  interstitial cell of Cajal network structure and function [@Sathar2014Biophysically; @Gao2014Developmental; @Sathar2015Tissue; @Sathar2015Comparison; @Sathar2015Multiscale; @Gao2015Stochastic].

Publications enabled by Cell-based Chaste have focused on: the cellular mechanisms and dynamics of intestinal homeostasis and carcinogenesis [@Dunn2013Computational; @Hu2014Epidermal; @Baker2014Quantification; @Osborne2015Multiscale; @Dunn2016Combined; @Langlands2016Paneth; @Carroll2017Interkinetic; @Almet2018Multicellular; @Muraro2018TNF];
the mechanisms underlying vascular tumour growth and response to therapy in the Microvessel Chaste project [@Grogan2017Microvessel; @Grogan2017Predicting; @Grogan2018Importance];
the biomechanical characterization of skin lesions [@Franzetti2015Combined];
the organisation and proliferation of stem and pluripotent cells in development [@Atwell2015Mechano-logical; @Koke2014Computational; @Godwin2017Extended];
the dynamics of developing epithelial tissues [@Kursawe2015Capabilities; @Tetley2016Unipolar; @Abdullah2017Universal; @Finegan2018Tissue; @Waites2018Information];
the spread of sexually-transmitted infections [@Nelson2014STI-GMaS];
vascular remodelling [@Osborne2018Fully];
the similarities and differences between competing cell-based modelling approaches [@Figueredo2013On-lattice; @Davit2013Validity; @Fletcher2013Implementing; @Osborne2017Comparing];
the calibration and parameterisation of such models [@Cooper2013Connecting; @Kursawe2018Approximate];
and their efficient numerical solution [@Harvey2015Parallel; @Rubinacci2015Cognac; @Kursawe2017Impact; @Cooper2017Numerical].

Papers on Lung Chaste describe its use for patient-specific airway tree generation and flow modelling  [@bordas2015development; @soares2017evaluation; @burrowes2017combined].

# Acknowledgements

F.R.C. was supported by the Engineering and Physical Sciences Research Council [grant number EP/G03706X/1];
R.E.B. is a Royal Society Wolfson Research Merit Award holder, would like to thank the Leverhulme Trust for a Research Fellowship, and also acknowledges the BBSRC for funding [grant number BB/R000816/1];
B.D.E. was generously supported by the Wellcome Trust Institutional Strategic Support Award [grant number 204909/Z/16/Z];
A.G.F. was supported by a University of Sheffield Vice-Chancellor's Fellowship and the Biotechnology and Biological Sciences Research Council [grant number BB/R016925/1];
J.K.   was supported by the Engineering and Physical Sciences Research Council [grant number EP/N509711/1];
G.R.M. was supported by a Sir Henry Dale Fellowship jointly funded by the Wellcome Trust and Royal Society [Wellcome Trust grant number 101222/Z/13/Z] and a Wellcome Trust Senior Research Fellowship [grant number 212203/Z/18/Z];
A.B.O. was supported by a British Heart Foundation Intermediate Basic Science Fellowship [grant number FS/17/22/32644];
B.R. was supported by a Wellcome Trust Senior Research Fellowship [grant numbers 100246/Z/12/Z and 214290/Z/18/Z];
creation of the virtual airway structures was funded partially through the EU FP 7 AirPROM project.

# References
