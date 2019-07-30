/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTRUNNINGBIDOMAINSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGBIDOMAINSIMULATIONSTUTORIAL_HPP_

/*
 * HOW_TO_TAG Cardiac/Problem definition
 * Set up and run basic bidomain simulations
 */

/*
 * = An example showing how to run bidomain simulations =
 *
 * == Introduction ==
 *
 * In this tutorial we show how Chaste is used to run a standard bidomain simulation.
 * Note that monodomain simulations are run very similarly.
 *
 * The first thing that needs to be done, when writing any Chaste test,
 * is to include the following header.
 */

#include <cxxtest/TestSuite.h>
/* The main class to be used for running bidomain simulations is {{{BidomainProblem}}}. */
#include "BidomainProblem.hpp"
/* The type of intracellular stimulus we'll apply. */
#include "SimpleStimulus.hpp"
/* All tests which run cardiac simulations (which use Petsc) should include
 * {{{PetscSetupAndFinalize.hpp}}}.  This class ensures that {{{PetscInitialise()}}}
 * is called with the appropriate arguments before any tests in the suite are run. */
#include "PetscSetupAndFinalize.hpp"

/* HOW_TO_TAG Cardiac/Problem definition
 * Use different cell models, defined using CellML files
 */

/* The above files are contained in the source release and can be located and studied. Cardiac cell
 * models are different: the C++ code is automatically generated from CellML files. To use a particular
 * CellML file, place it in `heart/src/odes/cellml` (there are several in here already). If the CellML
 * is called `<CELLMODEL>.cellml`, a file `<CELLMODEL>.hpp` will be automatically generated, which will define
 * a class called `Cell<CELLMODEL>FromCellML`. So to use a particular cell model in a tissue simulation,
 * given the CellML, you just have to do two things: include this `.hpp` file, and then use the class.
 * For example, we will use the !LuoRudy1991 model, so we have to include the following, and
 * later on use {{{CellLuoRudy1991FromCellML}}} as the cell model class.
 * See ["ChasteGuides/CodeGenerationFromCellML"] for more information on this process.
 */
#include "LuoRudy1991.hpp"


/*
 * == Defining a cell factory ==
 *
 * All mono/bidomain simulations need a ''cell factory'' as input. This is a class
 * which tells the problem class what type of cardiac cells to create. The cell-factory
 * class has to inherit from {{{AbstractCardiacCellFactory<DIM>}}}, which means it must
 * implement the method {{{CreateCardiacCellForTissueNode(Node<DIM>*)}}}, which returns
 * a pointer to an {{{AbstractCardiacCell}}}. Note, some concrete cell factories have
 * been defined, such as the {{{PlaneStimulusCellFactory}}} (see later tutorials), which
 * could be used in the simulation, but for completeness we create our own cell factory in
 * this test. For complicated problems with, say, heterogeneous cell types or particular stimuli,
 * a new cell factory will have to be defined by the user for their particular problem.
 *
 * This cell factory is a simple cell factory where every cell is a Luo-Rudy 91 cell,
 * and only the cell at position (0,0) is given a non-zero stimulus.
 */
class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
/* Declare (smart) pointer to a {{{SimpleStimulus}}} for the cell which is stimulated.
 * Note that {{{AbstractCardiacCellFactory}}} also has as protected members: {{{mpZeroStimulus}}}
 * of type {{{boost::shared_ptr<ZeroStimulus>}}}; {{{mpMesh}}}, a pointer to the mesh used (the problem
 * class will set this before it calls {{{CreateCardiacCellForTissueNode}}}, so it can be used
 * in that method); {{{mTimestep}}}, a double (see below); and {{{boost::shared_ptr<mpSolver>}}}
 * a forward euler ode solver (see below). */
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    /* Our contructor takes in nothing. It calls the constructor of {{{AbstractCardiacCellFactory}}}
     * and we also initialise the stimulus to have magnitude -500000 uA/cm^3 and duration 0.5 ms.
     */
    PointStimulus2dCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-5e5, 0.5))
    {
    }

    /* Now we implement the pure method which needs to be implemented. We return
     * a LR91 cell for each node, with the nodes in a 0.2mm block given the non-zero stimulus,
     * and all other nodes given the zero stimulus. Note that we use {{{mpMesh}}},
     * {{{mTimestep}}}, {{{mpZeroStimulus}}} and {{{mpSolver}}} which are all
     * members of the base class. The timestep and solver are defined in the base
     * class just so that the user doesn't have to create them here. */
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        if (x<0.02+1e-6 && y<0.02+1e-6) // ie if x<=0.02 and y<=0.02 (and we are assuming here x,y>=0).
        {
            /* Create a LR91 cell with the non-zero stimulus. This is a volume stimulus, ie
             * the function on the right-hand side of the first of the two bidomain equations.
             * An equal and opposite extra-cellular stimulus is implicitly enforced by the code,
             * which corresponds to having zero on the right-hand side of the second of the
             * bidomain equations.
             */
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            /* The other cells have zero stimuli. */
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    /* We have no need for a destructor, since the problem class deals with deleting the cells. */
};

/*
 * == Running the bidomain simulation ==
 *
 * Now we can define the test class, which must inherit from {{{CxxTest::TestSuite}}}
 * as described in the writing basic tests tutorial. */
class TestRunningBidomainSimulationsTutorial : public CxxTest::TestSuite
{
/* Tests should be public... */
public:
    /* Define the test. Note the {{{}}} - without this exception messages
     * might not get printed out.
     */
    void TestSimpleSimulation()
    {
        /* The {{{HeartConfig}}} class is used to set various parameters (see the main ChasteGuides page
         * for information on default parameter values. Parameters in this file can be re-set
         * with {{{HeartConfig}}} if the user wishes, and other parameters such as end time must be set
         * using {{{HeartConfig}}}. Let us begin by setting the end time (in ms), the mesh to use, and the
         * output directory and filename-prefix. Note that the spatial units in cardiac Chaste is CENTIMETRES,
         * so that mesh 2D_0_to_1mm_800_elements is a mesh over [0,0.1]x[0,0.1].
         */
        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_800_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainTutorial");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        /* There is an alternate method of loading a mesh that can be seen in [wiki:UserTutorials/Monodomain3dExample Monodomain3dExample],
         *  using `DistributedTetrahedralMesh`.
         *
         * It is possible to over-ride the default visualisation output (which is done during simulation
         * post-processing).
         */
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithCmgui(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        /* If the mesh is a DistributedTetrahedralMesh then we can use parallel VTK files (.pvtu)*/
        //HeartConfig::Instance()->SetVisualizeWithParallelVtk(true);

        /* Next, we have to create a cell factory of the type we defined above. */
        PointStimulus2dCellFactory cell_factory;

        /* Now we create a problem class using (a pointer to) the cell factory. */
        BidomainProblem<2> bidomain_problem( &cell_factory );

        /* This is enough setup to run a simulation: we could now call {{{Initialise()}}}
         * and {{{Solve()}}} to run... */
        // bidomain_problem.Initialise();
        // bidomain_problem.Solve();

        /* ..however, instead we show how to set a few more parameters. To set the conductivity values
         *  in the principal fibre, sheet and normal directions do the following.
         * Note that {{{Create_c_vector}}} is just a helper method for creating a {{{c_vector<double,DIM>}}}
         * of the correct size (2, in this case). Make sure these methods are called before
         * {{{Initialise()}}}.
         */
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.19));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(6.2, 2.4));

        /* This is how to reset the surface-area-to-volume ratio and the capacitance.
         * (Here, we are actually just resetting them to their default values). */
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400); // 1/cm
        HeartConfig::Instance()->SetCapacitance(1.0); // uF/cm^2

        /* This is how to set the ode timestep (the timestep used to solve the cell models)
         * the pde timestep (the timestep used in solving the bidomain PDE), and the
         * printing timestep (how often the output is written to file). The defaults are
         * all 0.01, here we increase the printing timestep.
         */
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);

        /* Now we call {{{Initialise()}}}... */
        bidomain_problem.Initialise();

        /* Now we call Solve() to run the simulation. The output will be written to
         * `/tmp/$USER/testoutput/BidomainTutorial` in HDF5 format.  The
         * output will also be converted to selected visualiser formats at the end of the simulation.
         * Note that if you want to view the progress of longer simulations
         * go to the the output directory and look at the file
         * {{{progress_status.txt}}}, which will say the percentage of the
         * simulation run. */
        bidomain_problem.Solve();

        /*
         * == Examining the output ==
         * In order to visualise the results, go to one of the sub-folders
         *  * `/tmp/$USER/testoutput/BidomainTutorial/output` for Meshalyzer
         *  * `/tmp/$USER/testoutput/BidomainTutorial/cmgui_output` for Cmgui
         *  * `/tmp/$USER/testoutput/BidomainTutorial/vtk_output` for Paraview (VTK)
         * where you should find the geometric mesh data and simulation output.
         *
         * Please see ChasteGuides/VisualisationGuides for details of using !Meshalyzer/Cmgui/Paraview.
         *
         * Note: the easiest way to look at the resultant voltage values from the code
         * (for the last timestep - the data for the previous timesteps is written to file
         * but not retained) is to use a {{{ReplicatableVector}}}.
         * {{{bidomain_problem.GetSolution())}}} returns a !PetSc vector
         * of the form (V_0, phi_0, V_1, phi_e_1, ... V_n, phi_e_n), and we can create a
         * {{{ReplicatableVector}}} for easy access to this !PetSc vector's data.
         * (This won't be very efficient with huge problems in parallel - the next tutorial
         * will mention how to do parallel access).
         */
        ReplicatableVector res_repl(bidomain_problem.GetSolution());
        for (unsigned i=0; i<res_repl.GetSize(); i++)
        {
        //    std::cout << res_repl[i] << "\n";
        }

        /* Behind the scenes there are some logging routines which find out how much time
         * has been spent in the major parts of the code (solving ODEs, assembling matrices etc.)
         * The logging routines are in {{{HeartEventHandler}}} which is enabled by default.
         * If you think this is getting in the way, you can turn it off at the top of your test with
         * {{{HeartEventHandler::Disable()}}}.
         * In this test, we want to get information out of the {{{HeartEventHandler}}}.
         */
        /* {{{Headings()}}} prints a single (very long) line reminding us what catagories of events are being instrumented.*/
        HeartEventHandler::Headings();
        /* {{{Report()}}} prints a single line with times spent in each catagory.  When run in parallel it prints one line of times per process and also lines for average
         * and maximum times.  (This can be useful if you need to identify a load imbalance.)*/
        HeartEventHandler::Report();
    }
};

#endif /*TESTRUNNINGBIDOMAINSIMULATIONSTUTORIAL_HPP_*/
