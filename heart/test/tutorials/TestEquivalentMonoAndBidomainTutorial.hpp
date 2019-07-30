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
#ifndef TESTEQUIVALENTMONOANDBIDOMAINTUTORIAL_HPP_
#define TESTEQUIVALENTMONOANDBIDOMAINTUTORIAL_HPP_

/* [[PageOutline]]
 *
 * = How to run a Bidomain simulation and its equivalent Monodomain reduction =
 *
 * == Introduction ==
 *
 * In this tutorial we show how Chaste is used to run a standard mono and a standard bidomain simulation.
 * With equivalent parameters so that the bidomain could be reduced to the monodomain case.
 *
 * The bulk of this tutorial is the same as UserTutorials/RunningBidomainSimulations, so for details of each line see that page.
 */

#include <cxxtest/TestSuite.h>
// The main classes to be used for running simulations.
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "SimpleStimulus.hpp"
///* All tests which run cardiac simulations (which use Petsc) should include
// * {{{PetscSetupAndFinalize.hpp}}}.  This class ensures that {{{PetscInitialise()}}}
// * is called with the appropriate arguments before any tests in the suite are run. */
#include "PetscSetupAndFinalize.hpp"

///* The above files are contained in the source release and can be located and studied. Cardiac cell
// * models are different: the C++ code is automatically generated from cellml files. To use a particular
// * cellml file, place it in `heart/src/odes/cellml` (there are several in here already). If the cellml
// * is called `<CELLMODEL>.cellml`, a file `<CELLMODEL>.hpp` will be automatically generated, which will define
// * a class called `Cell<CELLMODEL>FromCellML`. So to use a particular cell model in a tissue simulation,
// * given the cellml, you just have to do two things: include this `.hpp` file, and then use the class.
// * For example, we will use the !LuoRudy1991 model, so we have to include the following, and
// * later on use {{{CellLuoRudy1991FromCellML}}} as the cell model class.
// * See ["ChasteGuides/CodeGenerationFromCellML"] for more information on this process.
// */
#include "LuoRudy1991.hpp"


// * == Defining a cell factory ==
// *
// * All mono/bidomain simulations need a ''cell factory'' as input. This is a class
// * which tells the problem class what type of cardiac cells to create. The cell-factory
// * class has to inherit from {{{AbstractCardiacCellFactory<DIM>}}}, which means it must
// * implement the method {{{CreateCardiacCellForTissueNode(Node<DIM>*)}}}, which returns
// * a pointer to an {{{AbstractCardiacCell}}}. Note, some concrete cell factories have
// * been defined, such as the {{{PlaneStimulusCellFactory}}} (see later tutorials), which
// * could be used in the simulation, but for completeness we create our own cell factory in
// * this test. For complicated problems with, say, heterogeneous cell types or particular stimuli,
// * a new cell factory will have to be defined by the user for their particular problem.
// *
// * This cell factory is a simple cell factory where every cell is a Luo-Rudy 91 cell,
// * and only the cell at position (0) is given a non-zero stimulus.
// */
class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
///* Declare (smart) pointer to a {{{SimpleStimulus}}} for the cell which is stimulated.
// * Note that {{{AbstractCardiacCellFactory}}} also has as protected members: {{{mpZeroStimulus}}}
// * of type {{{boost::shared_ptr<ZeroStimulus>}}}; {{{mpMesh}}}, a pointer to the mesh used (the problem
// * class will set this before it calls {{{CreateCardiacCellForTissueNode}}}, so it can be used
// * in that method); {{{mTimestep}}}, a double (see below); and {{{boost::shared_ptr<mpSolver>}}}
// * a forward euler ode solver (see below). */
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    // Our contructor takes in nothing. It calls the constructor of {{{AbstractCardiacCellFactory}}}
    // and we also initialise the stimulus to have magnitude -500000 uA/cm^3 and duration 0.5 ms.
    PointStimulus2dCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-5e5, 0.5))
    {
    }

//     * Now we implement the pure method which needs to be implemented. We return
//     * a LR91 cell for each node, with the nodes in a 0.2mm block given the non-zero stimulus,
//     * and all other nodes given the zero stimulus. Note that we use {{{mpMesh}}},
//     * {{{mTimestep}}}, {{{mpZeroStimulus}}} and {{{mpSolver}}} which are all
//     * members of the base class. The timestep and solver are defined in the base
//     * class just so that the user doesn't have to create them here.
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        if (x<0.02+1e-6 && y<0.02+1e-6)
        {
//            * Create a LR91 cell with the non-zero stimulus. This is a volume stimulus, ie
//            * the function on the right-hand side of the first of the two bidomain equations.
//            * An equal and opposite extra-cellular stimulus is implicitly enforced by the code,
//            * which corresponds to having zero on the right-hand side of the second of the
//            * bidomain equations.
//            *
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            // The other cells have zero stimuli.
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
    // We have no need for a destructor, since the problem class deals with deleting the cells.
};

///*
// * == Running the mono- and bi-domain simulations ==
// *
// * Now we can define the test class, which must inherit from {{{CxxTest::TestSuite}}}
// * as described in the writing basic tests tutorial. */
class TestEquivalentMonoAndBidomainTutorial : public CxxTest::TestSuite
{
// Tests should be public...
public:
    // Define the test. Note the {{{}}} - without this exception messages
    // might not get printed out.
    void TestCompareMonoAndBidomain()
    {
        /* The {{{HeartConfig}}} class is used to set various parameters (see the main ChasteGuides page
         * for information on default parameter values.
         *
         * See UserTutorials/RunningBidomainSimulations for more details.
         */
        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_800_elements");
        HeartConfig::Instance()->SetOutputDirectory("EquivalentMonoAndBidomainTutorial");

        /* This is how to reset the surface-area-to-volume ratio and the capacitance.
         * (Here, we are actually just resetting them to their default values). */
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400); // 1/cm
        HeartConfig::Instance()->SetCapacitance(1.0); // uF/cm^2

        // This is how to set the ode timestep (the timestep used to solve the cell models)
        // the pde timestep (the timestep used in solving the bidomain PDE), and the
        // printing timestep (how often the output is written to file). The defaults are
        // all 0.01, here we increase the printing timestep.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);

        // Next, we have to create a cell factory of the type we defined above.
        PointStimulus2dCellFactory cell_factory;

        /*
         * == Setting Bidomain Conductivities ==
         */
        c_vector<double,2> intracellular_conductivities = Create_c_vector(1.75, 0.19);
        c_vector<double,2> extracellular_conductivities = Create_c_vector(7, 0.76);

        ReplicatableVector* p_bidomain_results;
        // Bidomain
        {
            HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_results");

            // Now we create a problem class using (a pointer to) the cell factory.
            BidomainProblem<2> bidomain_problem( &cell_factory );

            /*
             * Here we have conductivities that can be expressed as sigma_i = scalar * sigma_e.
             *
             * Then this is a special case, in which the bidomain equations can be reduced to the monodomain
             * equation. They are exactly equivalent.
             *
             * In this test we check that Chaste says this too when doing both simulations.
             *
             * This is how we set intracellular and extracellular conductivity directions.
             * Note that the extracellular is a constant scaling (4x) of the intracellular,
             * and hence a reduction to the monodomain equation can be made.
             *
             * For more information on this see e.g. Keener & Sneyd, Mathematical Physiology textbook.
             */
            HeartConfig::Instance()->SetIntracellularConductivities(intracellular_conductivities);
            HeartConfig::Instance()->SetExtracellularConductivities(extracellular_conductivities);

            /* Initialise and solve as normal */
            bidomain_problem.Initialise();
            bidomain_problem.Solve();

            /*
             * NB: the easiest way to look at the resultant voltage values from the code
             * (for the last timestep - the data for the previous timesteps is written to file
             * but not retained) is to use a {{{ReplicatableVector}}}.
             * {{{bidomain_problem.GetSolution())}}} returns a !PetSc vector
             * of the form (V_0, phi_0, V_1, phi_e_1, ... V_n, phi_e_n), and we can create a
             * {{{ReplicatableVector}}} for easy access to this !PetSc vector's data.
             * (This won't be very efficient with huge problems in parallel - the next tutorial
             * will mention how to do parallel access).
             */
            p_bidomain_results = new ReplicatableVector(bidomain_problem.GetSolution());
        }

        ReplicatableVector* p_monodomain_results;
        // Monodomain
        {
            /* HOW_TO_TAG Cardiac/Problem definition
             * Do a monodomain simulation with equivalent conductivities to a bidomain simulation.
             */

            /*
             * == Reduction to Monodomain  ==
             *
             * We now work out the equivalent conductivity to use in the monodomain. Note that
             * we set this with the `SetIntracellularConductivities()` method. Which is perhaps
             * better described as "Set first domain conductivities".
             *
             * So we calculate the equivalent conductivity according to (elementwise)
             *
             * sigma_monodomain = sigma_i sigma_e / (sigma_i + sigma_e)
             *
             */
            c_vector<double,2> monodomain_conductivities;

            // Just a little check that this reduction is valid in case you copy and paste this code!
            double x_scaling = extracellular_conductivities[0] / intracellular_conductivities[0];
            double y_scaling = extracellular_conductivities[1] / intracellular_conductivities[1];
            TS_ASSERT_DELTA(x_scaling, y_scaling, 1e-12);

            for (unsigned dim=0; dim<2; dim++)
            {
                monodomain_conductivities[dim] = intracellular_conductivities[dim]*extracellular_conductivities[dim]
                                               / (intracellular_conductivities[dim] + extracellular_conductivities[dim]);
            }

            HeartConfig::Instance()->SetIntracellularConductivities(monodomain_conductivities);

            /* Now we create a monodomain problem class in exactly the same way as bidomain above */
            HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain_results");
            MonodomainProblem<2> monodomain_problem( &cell_factory );
            monodomain_problem.Initialise();
            monodomain_problem.Solve();
            p_monodomain_results = new ReplicatableVector(monodomain_problem.GetSolution());
        }

        /* The bidomain solution includes extracellular (phi_e) so should be twice as big as monodomain solution.*/
        TS_ASSERT_EQUALS(p_bidomain_results->GetSize(),2*p_monodomain_results->GetSize());

        /* We then check that the voltage at each node at the end of the simulation is the same
         * whether we did a bidomain simulation, or the equivalent monodomain simulation.*/
        for (unsigned i=0; i<p_monodomain_results->GetSize(); i++)
        {
            TS_ASSERT_DELTA((*p_monodomain_results)[i], (*p_bidomain_results)[2u*i], 1e-6);
        }

        delete p_monodomain_results;
        delete p_bidomain_results;
    }
};

#endif /*TESTEQUIVALENTMONOANDBIDOMAINTUTORIAL_HPP_*/
