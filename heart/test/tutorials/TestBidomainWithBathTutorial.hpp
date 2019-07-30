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
#ifndef TESTBIDOMAINWITHBATHTUTORIAL_HPP_
#define TESTBIDOMAINWITHBATHTUTORIAL_HPP_

/* HOW_TO_TAG Cardiac/Problem definition
 * Run bidomain simulations with a perfusing bath, and apply shocks using electrodes
 */

/*
 * = An example showing how to run a bidomain simulation for tissue contained in a perfusing bath =
 *
 * In this tutorial we show how the changes the need to be made when running a simulation of
 * cardiac tissue contained in a bath.
 *
 */

/* The usual headers are included */
#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"

/* HOW_TO_TAG Cardiac/Solver
 * Using specialised Backward Euler implementation to solve the cell models (allows for much larger timesteps)
 *
 * Cell models can be solved using a specialised (for cardiac cell models) Backward Euler
 * implementation, with again the code being automatically generated from the cellml files.
 * Backward Euler allows much greater ODE timesteps to be used. For cellml models provided with Chaste,
 * using Backward Euler is trivial: just change the .hpp included as follows, and the class name as below.
 *
 * If you want to use Backward Euler for your own CellML file, there is an intermediate step - see
 * https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/CodeGenerationFromCellML#GeneratingBackwardEulercellmodels
 * (note: requires Maple). (This step creates the .out files that are already present in heart/src/odes/cellml)
 *
 */
#include "LuoRudy1991BackwardEuler.hpp"
#include "PetscSetupAndFinalize.hpp"
/* This test will show how to load a mesh in the test and pass it into the problem,
 * for which the following includes are needed */
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
/* This header is needed for the sqrt function. */
#include <cmath>

/* Define the test */
class TestBidomainWithBathTutorial : public CxxTest::TestSuite
{
public: // Tests should be public!

    void TestWithBathAndElectrodes()
    {
        /* First, set the end time and output info. In this simulation
         * we'll explicitly read the mesh, alter it, then pass it
         * to the problem class, so we don't set the mesh file name.
         */
        HeartConfig::Instance()->SetSimulationDuration(3.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainTutorialWithBath");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        /* Bath problems seem to require decreased ODE timesteps.
         */
        HeartConfig::Instance()->SetOdeTimeStep(0.001);  //ms

        /* Use the {{{PlaneStimulusCellFactory}}} to define a set
         * of Luo-Rudy cells. We pass the stimulus magnitude as 0.0
         * as we don't want any stimulated cells.
         */
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellMLBackwardEuler,2> cell_factory(0.0);

        /*
         * Now, we load up a rectangular mesh (in triangle/tetgen format), done as follows,
         * using {{{TrianglesMeshReader}}}.  Note that we use a distributed mesh, so the data
         * is shared among processes if run in parallel.
         */
        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_400_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        /*
         * In most simulations there is one valid tissue identifier and one valid bath identifier
         * (for elements).
         * One of these can be assigned to an element with
         *  * {{{mesh.GetElement(i)->SetAttribute(HeartRegionCode::GetValidTissueId());}}}
         *  * {{{mesh.GetElement(i)->SetAttribute(HeartRegionCode::GetValidBathId());}}}
         *
         * If we want heterogeneous conductivities outside the heart (for example for torso and blood)
         * then we will need different identifiers:
         */
        std::set<unsigned> tissue_ids;
        static unsigned tissue_id=0;
        tissue_ids.insert(tissue_id);

        std::set<unsigned> bath_ids;
        static unsigned bath_id1=1;
        bath_ids.insert(bath_id1);
        static unsigned bath_id2=2;
        bath_ids.insert(bath_id2);

        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);

        /* In bath problems, each element has an attribute which must be set
         * to 0 (cardiac tissue) or 1 (bath). This can be done by having an
         * extra column in the element file (see the file formats documentation,
         * or for example
         * mesh/test/data/1D_0_to_1_10_elements_with_two_attributes.ele,
         * and note that the header in this file has 1 at the end to indicate that
         * the file defines an attribute for each element). We have read in a mesh
         * without this type of information set up, so we set it up manually,
         * by looping over elements and setting those more than 2mm from the centre
         * as bath elements (by default, the others are cardiac elements).
         */
        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            double x = iter->CalculateCentroid()[0];
            double y = iter->CalculateCentroid()[1];
            if (sqrt((x-0.05)*(x-0.05) + (y-0.05)*(y-0.05)) > 0.02)
            {
                if (y<0.05)
                {
                    //Outside circle on the bottom
                    iter->SetAttribute(bath_id1);
                }
                else
                {
                    //Outside circle on the top
                    iter->SetAttribute(bath_id2);
                }
            }
            else
            {
                //IDs default to 0, but we want to be safe
                iter->SetAttribute(tissue_id);
            }
        }

        /* HOW_TO_TAG Cardiac/Problem definition
         * Tell Chaste that a mesh has been modified
         *
         * Since we have modified the mesh by setting element attributes, we need to inform Chaste of this fact.
         * If we do not, problems will arise when [wiki:UserTutorials/CardiacCheckpointingAndRestarting checkpointing],
         * since the code that saves the simulation state will assume that it can just reuse the original mesh files,
         * and thus won't save the new element attributes.
         *
         * (Some mesh modifications, that use methods on the mesh class directly, will automatically record that
         * the mesh has been modified.  Since we're just modifying elements, this information isn't propagated at
         * present.)
         */
        mesh.SetMeshHasChangedSinceLoading();

        /*
         * The external conductivity can set two ways:
         *  * the default conductivity in the bath is set with {{{SetBathConductivity(double)}}}
         *  * heterogeneous overides can be set with {{{SetBathMultipleConductivities(std::map<unsigned, double> )}}}
         */

        HeartConfig::Instance()->SetBathConductivity(7.0);  //bath_id1 tags will take the default value (actually 7.0 is the default)
        std::map<unsigned, double> multiple_bath_conductivities;
        multiple_bath_conductivities[bath_id2] = 6.5;  // mS/cm

        HeartConfig::Instance()->SetBathMultipleConductivities(multiple_bath_conductivities);


        /* Now we define the electrodes. First define the magnitude of the electrodes
         * (ie the magnitude of the boundary extracellular stimulus), and the duration
         * it lasts for. Currently, electrodes switch on at time 0 and have constant magnitude
         * until they are switched off. (Note that this test has a small range of
         * magnitudes that will work, perhaps because the electrodes are close to the tissue).
         */
        // For default conductivities and explicit cell model -1e4 is under threshold, -1.4e4 too high - crashes the cell model
        // For heterogeneous conductivities as given, -1e4 is under threshold
        double magnitude = -14.0e3; // uA/cm^2
        double start_time = 0.0;
        double duration = 1; //ms

        /* Electrodes work in two ways: the first electrode applies an input flux, and
         * the opposite electrode can either be grounded or apply an equal and opposite
         * flux (ie an output flux). The `false` here indicates the second electrode
         * is not grounded, ie has an equal and opposite flux. The "0" indicates
         * that the electrodes should be applied to the bounding surfaces in the x-direction
         * (1 would be y-direction, 2 z-direction), which are X=0.0 and X=0.1 in the given mesh.
         * (This explains why the full mesh ought to be rectangular/cuboid - the nodes on
         * x=xmin and x=xmax ought to be form two surfaces of equal area.
         */
        HeartConfig::Instance()->SetElectrodeParameters(false, 0, magnitude, start_time, duration);

        /* Now create the problem class, using the cell factory and passing
         * in `true` as the second argument to indicate we are solving a bath
         * problem..
         */
        BidomainProblem<2> bidomain_problem( &cell_factory, true );

        /* ..set the mesh and electrodes.. */
        bidomain_problem.SetMesh(&mesh);

        /* ..and solve as before. */
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        /* The results can be visualised as before. '''Note:''' The voltage is only
         * defined at cardiac nodes (a node contained in ''any'' cardiac element), but
         * for visualisation and computation a 'fake' value of ZERO is given for the
         * voltage at bath nodes.
         *
         * Finally, we can check that an AP was induced in any of the cardiac
         * cells. We use a `ReplicatableVector` as before, and make sure we
         * only check the voltage at cardiac cells.
         */
        Vec solution = bidomain_problem.GetSolution(); // the Vs and phi_e's, as a PetSc vector
        ReplicatableVector solution_repl(solution);

        bool ap_triggered = false;
        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            if (HeartRegionCode::IsRegionTissue( iter->GetRegion() ))
            {
                if (solution_repl[2*iter->GetIndex()] > 0.0) // 2*i, ie the voltage for this node (would be 2*i+1 for phi_e for this node)
                {
                    ap_triggered = true;
                }
            }
        }
        TS_ASSERT(ap_triggered);
    }
};

#endif /*TESTBIDOMAINWITHBATHTUTORIAL_HPP_*/
