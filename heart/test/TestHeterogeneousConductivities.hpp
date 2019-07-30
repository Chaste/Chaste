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
#ifndef TESTHETEROGENEOUSCONDUCTIVITIES_HPP_
#define TESTHETEROGENEOUSCONDUCTIVITIES_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <fstream>

#include "BidomainProblem.hpp"
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "TrianglesMeshReader.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "LuoRudy1991.hpp"
#include "PetscSetupAndFinalize.hpp"

using std::ofstream;

/* test class*/
class TestHeterogeneousConductivities : public CxxTest::TestSuite
{
public:
    void TestSimpleSimulation()
    {
        /*Simulation parameters*/
        HeartConfig::Instance()->SetSimulationDuration(0.7); //ms (falls over after this)
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-6);
        //HeartConfig::Instance()->SetOdeTimeStep(0.01);

        const double width = 0.1;
        const double height = 0.1;
        const double depth = 0.1;

        const unsigned num_elem_x = 8;
        const double space_step = width/num_elem_x;

        /* Make the mesh*/
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructRegularSlabMesh(space_step, width, height, depth);

        /*Create a cell factory of the type we defined above. */
        GeneralPlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(num_elem_x, width);

        /* monodomain problem class using (a pointer to) the cell factory */
        BidomainProblem<3> problem( &cell_factory );
        problem.SetMesh(&mesh);

        /*
        * HOW_TO_TAG Cardiac/Problem definition
        * Set discrete '''cuboid''' areas to have heterogeneous (intra- and/or extra-cellular) conductivity tensors.
        */
        std::vector<ChasteCuboid<3> > input_areas;
        std::vector< c_vector<double,3> > intra_conductivities;
        std::vector< c_vector<double,3> > extra_conductivities;
        ChastePoint<3> corner_a(width/2, 0, 0);
        ChastePoint<3> corner_b(width, height, depth);

        input_areas.push_back(ChasteCuboid<3> (corner_a, corner_b));
        //within the cuboid
        intra_conductivities.push_back( Create_c_vector(0.1, 0.1, 0.1) );
        extra_conductivities.push_back( Create_c_vector(0.0, 0.0, 0.0) );
        //This test should *fail* if you comment out the following line
        //(which blocks conductivity on the RHS of the slab).
        HeartConfig::Instance()->SetConductivityHeterogeneities(input_areas, intra_conductivities, extra_conductivities);

        //elsewhere
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.2, 1.2, 1.2));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(1.2, 1.2, 1.2));

        /* set  parameters*/
        // HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        // HeartConfig::Instance()->SetCapacitance(1.0);

         /* Output Directory and prefix (for the hdf5 file), relative to CHASTE_TEST_OUTPUT*/
        HeartConfig::Instance()->SetOutputDirectory("slab_results_het_halfcond");
        HeartConfig::Instance()->SetOutputFilenamePrefix("Slab_small");

        /* Initialise the problem*/
        problem.Initialise();

        /* Solve the PDE monodomain equaion*/
        problem.Solve();

        ReplicatableVector voltage_replicated(problem.GetSolution());
        TS_ASSERT_EQUALS(mesh.GetNumNodes() * 2, voltage_replicated.GetSize());
        unsigned lo, hi;
        lo = mesh.GetDistributedVectorFactory()->GetLow();
        hi = mesh.GetDistributedVectorFactory()->GetHigh();

        for (unsigned i=lo; i<hi; i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            if (x<width/2)
            {
                 //Left side is stimulated
                 TS_ASSERT_LESS_THAN(-71.0,voltage_replicated[2 * i]);
            }
            else if (x>width/2)
            {
                //Right side is blocked
                TS_ASSERT_LESS_THAN(voltage_replicated[2 * i],-82.0);
            }
        }
    }
};

#endif /*TESTHETEROGENEOUSCONDUCTIVITIES_HPP_*/
