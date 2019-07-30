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


#ifndef TESTBIDOMAIN3D_HPP_
#define TESTBIDOMAIN3D_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "DistributedVector.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"


class TestBidomain3D :  public CxxTest::TestSuite
{
public:
    void TestBidomain3d()
    {
        HeartEventHandler::Reset();
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        HeartConfig::Instance()->SetSimulationDuration(4.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");
        HeartConfig::Instance()->SetOutputDirectory("Bidomain3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain3d");

        // Check the linear system can be solved to a low tolerance (in particular, checks the null space
        // stuff was implemented correctly
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-14);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> bidomain_cell_factory(-600.0*1000);

        BidomainProblem<3> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.Initialise();

        //bidomain_problem.SetNodeForAverageOfPhiZeroed(1330);
        //bidomain_problem.SetFixedExtracellularPotentialNodes(1330);

        bidomain_problem.Solve();

        Vec voltage=bidomain_problem.GetSolution();
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(voltage);

        /*
         * Test the top right node against the right one in the 1D case,
         * comparing voltage, and then test all the nodes on the right hand
         * face of the cube against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage=0;

        need_initialisation = true;

        // Test the RHF of the mesh
        for (AbstractTetrahedralMesh<3,3>::NodeIterator it = bidomain_problem.rGetMesh().GetNodeIteratorBegin();
             it != bidomain_problem.rGetMesh().GetNodeIteratorEnd();
             ++it)
        {
            if (it->GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to
                //       to test the value of the nodes

                if (need_initialisation)
                {
                    probe_voltage = voltage_replicated[2*it->GetIndex()];
                    need_initialisation = false;
                }
                else
                {
                    // the voltage at the end face varies a little because
                    // of drift due to the orientation of the tets in the mesh,
                    // hence the tolerance of 0.2
                    TS_ASSERT_DELTA(voltage_replicated[2*it->GetIndex()], probe_voltage, 0.2);
                }

                // if a 1D simulation is run for 4ms on the 0_1mm_10elements mesh
                // the result at the end node is 20.0755
                TS_ASSERT_DELTA(voltage_replicated[2*it->GetIndex()], 20.0755, 1.3);
            }
        }

        /*
         * HOW_TO_TAG Cardiac/Output
         * Collect and print timings to benchmark different parts of the cardiac code.
         *
         * N.B. You may also want to use HeartEventHandler::Reset() if you call these again.
         */
        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }

    ////////////////////////////////////////////////////////////
    // Compare Mono and Bidomain Simulations
    ////////////////////////////////////////////////////////////
    void TestCompareBidomainProblemWithMonodomain3D()
    {
        // the bidomain equations reduce to the monodomain equations
        // if sigma_e is infinite (equivalent to saying the extra_cellular
        // space is grounded. sigma_e is set to be very large here:
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(17500, 17500, 17500));
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");
        HeartConfig::Instance()->SetOutputDirectory("Monodomain3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain3d");

        ///////////////////////////////////////////////////////////////////
        // monodomain
        ///////////////////////////////////////////////////////////////////
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-600.0*1000);
        MonodomainProblem<3> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        ///////////////////////////////////////////////////////////////////
        // bidomain
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("Bidomain3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain3d");

        BidomainProblem<3> bidomain_problem( &cell_factory );
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector monodomain_voltage = monodomain_problem.GetSolutionDistributedVector();
        DistributedVector bidomain_solution = bidomain_problem.GetSolutionDistributedVector();
        DistributedVector::Stripe bidomain_voltage(bidomain_solution,0);
        DistributedVector::Stripe extracellular_potential(bidomain_solution,1);
        for (DistributedVector::Iterator index = bidomain_solution.Begin();
             index != bidomain_solution.End();
             ++index)
        {
            TS_ASSERT_DELTA(monodomain_voltage[index], bidomain_voltage[index], 0.5);
            TS_ASSERT_DELTA(extracellular_potential[index], 0, 1.0);
        }
    }
};

#endif /*TESTBIDOMAIN3D_HPP_*/
