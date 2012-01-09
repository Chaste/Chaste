/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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

    void TestBidomain3d() throw (Exception)
    {
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
