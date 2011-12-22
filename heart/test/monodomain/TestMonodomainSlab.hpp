/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef _TESTMONODOMAINSLAB_HPP_
#define _TESTMONODOMAINSLAB_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "petscvec.h"
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "ReplicatableVector.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"

class TestMonodomainSlab : public CxxTest::TestSuite
{

public:

    // Solve on a 3D 1mm by 1mm by 1mm mesh (space step = 0.1mm), stimulating
    // the left face.
    // Should behave like the 1D case, extrapolated.
    // See also TestMonodomainSlab.hpp
    void TestMonodomain3DWithFaceStimulus( void ) throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetSimulationDuration(4); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoDg03dWithFaceStimulus");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_3dWithFaceStimulus");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-600.0*1000);

        MonodomainProblem<3> monodomain_problem(&cell_factory);

        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars(monodomain_problem);

        /*
         * Test the top right node against the right one in the 1D case,
         * comparing voltage, and then test all the nodes on the right hand
         * face of the cube against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage=0.0;
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        need_initialisation = true;

        // Test the RHF of the mesh
        for (unsigned i = 0; i < monodomain_problem.rGetMesh().GetNumNodes(); i++)
        {
            if (monodomain_problem.rGetMesh().GetNode(i)->GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to
                //       to test the value of the nodes

                if (need_initialisation)
                {
                    probe_voltage = voltage_replicated[i];
                    need_initialisation = false;
                }
                else
                {
                    TS_ASSERT_DELTA(voltage_replicated[i], probe_voltage, 1.1);
                    // std::cout << "y=" << monodomain_problem.mMesh.GetNode(i)->GetPoint()[1] << std::endl;
                }

                // Check against 1d case - the TestMonodomainDg01D test is run
                // for 4ms

                TS_ASSERT_DELTA(voltage_replicated[i], 21.3, 1);

            }
        }
    }
};


#endif //_TESTMONODOMAINSLAB_HPP_
