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


#ifndef _TESTMONODOMAINFITZHUGHNAGUMO_HPP_
#define _TESTMONODOMAINFITZHUGHNAGUMO_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>


#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"

#include "SimpleStimulus.hpp"

#include "FitzHughNagumo1961OdeSystem.hpp"
#include "ReplicatableVector.hpp"


class FhnEdgeStimulusCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
public:
    FhnEdgeStimulusCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-10.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        if (GetMesh()->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new FitzHughNagumo1961OdeSystem(mpSolver, mpStimulus);
        }
        else
        {
            return new FitzHughNagumo1961OdeSystem(mpSolver, mpZeroStimulus);
        }
    }
};

class TestMonodomainFitzHughNagumo : public CxxTest::TestSuite
{
public:

    void tearDown()
    {
        HeartConfig::Reset();
    }

    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating the left
    // edge.
    void TestMonodomainFitzHughNagumoWithEdgeStimulus( void ) throw (Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.01, 0.01));
        HeartConfig::Instance()->SetSimulationDuration(1.2); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        HeartConfig::Instance()->SetOutputDirectory("FhnWithEdgeStimulus");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainFhn_2dWithEdgeStimulus");

        FhnEdgeStimulusCellFactory cell_factory;

        // using the criss-cross mesh so wave propagates properly
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);


        monodomain_problem.Solve();


        /*
        * Test the top right node against the right one in the 1D case,
        * comparing voltage, and then test all the nodes on the right hand
        * side of the square against the top right one, comparing voltage.
        */
        bool need_initialisation = true;
        double probe_voltage=0.0;

        DistributedVector voltage = monodomain_problem.GetSolutionDistributedVector();
        need_initialisation = true;

        // Test the RHS of the mesh
        for (DistributedVector::Iterator node_index = voltage.Begin();
             node_index != voltage.End();
             ++node_index)
        {
            if (monodomain_problem.rGetMesh().GetNode(node_index.Global)->GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to
                //       to test the value of the nodes

                if (need_initialisation)
                {
                    probe_voltage = voltage[node_index];
                    need_initialisation = false;
                }
                else
                {
                    // Tests the final voltages for all the RHS edge nodes
                    // are close to each other.
                    // This works as we are using the 'criss-cross' mesh,
                    // the voltages would vary more with a mesh with all the
                    // triangles aligned in the same direction.

                    TS_ASSERT_DELTA(voltage[node_index], probe_voltage, 2e-4);
                }

                TS_ASSERT_DELTA(voltage[node_index], 0.139426, 2e-3);
             }
        }
    }
};
#endif //_TESTMONODOMAINFITZHUGHNAGUMO_HPP_
