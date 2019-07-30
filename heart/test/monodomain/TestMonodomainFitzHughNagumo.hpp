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


#ifndef _TESTMONODOMAINFITZHUGHNAGUMO_HPP_
#define _TESTMONODOMAINFITZHUGHNAGUMO_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "SimpleStimulus.hpp"
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "ReplicatableVector.hpp"
#include "PetscSetupAndFinalize.hpp"


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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        if (pNode->GetPoint()[0] == 0.0)
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
    void TestMonodomainFitzHughNagumoWithEdgeStimulus( void )
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
