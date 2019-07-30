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


#ifndef _TESTMONODOMAINPROBLEMLONG_HPP_
#define _TESTMONODOMAINPROBLEMLONG_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include <petscvec.h>
#include <vector>

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "ReplicatableVector.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "SimpleStimulus.hpp"

class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    unsigned mNodeNum;
public:
    PointStimulus2dCellFactory(unsigned nodeNum)
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-6000.0, 0.5)),
          mNodeNum(nodeNum)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        unsigned node_index = pNode->GetIndex();
        if (node_index == mNodeNum)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};


class TestMonodomainProblemLong : public CxxTest::TestSuite
{
public:

    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating in the
    // very centre of the mesh.
    // We run for 500 ms and then check that all the voltages at the final time
    // have returned to the resting potential of -84.5
    // test should take about 30mins (or less)
    //
    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainProblem2DWithPointStimulusInTheVeryCentreOfTheMesh( void )
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetPrintingTimeStep(1.0); //ms
        HeartConfig::Instance()->SetSimulationDuration(500); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem2dWithPointStimulusLong");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_2dWithPointStimulusLong");

        PointStimulus2dCellFactory cell_factory(60); // Central node

        MonodomainProblem<2> monodomain_problem(&cell_factory);

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        CheckMonoLr91Vars(monodomain_problem);

        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());

        /*
         * Test that corners are 'equal', and centres of sides.
         * Irregularities in which way the triangles are oriented make
         * this rather difficult, especially since the edges are sampled
         * during the upstroke.
         */

        // corners
        TS_ASSERT_DELTA(voltage_replicated[0], voltage_replicated[10],  0.1);
        TS_ASSERT_DELTA(voltage_replicated[0], voltage_replicated[110], 0.1);
        TS_ASSERT_DELTA(voltage_replicated[0], voltage_replicated[120], 0.1);

        // centres of edges
        TS_ASSERT_DELTA(voltage_replicated[5], voltage_replicated[55],  0.1);
        TS_ASSERT_DELTA(voltage_replicated[5], voltage_replicated[65],  0.1);
        TS_ASSERT_DELTA(voltage_replicated[5], voltage_replicated[115], 0.1);

        int num_nodes = monodomain_problem.rGetMesh().GetNumNodes();
        // test final voltages have returned to the resting potential
        for (int i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(voltage_replicated[i], -84.5, 1);
        }

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};
#endif //_TESTMONODOMAINPROBLEMLONG_HPP_
