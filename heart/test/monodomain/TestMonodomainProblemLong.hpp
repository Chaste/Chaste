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

#include <ctime>

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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        if (node == mNodeNum)
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
