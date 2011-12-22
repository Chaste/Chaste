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


#ifndef _TESTNEUMANNSTIMULUS_HPP_
#define _TESTNEUMANNSTIMULUS_HPP_


#include <cxxtest/TestSuite.h>
#include <vector>
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "ReplicatableVector.hpp"
#include "SimpleStimulus.hpp"
#include "StimulusBoundaryCondition.hpp"
#include "ConstBoundaryCondition.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestNeumannStimulus : public CxxTest::TestSuite
{
public:
    void tearDown()
    {
        HeartConfig::Reset();
    }

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    void TestMonodomainConstantStimulus() throw(Exception)
    {
        // this parameters are a bit arbitrary, and chosen to get a good spread of voltages
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75));
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoNeumannConst");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");

        ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1*1.75/0.0005);

        // create boundary conditions container
        boost::shared_ptr<BoundaryConditionsContainer<1,1,1> > p_bcc(new BoundaryConditionsContainer<1,1,1>);
        ConstBoundaryCondition<1>* p_bc_stim = new ConstBoundaryCondition<1>(2*1.75/0.0005);

        // get mesh
        AbstractTetrahedralMesh<1,1> &mesh = monodomain_problem.rGetMesh();
        // loop over boundary elements
        AbstractTetrahedralMesh<1, 1>::BoundaryElementIterator iter;
        iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            // if the element is on the left of the mesh, add a stimulus to the bcc
            if (((*iter)->GetNodeLocation(0))[0]==0.0)
            {
                p_bcc->AddNeumannBoundaryCondition(*iter, p_bc_stim);
            }
            iter++;
        }

        // pass the bcc to the monodomain problem
        monodomain_problem.SetBoundaryConditionsContainer(p_bcc);

        monodomain_problem.Solve();

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol=5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 94.6426, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 49.7867, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 30.5954, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 21.6782, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -33.9983, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -52.2396, atol);

    }

    void TestMonodomainSquareWaveStimulus() throw(Exception)
    {
        // this parameters are a bit arbitrary, and chosen to get a good spread of voltages
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75));
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoNeumannSquare");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");

        ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1*1.75/0.0005);

        // create boundary conditions container
        boost::shared_ptr<BoundaryConditionsContainer<1,1,1> > p_bcc(new BoundaryConditionsContainer<1,1,1>);
        SimpleStimulus stim(4*1.75/0.0005, 0.5);
        StimulusBoundaryCondition<1>* p_bc_stim = new StimulusBoundaryCondition<1>(&stim);

        // get mesh
        AbstractTetrahedralMesh<1,1> &mesh = monodomain_problem.rGetMesh();
        // loop over boundary elements
        AbstractTetrahedralMesh<1, 1>::BoundaryElementIterator iter;
        iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            // if the element is on the left of the mesh, add a stimulus to the bcc
            if (((*iter)->GetNodeLocation(0))[0]==0.0)
            {
                p_bcc->AddNeumannBoundaryCondition(*iter, p_bc_stim);
            }
            iter++;
        }

        // pass the bcc to the monodomain problem
        monodomain_problem.SetBoundaryConditionsContainer(p_bcc);

        monodomain_problem.Solve();

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol=8e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 22.4940, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 22.6008, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 23.3054, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.4932, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], 14.5184, atol);
        TS_ASSERT_DELTA(voltage_replicated[10],3.7081, atol);
    }

/// Note: Neumann intracellular stimuli won't work with bidomain (no solution to PDEs) (unless nodes are grounded)


};

#endif //_TESTNEUMANNSTIMULUS_HPP_
