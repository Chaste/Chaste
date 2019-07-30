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

#ifndef TESTMONODOMAINPURKINJEPROBLEMLONG_HPP_
#define TESTMONODOMAINPURKINJEPROBLEMLONG_HPP_

#include <cxxtest/TestSuite.h>

#include "MonodomainPurkinjeVolumeAssembler.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "MixedDimensionMesh.hpp"
#include "PetscTools.hpp"
#include "LuoRudy1991.hpp"
#include "DiFrancescoNoble1985.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainPurkinjeCableAssembler.hpp"
#include "PetscMatTools.hpp"
#include "MonodomainPurkinjeSolver.hpp"
#include "MonodomainSolver.hpp"
#include "TetrahedralMesh.hpp"
#include "MonodomainProblem.hpp"
#include "MonodomainPurkinjeProblem.hpp"

class NonPurkinjeCellFactory3d : public AbstractCardiacCellFactory<3>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    NonPurkinjeCellFactory3d()
        : AbstractCardiacCellFactory<3>(),
          mpStimulus(new SimpleStimulus(-1e6, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        ChastePoint<3> location = pNode->GetPoint();

        if (fabs(location[2])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};



class PurkinjeCellFactory3d : public AbstractPurkinjeCellFactory<3> //inherits from different base
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PurkinjeCellFactory3d()
        : AbstractPurkinjeCellFactory<3>(),
          mpStimulus(new SimpleStimulus(-1e6, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        ChastePoint<3> location = pNode->GetPoint();

        if (fabs(location[2])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    AbstractCardiacCell* CreatePurkinjeCellForTissueNode(Node<3>* pNode,
                                                         AbstractCardiacCellInterface* pCardiacCell)
    {
        ChastePoint<3> location = pNode->GetPoint();

        if (fabs(location[2])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};


class TestMonodomainPurkinjeProblemLong : public CxxTest::TestSuite
{
public:

    // Solve on a cylinder (fibre running down the centre) - run with just a monodomain problem, then
    // run with
    void TestIn3d()
    {
        // Settings common to both problems
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(1.0);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        // make both problems use the same surface area
        HeartConfig::Instance()->SetPurkinjeSurfaceAreaToVolumeRatio(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());


        ReplicatableVector soln_repl;
        ReplicatableVector soln_mono_repl;

        // Need to load the mesh ourselves at present, as AbstractCardiacProblem::Initialise will assume a DistributedTetrahedralMesh.
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/cylinder_refined");
        TrianglesMeshReader<3,3> reader(mesh_base);
        MixedDimensionMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        mesh.Scale(0.1, 0.1, 0.1); // so spatial resolution of mesh is appropriate for cardiac simulation

        {
            HeartEventHandler::Reset();

            // Set up normal monodomain
            NonPurkinjeCellFactory3d cell_factory_for_just_monodomain;
            HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinje3d_normal");
            MonodomainProblem<3,3> monodomain_problem(&cell_factory_for_just_monodomain);
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();
            monodomain_problem.SetWriteInfo();
            // Solve
            monodomain_problem.Solve();
            soln_mono_repl.ReplicatePetscVector(monodomain_problem.GetSolution());

            HeartEventHandler::Headings();
            HeartEventHandler::Report();
        }

        // Set a radius for each purkinje element
        for (MixedDimensionMesh<3,3>::CableElementIterator iter = mesh.GetCableElementIteratorBegin();
             iter != mesh.GetCableElementIteratorEnd();
             ++iter)
        {
            (*iter)->SetAttribute(0.01); // 50 microns
        }

        {
            HeartEventHandler::Reset();

            // Set up Purkinje problem
            HeartConfig::Instance()->SetOutputDirectory("TestMonodomainPurkinje3d_purkinje");
            PurkinjeCellFactory3d cell_factory;
            MonodomainPurkinjeProblem<3,3> purkinje_problem(&cell_factory);
            purkinje_problem.SetMesh(&mesh);
            purkinje_problem.Initialise();
            purkinje_problem.SetWriteInfo();

            // Solve
            purkinje_problem.Solve();
            soln_repl.ReplicatePetscVector(purkinje_problem.GetSolution());

            HeartEventHandler::Headings();
            HeartEventHandler::Report();
        }

        // Compare solutions at final time - note that the range is [-78.122, 25.5206] for
        // the myocardium voltage, so quite a spread.
        //
        // Note: checked (by looking at SetWriteInfo() output that results are the same at the
        // other times
        for (AbstractTetrahedralMesh<3,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
             node_iter != mesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            unsigned i = node_iter->GetIndex();

            // Purkinje cable runs along x=0,y=0
            if (   (fabs(node_iter->rGetLocation()[0]) < 1e-6)
                && (fabs(node_iter->rGetLocation()[1]) < 1e-6) )
            {
                // check the purkinje voltage is not exactly zero. Although we would expect
                // the purkinje voltage to match the myocardium voltage at these nodes, this
                // isn't the case - but have looked at this and decided it can be put down to mesh
                // error (incl. the fact that although we are stimulating a surface, since
                // p/w linear interpolation is used, the stimulus is interpolated a bit into the
                // interior).
                TS_ASSERT(fabs(soln_repl[2*i+1])>1e-6);

                // hardcoded test for purkinje cell at opposite face (coordinates (0,0,0.1))
                if( fabs(node_iter->rGetLocation()[2]-0.1)<1e-6 )
                {
                    TS_ASSERT_DELTA(soln_repl[2*i+1], -75.9912, 1e-4);
                }
            }
            else
            {
                TS_ASSERT_DELTA(soln_repl[2*i+1], 0.0, 1e-4);
            }

            TS_ASSERT_DELTA(soln_repl[2*i], soln_mono_repl[i], 1e-5);
        }
    }
};

#endif // TESTMONODOMAINPURKINJEPROBLEMLONG_HPP_
