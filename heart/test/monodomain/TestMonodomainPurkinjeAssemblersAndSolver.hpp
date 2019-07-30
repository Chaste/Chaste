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

#ifndef TESTMONODOMAINPURKINJEASSEMBLERSANDSOLVER_HPP_
#define TESTMONODOMAINPURKINJEASSEMBLERSANDSOLVER_HPP_

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

typedef MonodomainPurkinjeSolver<2,2> MonodomainPurkinjeSolver2d;
typedef MonodomainPurkinjeCableAssembler<2,2> MonodomainPurkinjeCableAssembler2d;

class NonPurkinjeCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    NonPurkinjeCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-1e6, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        ChastePoint<2> location = pNode->GetPoint();

        if (fabs(location[0])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};


class PurkinjeCellFactory : public AbstractPurkinjeCellFactory<2> //inherits from different base
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PurkinjeCellFactory()
        : AbstractPurkinjeCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-1e6, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        ChastePoint<2> location = pNode->GetPoint();

        if (fabs(location[0])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    AbstractCardiacCell* CreatePurkinjeCellForTissueNode(Node<2>* pNode,
                                                         AbstractCardiacCellInterface* pCardiacCell)
    {
        ChastePoint<2> location = pNode->GetPoint();

        if (fabs(location[0])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};



class TestMonodomainPurkinjeAssemblersAndSolvers : public CxxTest::TestSuite
{
public:
    void TestMonodomainPurkinjeVolumeAssembler()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
        /* The assembly requires 1/time-step,
         * here we are providing enough information without starting a whole simulation */
        PdeSimulationTime::SetTime(0.0);
        PdeSimulationTime::SetPdeTimeStepAndNextTime(0.01, 0.01);

        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.05, 0.1, 0.1);
        unsigned connectivity = mesh.CalculateMaximumNodeConnectivityPerProcess();
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(connectivity, 9u); //Not counting Purkinje interactions
        }

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;
        cell_factory.SetMesh(&mesh);
        MonodomainTissue<2> tissue(&cell_factory);

        // Make sure that a 2Nx2N matrix is partitioned in the same place as an NxN matrix.
        unsigned num_local_nodes = mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        Mat purkinje_mat;
        //Row 8 has 18 non-zeros
        PetscTools::SetupMat(purkinje_mat, 2*mesh.GetNumNodes(), 2*mesh.GetNumNodes(), 2*connectivity, 2*num_local_nodes, 2*num_local_nodes);
        Mat normal_mat;
        PetscTools::SetupMat(normal_mat, mesh.GetNumNodes(), mesh.GetNumNodes(), connectivity, num_local_nodes, num_local_nodes);

        MonodomainPurkinjeVolumeAssembler<2,2> purkinje_vol_assembler(&mesh, &tissue);
        purkinje_vol_assembler.SetMatrixToAssemble(purkinje_mat, true);
        purkinje_vol_assembler.Assemble();

        MonodomainAssembler<2,2> normal_vol_assembler(&mesh, &tissue);
        normal_vol_assembler.SetMatrixToAssemble(normal_mat, true);
        normal_vol_assembler.Assemble();

        PetscMatTools::Finalise(purkinje_mat);
        PetscMatTools::Finalise(normal_mat);

        PetscInt lo, hi;
        PetscMatTools::GetOwnershipRange(purkinje_mat, lo, hi);

        //Check that the partitioning is exactly as expected
        TS_ASSERT_EQUALS((unsigned)lo, 2*mesh.GetDistributedVectorFactory()->GetLow());
        TS_ASSERT_EQUALS((unsigned)hi, 2*mesh.GetDistributedVectorFactory()->GetHigh());

        for (AbstractTetrahedralMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd(); ++node_iter)
        {
            unsigned i = node_iter->GetIndex();
            assert(lo<=(int)(2*i) && (int)(2*i)<hi);
            for(unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i,2*j),   PetscMatTools::GetElement(normal_mat,i,j), 1e-8);
                TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i,2*j+1), 0.0, 1e-8);
            }

            assert(lo<=(int)(2*i+1) && (int)(2*i+1)<hi);
            for(unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j),   0.0, 1e-8);
                TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0, 1e-8);
            }
        }

        PetscTools::Destroy(purkinje_mat);
        PetscTools::Destroy(normal_mat);

    }

    void TestMonodomainPurkinjeCableAssembler()
    {
        /* The assembly requires 1/time-step,
         * here we are providing enough information without starting a whole simulation */
        PdeSimulationTime::SetTime(0.0);
        PdeSimulationTime::SetPdeTimeStepAndNextTime(0.01, 0.01);

        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);

        // There are named indices in the test, so we need predictable numbering...
        MixedDimensionMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(reader);

        Mat purkinje_mat;
        unsigned num_local_nodes = mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        PetscTools::SetupMat(purkinje_mat, 2*mesh.GetNumNodes(), 2*mesh.GetNumNodes(), 9, 2*num_local_nodes, 2*num_local_nodes);

        MonodomainPurkinjeCableAssembler<2,2> purkinje_cable_assembler(&mesh);

        purkinje_cable_assembler.SetMatrixToAssemble(purkinje_mat, true);
        purkinje_cable_assembler.Assemble();
        PetscMatTools::Finalise(purkinje_mat);

        PetscInt lo, hi;
        PetscMatTools::GetOwnershipRange(purkinje_mat, lo, hi);
        for (AbstractTetrahedralMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
             node_iter != mesh.GetNodeIteratorEnd(); ++node_iter)
        {
            unsigned i = node_iter->GetIndex();
            assert(lo<=(int)(2*i) && (int)(2*i)<hi);
            for (unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i,2*j), 0 , 1e-8);
                TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i,2*j+1), 0.0, 1e-8);
            }

            assert(lo<=(int)(2*i+1) && (int)(2*i+1)<hi);
            for (unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                //Non-Purkinje are all zero
                TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j),   0.0, 1e-8);

                //Make sure that columns associated with cable node have non-zero Purkinje entries
                if ( (i>55) && (i<65) && (j>=i-1) && (j<=i+1))
                {
                    TS_ASSERT_DIFFERS( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0);
                }
                else if ((i==55) && (j>=55) && (j<=56) )
                {
                    TS_ASSERT_DIFFERS( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0);
                }
                else if ((i==65) && (j>=64) && (j<=65) )
                {
                    TS_ASSERT_DIFFERS( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0);
                }
                else
                {
                    //Other entries are zero
                    TS_ASSERT_DELTA(PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0 ,1.0e-8);
                }
            }
        }

        PetscTools::Destroy(purkinje_mat);

        // coverage
        MixedDimensionMesh<2,2>::CableElementIterator iter = mesh.GetCableElementIteratorBegin();
        if (iter != mesh.GetCableElementIteratorEnd())
        {
            (*iter)->SetAttribute(0.0);
            TS_ASSERT_THROWS_CONTAINS(MonodomainPurkinjeCableAssembler2d another_assembler(&mesh), "Radii not provided for all Purkinje elements");
        }
    }

    // Solve a Purkinje problem and check that the solution for the myocardium nodes is exactly the same as
    // the solution of an equivalent monodomain-only problem, that the purkinje voltage for
    // purkinje nodes doesn't change (no purkinje stimulus is given), and is 0 for non-purkinje nodes.
    void TestMonodomainPurkinjeSolver()
    {
        /* The assembly requires 1/time-step,
         * here we are providing enough information without starting a whole simulation */
        PdeSimulationTime::SetTime(0.0);
        PdeSimulationTime::SetPdeTimeStepAndNextTime(0.01, 0.01);

        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-12);
        HeartConfig::Instance()->SetPurkinjeSurfaceAreaToVolumeRatio(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());

        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(reader);

        //The mixed dimension mesh specifies fibre radii, we override these for the purposes of the
        //test
        for (MixedDimensionMesh<2, 2>::CableElementIterator iter = mesh.GetCableElementIteratorBegin();
             iter != mesh.GetCableElementIteratorEnd();
             ++iter)
        {
            (*iter)->SetAttribute(0.56419); //1/sqrt(pi), so that the fibre cross section area is 1.0
        }


        std::string mesh_base2("mesh/test/data/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader2(mesh_base2);
        TetrahedralMesh<2,2> mesh_just_monodomain;
        mesh_just_monodomain.ConstructFromMeshReader(reader2);

        PurkinjeCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        NonPurkinjeCellFactory cell_factory_for_just_monodomain;
        cell_factory_for_just_monodomain.SetMesh(&mesh_just_monodomain);

        MonodomainTissue<2> tissue( &cell_factory );
        MonodomainTissue<2> tissue_for_just_monodomain( &cell_factory_for_just_monodomain );

        // Create an empty BCC - zero Neumann BCs will be applied everywhere
        BoundaryConditionsContainer<2,2,2> bcc;
        BoundaryConditionsContainer<2,2,1> bcc_for_just_monodomain;

        MonodomainPurkinjeSolver<2,2> solver(&mesh, &tissue, &bcc);

        MonodomainSolver<2,2> solver_just_monodomain(&mesh_just_monodomain, &tissue_for_just_monodomain, &bcc_for_just_monodomain);


        //Create an initial condition
        DistributedVectorFactory* p_vector_factory = mesh.GetDistributedVectorFactory();
        Vec init_cond = p_vector_factory->CreateVec(2);
        Vec init_cond_just_monodomain = p_vector_factory->CreateVec(1);

        // get the voltage stripes
        DistributedVector ic = mesh.GetDistributedVectorFactory()->CreateDistributedVector(init_cond);
        DistributedVector ic2 = mesh.GetDistributedVectorFactory()->CreateDistributedVector(init_cond_just_monodomain);
        DistributedVector::Stripe volume_stripe = DistributedVector::Stripe(ic, 0);
        DistributedVector::Stripe cable_stripe = DistributedVector::Stripe(ic, 1);

        for (DistributedVector::Iterator index = ic.Begin();
             index != ic.End();
             ++index)
        {
            volume_stripe[index] = tissue.GetCardiacCell(index.Global)->GetVoltage();
            cable_stripe[index] = tissue.GetPurkinjeCell(index.Global)->GetVoltage();//doesn't matter if this is fake
            // make it zero in the cable stripe for the nodes that are not in purkinje ..
        }
        ic.Restore();

        for (DistributedVector::Iterator index = ic2.Begin();
             index != ic2.End();
             ++index)
        {
            ic2[index] = tissue.GetCardiacCell(index.Global)->GetVoltage();
        }
        ic2.Restore();


        double t_end = 1;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(0.01);
        solver.SetInitialCondition(init_cond);

        solver.SetOutputDirectoryAndPrefix("MonodomainPurkinje","results");
        solver.SetOutputToTxt(true);
        solver.SetPrintingTimestepMultiple(10);

        Vec solution = solver.Solve();

        // the following assumes Luo-Rudy!!
        Vec init_cond_for_just_monodomain = PetscTools::CreateAndSetVec(mesh.GetNumNodes(), -83.853);

        solver_just_monodomain.SetTimes(0, t_end);
        solver_just_monodomain.SetTimeStep(0.01);
        solver_just_monodomain.SetInitialCondition(init_cond_for_just_monodomain);

        Vec solution_just_monodomain = solver_just_monodomain.Solve();


        // Test that certain blocks of the matrix and rhs vector in the monodomain-purkinje solve
        // match the matrix and vector for a normal monodomain solve
        Vec& r_purk_rhs = solver.GetLinearSystem()->rGetRhsVector();
        Vec& r_mono_rhs = solver_just_monodomain.GetLinearSystem()->rGetRhsVector();

        Mat& r_purk_mat = solver.GetLinearSystem()->rGetLhsMatrix();
        Mat& r_mono_mat = solver_just_monodomain.GetLinearSystem()->rGetLhsMatrix();

        TS_ASSERT_EQUALS(PetscVecTools::GetSize(r_purk_rhs), 2*PetscVecTools::GetSize(r_mono_rhs));
        assert(PetscVecTools::GetSize(r_mono_rhs)==mesh.GetNumNodes());

        int lo, hi;
        VecGetOwnershipRange(r_mono_rhs, &lo, &hi);

        // We don't explicitly test the values of the rhs, as this is also tested by comparing the solutions.
        // As the system matrices are different, the results won't be identical (but will be close) for the same linear solver tolerance.
        // for(int i=lo; i<hi; i++)
        // {
        //     TS_ASSERT_DELTA(PetscVecTools::GetElement(r_purk_rhs, 2*i), PetscVecTools::GetElement(r_mono_rhs, i), 1e-8);
        // }

        for(int i=lo; i<hi; i++)
        {
            for(unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                // 'top-left' block
                TS_ASSERT_DELTA(PetscMatTools::GetElement(r_purk_mat, 2*i,2*j), PetscMatTools::GetElement(r_mono_mat, i,j), 1e-8);
                // 'off-diagonal' blocks
                TS_ASSERT_DELTA(PetscMatTools::GetElement(r_purk_mat, 2*i,2*j+1), 0.0, 1e-8);
                TS_ASSERT_DELTA(PetscMatTools::GetElement(r_purk_mat, 2*i+1,2*j), 0.0, 1e-8);
            }
        }


        ReplicatableVector soln_repl(solution);
        ReplicatableVector soln_mono_repl(solution_just_monodomain);

        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if(55<=i && i<=65) // purkinje nodes for this mesh
            {
                //The Purkinje domain is a 1D line embedded within the tissue.
                //It is stimulated in the same way as the tissue domain, therefore
                //the propagation velocity should be the same (within numerical error).
                TS_ASSERT_DELTA(soln_repl[2*i+1], soln_repl[2*i], 1e-1);
            }
            else
            {
                TS_ASSERT_DELTA(soln_repl[2*i+1], 0.0, 1e-4)
            }

            TS_ASSERT_DELTA(soln_repl[2*i], soln_mono_repl[i], 1e-5);
        }

        HeartConfig::Instance()->SetUseStateVariableInterpolation(true);
        TS_ASSERT_THROWS_THIS(MonodomainPurkinjeSolver2d bad_solver(&mesh, &tissue, &bcc),"State-variable interpolation is not yet supported with Purkinje");
        HeartConfig::Instance()->SetUseStateVariableInterpolation(false);


        PetscTools::Destroy(solution);
        PetscTools::Destroy(init_cond);
        PetscTools::Destroy(init_cond_for_just_monodomain);
        PetscTools::Destroy(solution_just_monodomain);
        PetscTools::Destroy(init_cond_just_monodomain);
    }
};

#endif // TESTMONODOMAINPURKINJEASSEMBLERSANDSOLVER_HPP_
