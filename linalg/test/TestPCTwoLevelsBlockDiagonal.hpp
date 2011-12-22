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

#ifndef TESTPCTWOLEVELSBLOCKDIAGONAL_HPP_
#define TESTPCTWOLEVELSBLOCKDIAGONAL_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearSystem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "Timer.hpp"
#include "DistributedVectorFactory.hpp"
#include <cstring>

/*
 * Warning: these tests do not inform PETSc about the nullspace of the matrix. Therefore, convergence might be
 * different compared to a real cardiac simulation. Do not take conclusions about preconditioner performance
 * based on these tests only.
 */
class TestPCTwoLevelsBlockDiagonal : public CxxTest::TestSuite
{
public:

    void TestBasicFunctionality() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        unsigned num_nodes = 441;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        const unsigned num_bath_nodes = 84;
        PetscInt bath_nodes[num_bath_nodes] = {0, 1, 19, 20, 21, 22, 40, 41, 42, 43, 61, 62, 63, 64, 82, 83, 84,
                                               85, 103, 104, 105, 106, 124, 125, 126, 127, 145, 146, 147, 148, 166,
                                               167, 168, 169, 187, 188, 189, 190, 208, 209, 210, 211, 229, 230, 231,
                                               232, 250, 251, 252, 253, 271, 272, 273, 274, 292, 293, 294, 295, 313,
                                               314, 315, 316, 334, 335, 336, 337, 355, 356, 357, 358, 376, 377, 378,
                                               379, 397, 398, 399, 400, 418, 419, 420, 421, 439, 440};
        boost::shared_ptr<std::vector<PetscInt> > p_bath(new std::vector<PetscInt>(bath_nodes, &bath_nodes[num_bath_nodes]));
        assert(p_bath->size() == num_bath_nodes);

        Mat system_matrix;
        PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/PP_system_with_bath.mat", parallel_layout);

        VecDestroy(parallel_layout);

        // Set rhs = A * [1 0 1 0 ... 1 0]'
        Vec one_zeros = factory.CreateVec(2);
        Vec rhs = factory.CreateVec(2);

        for (unsigned node_index=0; node_index<2*num_nodes; node_index+=2)
        {
            PetscVecTools::SetElement(one_zeros, node_index, 1);
            PetscVecTools::SetElement(one_zeros, node_index+1, 0);
        }
        PetscVecTools::Finalise(one_zeros);

        MatMult(system_matrix, one_zeros, rhs);
        VecDestroy(one_zeros);

        LinearSystem ls = LinearSystem(rhs, system_matrix);

        ls.SetAbsoluteTolerance(1e-9);
        ls.SetKspType("cg");
        ls.SetPcType("twolevelsblockdiagonal", p_bath);

        ls.AssembleFinalLinearSystem();

        Vec solution = ls.Solve();

        DistributedVector distributed_solution = factory.CreateDistributedVector(solution);
        DistributedVector::Stripe phi_i(distributed_solution, 0);
        DistributedVector::Stripe phi_e(distributed_solution, 1);

        // Create a std::set of bath node indices for convenience
        std::set<PetscInt> bath_nodes_set(p_bath->begin(), p_bath->end());

        double phi_e_at_tissue = DBL_MAX;

        for (DistributedVector::Iterator index = distributed_solution.Begin();
             index!= distributed_solution.End();
             ++index)
        {
            /*
             * Although we're trying to enforce the solution to be [1 0 ... 1 0], the system is singular and
             * therefore it has infinite solutions. I've (migb) found that the use of different preconditioners
             * lead to different solutions ([0.8 -0.2 ... 0.8 -0.2], [0.5 -0.5 ... 0.5 -0.5], ...)
             *
             * If we were using PETSc null space, it would find the solution that satisfies x'*v=0,
             * being v the null space of the system (v=[1 1 ... 1])
             */

            if (bath_nodes_set.find(index.Global) != bath_nodes_set.end())
            {
                // Bath node: phi_i in the bath is 1 because of the dummy equations we introdude x_i = b_i
                TS_ASSERT_DELTA(phi_i[index], 1.0, 1e-6);

                if (phi_e_at_tissue != DBL_MAX)
                {
                    // phi_e is constant accross the whole domain
                    TS_ASSERT_DELTA(phi_e[index], phi_e_at_tissue, 1e-6);
                }
            }
            else
            {
                // Tissue node
                TS_ASSERT_DELTA(phi_i[index] - phi_e[index], 1.0, 1e-6);

                if (phi_e_at_tissue == DBL_MAX)
                {
                    phi_e_at_tissue = phi_e[index];
                }

            }
        }

        // Coverage (setting PC type after first solve)
        ls.SetPcType("twolevelsblockdiagonal", p_bath);

        MatDestroy(system_matrix);
        VecDestroy(rhs);
        VecDestroy(solution);

#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        const PCType pc;
#else
        PCType pc;
#endif
        PC prec;
        KSPGetPC(ls.mKspSolver, &prec);
        PCGetType(prec, &pc);
        // Although we call it "twolevelsblockdiagonal", PETSc considers this PC a generic SHELL preconditioner
        TS_ASSERT( strcmp(pc,"shell")==0 );

    }

    void TestBetterThanNoPreconditioning()
    {
        EXIT_IF_PARALLEL;

        unsigned num_nodes = 441;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        const unsigned num_bath_nodes = 84;
        PetscInt bath_nodes[num_bath_nodes] = {0, 1, 19, 20, 21, 22, 40, 41, 42, 43, 61, 62, 63, 64, 82, 83, 84,
                                               85, 103, 104, 105, 106, 124, 125, 126, 127, 145, 146, 147, 148, 166,
                                               167, 168, 169, 187, 188, 189, 190, 208, 209, 210, 211, 229, 230, 231,
                                               232, 250, 251, 252, 253, 271, 272, 273, 274, 292, 293, 294, 295, 313,
                                               314, 315, 316, 334, 335, 336, 337, 355, 356, 357, 358, 376, 377, 378,
                                               379, 397, 398, 399, 400, 418, 419, 420, 421, 439, 440};
        boost::shared_ptr<std::vector<PetscInt> > p_bath(new std::vector<PetscInt>(bath_nodes, &bath_nodes[num_bath_nodes]));
        assert(p_bath->size() == num_bath_nodes);

        unsigned point_jacobi_its;
        unsigned block_diag_its;

        Timer::Reset();
        {
            Mat system_matrix;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/PP_system_with_bath.mat", parallel_layout);

            Vec system_rhs;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/PP_system_with_bath.vec", parallel_layout);

            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetAbsoluteTolerance(1e-9);
            ls.SetKspType("cg");
            ls.SetPcType("none");

            Vec solution = ls.Solve();

            point_jacobi_its = ls.GetNumIterations();

            MatDestroy(system_matrix);
            VecDestroy(system_rhs);
            VecDestroy(solution);
        }
        Timer::PrintAndReset("No preconditioning");

        {
            Mat system_matrix;

            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/PP_system_with_bath.mat", parallel_layout);

            Vec system_rhs;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/PP_system_with_bath.vec", parallel_layout);

            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetAbsoluteTolerance(1e-9);
            ls.SetKspType("cg");
            ls.SetPcType("twolevelsblockdiagonal", p_bath);

            Vec solution = ls.Solve();

            block_diag_its = ls.GetNumIterations();

            // Coverage (setting PC type after using blockdiagonal solve)
            ls.SetPcType("twolevelsblockdiagonal", p_bath);

            MatDestroy(system_matrix);
            VecDestroy(system_rhs);
            VecDestroy(solution);
        }
        Timer::Print("Block diagonal preconditioner");

        std::cout << block_diag_its << " " << point_jacobi_its << std::endl;
        TS_ASSERT_LESS_THAN_EQUALS(block_diag_its, point_jacobi_its);

        VecDestroy(parallel_layout);
    }

    /*
     * Use this test to generate a matrix and vector coming from a PP simulation with bath and to
     * work out the indices of the nodes in the bath.
     *
     *  In order to get it to work you need to comment out the following lines:
     *
     *    delete mpSolver;
     *    mpSolver = NULL;
     *
     * in AbstractCardiacProblem.cpp (469 and 470 at the time of writing this) so that you don't
     * destroy the cardiac solver and the end of Solve(). Otherwise you will get a segfault when
     * trying to get the linear system out of the solver object.
     *
     * Add the following includes as well:
     *
     *   #include "BidomainParaParaProblem.hpp"
     *   #include "LuoRudy1991BackwardEuler.hpp"
     *   #include "SimpleBathProblemSetup.hpp"
     *
     */
//    void TestWriteOutPPWithBathLinearSystem()
//    {
//        char output_folder[] = "LAParaParaBath2d";
//
//        HeartConfig::Instance()->SetSimulationDuration(2.0);  //ms
//        HeartConfig::Instance()->SetOutputDirectory(output_folder);
//        HeartConfig::Instance()->SetOutputFilenamePrefix("parapara2d");
//
//        // need to create a cell factory but don't want any intra stim, so magnitude
//        // of stim is zero.
//        c_vector<double,2> centre;
//        centre(0) = 5; // cm
//        centre(1) = 5; // cm
//        BathCellFactory<2,CellLuoRudy1991FromCellMLBackwardEuler> cell_factory( 0.0, centre);
//
//        BidomainParaParaProblem<2> para_para_problem( &cell_factory, true );
//
//        TetrahedralMesh<2,2> mesh;
//        mesh.ConstructRegularSlabMesh(0.05, 1, 1);
//
//        std::cout << "num_nodes: " << mesh.GetNumNodes() << std::endl;
//
//        // set the x<0.15 and x>0.85 regions as the bath region
//        for (unsigned i=0; i<mesh.GetNumElements(); i++)
//        {
//            double x = mesh.GetElement(i)->CalculateCentroid()[0];
//            if ( (x<0.1) || (x>0.9) )
//            {
//                mesh.GetElement(i)->SetRegion(HeartRegionCode::GetValidBathId());
//            }
//        }
//
//        para_para_problem.SetMesh(&mesh);
//
//        double boundary_flux = 20.0e3;
//        double start_time = 0.1;
//        double duration = 0.5; // of the stimulus, in ms
//
//        HeartConfig::Instance()->SetElectrodeParameters(false,0,boundary_flux, start_time, duration);
//
//        para_para_problem.Initialise();
//
//        // Nodes are not labeled until Initialise() is called.
//        std::cout << "bath nodes: ";
//        unsigned num_bath_nodes = 0;
//        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
//        {
//            if (HeartRegionCode::IsRegionBath( mesh.GetNode(i)->GetRegion() ))
//            {
//                std::cout << i << ", ";
//                num_bath_nodes++;
//            }
//        }
//        std::cout << std::endl;
//        std::cout << "num_bath_nodes: " << num_bath_nodes << std::endl;
//
//
//        para_para_problem.Solve();
//
//        OutputFileHandler handler(output_folder, false);
//        PetscTools::DumpPetscObject(para_para_problem.mpSolver->GetLinearSystem()->GetLhsMatrix(), handler.GetOutputDirectoryFullPath() + "PP_system_with_bath.mat");
//        PetscTools::DumpPetscObject(para_para_problem.mpSolver->GetLinearSystem()->GetRhsVector(), handler.GetOutputDirectoryFullPath() + "PP_system_with_bath.vec");
//    }
};

#endif /*TESTPCBLOCKDIAGONAL_HPP_*/
