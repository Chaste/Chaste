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

#ifndef TESTBIDOMAINDISTRIBUTEDMESH_HPP_
#define TESTBIDOMAINDISTRIBUTEDMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "LuoRudy1991.hpp"
#include "BidomainProblem.hpp"
#include "DistributedVector.hpp"
#include "HeartConfig.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "MemfemMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestBidomainDistributedMesh : public CxxTest::TestSuite
{
public:

    void TestBidomainProblemWithDistributedMesh2D()
    {
        HeartConfig::Instance()->SetSimulationDuration(1);  //ms
        HeartConfig::Instance()->SetOutputDirectory("DistributedMesh2d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("tetrahedral2d");

        // The default stimulus in PlaneStimulusCellFactory is not enough to generate propagation
        // here, increasing it an order of magnitude
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-6000);

        // To avoid an issue with the Event handler only one simulation should be
        // in existance at a time: therefore monodomain simulation is defined in a block
        double seq_ave_voltage;
        {
            ///////////////////////////////////////////////////////////////////
            // TetrahedralMesh
            ///////////////////////////////////////////////////////////////////
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            BidomainProblem<2> nondistributed_problem( &cell_factory );
            nondistributed_problem.SetMesh(&mesh);
            nondistributed_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            nondistributed_problem.Solve();

            DistributedVector dist_nondistributed_voltage = nondistributed_problem.GetSolutionDistributedVector();
            DistributedVector::Stripe nondistributed_voltage(dist_nondistributed_voltage, 0);
            DistributedVector::Stripe nondistributed_potential(dist_nondistributed_voltage, 1);

            double seq_local_ave_voltage = 0.0;

            for (DistributedVector::Iterator index = dist_nondistributed_voltage.Begin();
                 index != dist_nondistributed_voltage.End();
                 ++index)
            {
                if (index.Global==0)
                {
                    TS_ASSERT_LESS_THAN(0, nondistributed_voltage[index]);
                }

                seq_local_ave_voltage += nondistributed_voltage[index];
            }

            MPI_Reduce(&seq_local_ave_voltage, &seq_ave_voltage, 1, MPI_DOUBLE, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);
            seq_ave_voltage /= mesh.GetNumNodes();
        }


        ///////////////////////////////////////////////////////////////////
        // DistributedTetrahedralMesh
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputFilenamePrefix("distributed2d");

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);

        mesh.ConstructFromMeshReader(mesh_reader);

        BidomainProblem<2> distributed_problem( &cell_factory );

        distributed_problem.SetMesh(&mesh);

        distributed_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        distributed_problem.Solve();

        DistributedVector dist_distributed_voltage = distributed_problem.GetSolutionDistributedVector();
        DistributedVector::Stripe distributed_voltage(dist_distributed_voltage, 0);
        DistributedVector::Stripe distributed_potential(dist_distributed_voltage, 1);

        double para_local_ave_voltage = 0.0;

        for (DistributedVector::Iterator index = dist_distributed_voltage.Begin();
             index != dist_distributed_voltage.End();
             ++index)
        {
            if (index.Global==0)
            {
                TS_ASSERT_LESS_THAN(0, distributed_voltage[index]);
            }

            para_local_ave_voltage += distributed_voltage[index];
        }


        double para_ave_voltage;
        MPI_Reduce(&para_local_ave_voltage, &para_ave_voltage, 1, MPI_DOUBLE, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);
        para_ave_voltage /= mesh.GetNumNodes();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        if (PetscTools::AmMaster())
        {
            std::cout << seq_ave_voltage << "  " << para_ave_voltage << std::endl;
            TS_ASSERT_DELTA(seq_ave_voltage, para_ave_voltage, 1.0);
        }
    }

    void TestBidomainProblemWithDistributedMesh2DParMetis()
    {
        HeartConfig::Instance()->SetSimulationDuration(1);  //ms
        HeartConfig::Instance()->SetOutputDirectory("DistributedMesh2d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("tetrahedral2d");

        // The default stimulus in PlaneStimulusCellFactory is not enough to generate propagation
        // here, increasing it an order of magnitude
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-6000, 0.5);

        // To avoid an issue with the Event handler only one simulation should be
        // in existance at a time: therefore monodomain simulation is defined in a block
        double seq_ave_voltage=0.0;
        {
            ///////////////////////////////////////////////////////////////////
            // TetrahedralMesh
            ///////////////////////////////////////////////////////////////////
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            BidomainProblem<2> nondistributed_problem( &cell_factory );
            nondistributed_problem.SetMesh(&mesh);
            nondistributed_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            nondistributed_problem.Solve();

            DistributedVector dist_nondistributed_voltage = nondistributed_problem.GetSolutionDistributedVector();
            DistributedVector::Stripe nondistributed_voltage(dist_nondistributed_voltage, 0);
            DistributedVector::Stripe nondistributed_potential(dist_nondistributed_voltage, 1);

            double seq_local_ave_voltage = 0.0;

            for (DistributedVector::Iterator index = dist_nondistributed_voltage.Begin();
                 index != dist_nondistributed_voltage.End();
                 ++index)
            {
                if (index.Global==0)
                {
                    TS_ASSERT_LESS_THAN(0, nondistributed_voltage[index]);
                }

                seq_local_ave_voltage += nondistributed_voltage[index];
            }

            MPI_Reduce(&seq_local_ave_voltage, &seq_ave_voltage, 1, MPI_DOUBLE, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);
            seq_ave_voltage /= mesh.GetNumNodes();
        }

        ///////////////////////////////////////////////////////////////////
        // DistributedTetrahedralMesh
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputFilenamePrefix("distributed2d");
        HeartConfig::Instance()->SetMeshPartitioning("parmetis");
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");

//        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
//        DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY);
//        mesh.ConstructFromMeshReader(mesh_reader);

        BidomainProblem<2> distributed_problem( &cell_factory );

        //distributed_problem.PrintOutput(false);

//        distributed_problem.SetMesh(&mesh);

        distributed_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        distributed_problem.Solve();

        DistributedVector dist_distributed_voltage = distributed_problem.GetSolutionDistributedVector();
        DistributedVector::Stripe distributed_voltage(dist_distributed_voltage, 0);
        DistributedVector::Stripe distributed_potential(dist_distributed_voltage, 1);

        double para_local_ave_voltage = 0.0;

        for (DistributedVector::Iterator index = dist_distributed_voltage.Begin();
             index != dist_distributed_voltage.End();
             ++index)
        {
            if (index.Global==0)
            {
                TS_ASSERT_LESS_THAN(0, distributed_voltage[index]);
            }

            para_local_ave_voltage += distributed_voltage[index];
        }


        double para_ave_voltage;
        MPI_Reduce(&para_local_ave_voltage, &para_ave_voltage, 1, MPI_DOUBLE, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);
        para_ave_voltage /= distributed_problem.rGetMesh().GetNumNodes();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        if (PetscTools::AmMaster())
        {
            std::cout << seq_ave_voltage << "  " << para_ave_voltage << std::endl;
            TS_ASSERT_DELTA(seq_ave_voltage, para_ave_voltage, 1.0);
        }
    }
    void TestBidomainProblemWithDistributedMeshFromMemfem3DParMetis()
    {
        HeartConfig::Instance()->SetSimulationDuration(1);  //ms
        HeartConfig::Instance()->SetOutputDirectory("DistributedMesh3dRepViaTri");
        HeartConfig::Instance()->SetOutputFilenamePrefix("tetrahedral3d");

        // The default stimulus in PlaneStimulusCellFactory is not enough to generate propagation
        // here, increasing it an order of magnitude
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-6000);

        // To avoid an issue with the Event handler only one simulation should be
        // in existance at a time: therefore monodomain simulation is defined in a block
        double seq_ave_voltage=0.0;
        {
            ///////////////////////////////////////////////////////////////////
            // TetrahedralMesh from Triangles
            ///////////////////////////////////////////////////////////////////
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/SlabFromMemfem");
            TetrahedralMesh<3,3> mesh;

            mesh.ConstructFromMeshReader(mesh_reader);
            TS_ASSERT_DELTA(mesh.GetVolume(), 1.5625, 1e-6);
            TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 381u);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 1030u);
            TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 758u);


            BidomainProblem<3> nondistributed_problem( &cell_factory );
            nondistributed_problem.SetMesh(&mesh);
            nondistributed_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            nondistributed_problem.Solve();

            DistributedVector dist_nondistributed_voltage = nondistributed_problem.GetSolutionDistributedVector();
            DistributedVector::Stripe nondistributed_voltage(dist_nondistributed_voltage, 0);
            DistributedVector::Stripe nondistributed_potential(dist_nondistributed_voltage, 1);

            double seq_local_ave_voltage = 0.0;

            for (DistributedVector::Iterator index = dist_nondistributed_voltage.Begin();
                 index != dist_nondistributed_voltage.End();
                 ++index)
            {
                if (index.Global==0)
                {
                    TS_ASSERT_LESS_THAN(0, nondistributed_voltage[index]);
                }

                seq_local_ave_voltage += nondistributed_voltage[index];
            }

            MPI_Reduce(&seq_local_ave_voltage, &seq_ave_voltage, 1, MPI_DOUBLE, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);
            seq_ave_voltage /= mesh.GetNumNodes();
        }

        ///////////////////////////////////////////////////////////////////
        // DistributedTetrahedralMesh from Memfem
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("DistributedMesh3dDistViaMem");

        BidomainProblem<3> distributed_problem( &cell_factory );

        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/Memfem_slab");
        distributed_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        distributed_problem.Solve();

        DistributedVector dist_distributed_voltage = distributed_problem.GetSolutionDistributedVector();
        DistributedVector::Stripe distributed_voltage(dist_distributed_voltage, 0);
        DistributedVector::Stripe distributed_potential(dist_distributed_voltage, 1);

        double para_local_ave_voltage = 0.0;

        for (DistributedVector::Iterator index = dist_distributed_voltage.Begin();
             index != dist_distributed_voltage.End();
             ++index)
        {
            if (index.Global==0)
            {
                TS_ASSERT_LESS_THAN(0, distributed_voltage[index]);
            }

            para_local_ave_voltage += distributed_voltage[index];
        }


        double para_ave_voltage;
        MPI_Reduce(&para_local_ave_voltage, &para_ave_voltage, 1, MPI_DOUBLE, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);
        para_ave_voltage /= distributed_problem.rGetMesh().GetNumNodes();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        if (PetscTools::AmMaster())
        {
            std::cout << seq_ave_voltage << "  " << para_ave_voltage << std::endl;
            TS_ASSERT_DELTA(seq_ave_voltage, para_ave_voltage, 1.0);
        }
    }
};

#endif /*TESTBIDOMAINDISTRIBUTEDMESH_HPP_*/
