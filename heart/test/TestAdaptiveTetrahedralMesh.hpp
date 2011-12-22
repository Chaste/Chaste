/*

Copyright (C) Fujitsu Laboratories of Europe, 2009-2010

*/

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



#ifndef TESTADAPTIVETETRAHEDRALMESH_HPP_
#define TESTADAPTIVETETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <fstream>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#ifdef CHASTE_ADAPTIVITY
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkTriangle.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellDerivatives.h>
#include <vtkDataSetToUnstructuredGridFilter.h>
#include <vtkCellDataToPointData.h>
#include <vtkExtractVectorComponents.h>
#include <vtkPointData.h>
#endif //CHASTE_ADAPTIVITY

#include "AdaptiveTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "VtkMeshReader.hpp"
#include "TetrahedralMesh.hpp"
//#include "DistributedTetrahedralMesh.hpp"
#include "HeartConfig.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestAdaptiveTetrahedralMesh : public CxxTest::TestSuite
{
public:

    void tearDown()
    {
        HeartConfig::Reset();
    }

    void TestConstructAdaptiveTetrahedralMeshFromVtu(void) throw(Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        AdaptiveTetrahedralMesh *adaptive_mesh = new AdaptiveTetrahedralMesh;
        adaptive_mesh->ConstructFromVtuFile("heart/test/data/adaptivity/twin_flow.vtu");

        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveTetrahedralMesh");

        TS_ASSERT_EQUALS( adaptive_mesh->GetNumNodes(), 48706U );
        TS_ASSERT_EQUALS( adaptive_mesh->GetNumElements(), 241857U );

        delete adaptive_mesh;
#endif
        }

    void TestConstructAdaptiveTetrahedralMeshFromMesh(void) throw(Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        TrianglesMeshReader<3,3> triangles_reader("mesh/test/data/slab_138_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(triangles_reader);

        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveTetrahedralMesh");

        AdaptiveTetrahedralMesh adaptive_mesh;
        adaptive_mesh.ConstructFromMesh(&mesh);

        adaptive_mesh.MakeVerbose(false);   // Coverage - this is a no-op because the default is false.

        TS_ASSERT_EQUALS( adaptive_mesh.GetNumElements(), 138U );
#endif
    }

    void TestConstructAdaptiveTetrahedralMeshFromDistributedMesh(void) throw(Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        TrianglesMeshReader<3,3> triangles_reader("mesh/test/data/slab_138_elements");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(triangles_reader);

        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveTetrahedralMesh");

        AdaptiveTetrahedralMesh adaptive_mesh;
        adaptive_mesh.ConstructFromDistributedMesh(&mesh);

        TS_ASSERT_EQUALS( adaptive_mesh.GetNumNodes(), 59U );
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumElements(), 138U );

        TS_ASSERT_EQUALS( adaptive_mesh.GetNumLocalAndHaloNodes(), (unsigned) adaptive_mesh.GetVtkUnstructuredGrid()->GetNumberOfPoints() )
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumLocalElements(), (unsigned) adaptive_mesh.GetVtkUnstructuredGrid()->GetNumberOfCells() )

        unsigned num_local_nodes = adaptive_mesh.GetNumLocalNodes();
        if (PetscTools::IsParallel())
        {
            TS_ASSERT_LESS_THAN( num_local_nodes, adaptive_mesh.GetNumLocalAndHaloNodes() ); // should be at least one halo node
        }

        unsigned total_nodes;
        MPI_Allreduce( &num_local_nodes, &total_nodes, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD );
        TS_ASSERT_EQUALS( total_nodes, 59u ); // Each node is privately owned by exactly one process.
#endif
    }

    void TestAddPointDataAndWriteToSequentialFile(void) throw(Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        TrianglesMeshReader<3,3> reader("heart/test/data/halfheart");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveTetrahedralMesh");
        OutputFileHandler file_handler(HeartConfig::Instance()->GetOutputDirectory(), false);
        std::string output_directory = file_handler.GetOutputDirectoryFullPath();

        AdaptiveTetrahedralMesh adaptive_mesh;
        adaptive_mesh.ConstructFromMesh(&mesh);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        adaptive_mesh.AddPointData("Distance from origin", distance);

        adaptive_mesh.WriteMeshToFile( output_directory, "halfheart_with_data.vtu" );
#endif
    }

    void TestAddPointDataAndWriteToDistributedFile(void) throw(Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        TrianglesMeshReader<3,3> reader("mesh/test/data/slab_138_elements");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveTetrahedralMesh");
        OutputFileHandler file_handler(HeartConfig::Instance()->GetOutputDirectory(), false);
        std::string output_directory = file_handler.GetOutputDirectoryFullPath();

        AdaptiveTetrahedralMesh adaptive_mesh;
        adaptive_mesh.ConstructFromDistributedMesh(&mesh);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        //Privately owned nodes
        for (DistributedTetrahedralMesh<3,3>::NodeIterator it=mesh.GetNodeIteratorBegin();
             it != mesh.GetNodeIteratorEnd();
             ++it)
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(it->rGetLocation()));
        }
        // Halo nodes
        std::vector<unsigned> halo_nodes;
        mesh.GetHaloNodeIndices(halo_nodes);
        for(unsigned i=0; i<halo_nodes.size(); i++)
        {
            distance.push_back(norm_2(mesh.GetNodeOrHaloNode(halo_nodes[i])->rGetLocation()));
        }

        adaptive_mesh.AddPointData("Distance from origin", distance);

        adaptive_mesh.WriteMeshToDistributedFile( output_directory, "slab_138_elements.pvtu" );
#endif
    }

    void TestConvertVtuToPvtu(void) throw(Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        VtkMeshReader<3,3> reader("heart/test/data/adaptivity/coarse_slab_neumann_no_adapt0002.vtu");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB); // No re-ordering, since this will mess up point data
        mesh.ConstructFromMeshReader(reader);

        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveTetrahedralMesh");
        OutputFileHandler file_handler(HeartConfig::Instance()->GetOutputDirectory(), false);
        std::string output_directory = file_handler.GetOutputDirectoryFullPath();

        AdaptiveTetrahedralMesh adaptive_mesh;
        adaptive_mesh.ConstructFromDistributedMesh(&mesh);

        std::vector<unsigned> halo_nodes;
        mesh.GetHaloNodeIndices(halo_nodes);

        std::vector<double> global_vm;
        reader.GetPointData("Vm", global_vm);
        std::vector<double> local_vm;

        for (DistributedTetrahedralMesh<3,3>::NodeIterator it=mesh.GetNodeIteratorBegin();
             it != mesh.GetNodeIteratorEnd();
             ++it)
        {
            local_vm.push_back(global_vm[it->GetIndex()]);
        }
        for(unsigned i=0; i<halo_nodes.size(); i++)
        {
            local_vm.push_back(global_vm[halo_nodes[i]]);
        }

        TS_ASSERT_EQUALS( local_vm.size(), adaptive_mesh.GetNumLocalAndHaloNodes() );
        adaptive_mesh.AddPointData("Vm", local_vm);

        adaptive_mesh.WriteMeshToDistributedFile( output_directory, "coarse_slab_with_data.pvtu" );
#endif
    }

    void TestGetEdgeLengthDistribution(void) throw(Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        AdaptiveTetrahedralMesh adaptive_mesh;
        adaptive_mesh.ConstructFromVtuFile("heart/test/data/adaptivity/coarse_slab_neumann0002.vtu");

        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveTetrahedralMesh");

        TS_ASSERT_EQUALS( adaptive_mesh.GetNumNodes(), 1077U );
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumElements(), 4724U );

        double error = 1.0;
        double sigma = 0.01;
        double max_length = 0.04;
        double min_length = 0.005;
        double gradation = 1.3;
        int max_nodes = 1000;

        HeartConfig::Instance()->SetAdaptivityParameters( error, sigma, max_length, min_length,
                                                          gradation, max_nodes, 5 );

        adaptive_mesh.GetGeometryConstraints();
        adaptive_mesh.CalculateErrorMetric();
        adaptive_mesh.mpVtkUnstructuredGrid->GetPointData()->RemoveArray("mean_desired_lengths");
        adaptive_mesh.mpVtkUnstructuredGrid->GetPointData()->RemoveArray("desired_lengths");

//        double distribution;
//        double range = 0.0;

        TS_ASSERT_DELTA( adaptive_mesh.GetEdgeLengthDistribution(0.0), 1.0, 1e-3 );
        TS_ASSERT_DELTA( adaptive_mesh.GetEdgeLengthDistribution(1.0), 0.0, 1e-3 );

//        while (range < 1.05)
//        {
//            distribution = adaptive_mesh.GetEdgeLengthDistribution(range);
//            std::cout << "range = " << range << ", proportion of bad edges: " << distribution << std::endl;
//            range = range + 0.1;
//        }

//        adaptive_mesh.SetAdaptCriterion(0.8, 0.45);
//        adaptive_mesh.Adapt();
//
//        TS_ASSERT( (adaptive_mesh.GetAdaptSuccess()) );

#endif
}

    void xTestGetEdgeLengthDistributionWhenSolutionChangesSlowly(void) throw(Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        EXIT_IF_PARALLEL;

        AdaptiveTetrahedralMesh adaptive_mesh;
        adaptive_mesh.ConstructFromVtuFile("/work2/southern/VTK_files/adapting_oxford_heart/adapting_oxford_heart0400.vtu");
//        adaptive_mesh.ConstructFromVtuFile("testoutput/TestAdaptiveBidomainProblem/coarse_slab_neumann0035.vtu");

        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveTetrahedralMesh");

        TS_ASSERT_EQUALS( adaptive_mesh.GetNumNodes(), 65012U );
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumElements(), 283405U );

        double error = 2.0;
        double sigma = 0.01;
        double max_length = 0.5;
        double min_length = 0.01;    // mean edge length in full res Oxford heart mesh is 0.0125 (cm)
        double gradation = 1.5;
        int max_nodes = 1e6;        // number of nodes in full res Oxford heart mesh is ~4e6

        HeartConfig::Instance()->SetAdaptivityParameters( error, sigma, max_length, min_length,
                                                          gradation, max_nodes, 5 );

        adaptive_mesh.CalculateSENListAndSids();
        adaptive_mesh.GetGeometryConstraints();
        adaptive_mesh.CalculateErrorMetric();
        adaptive_mesh.mpVtkUnstructuredGrid->GetPointData()->RemoveArray("mean_desired_lengths");
        adaptive_mesh.mpVtkUnstructuredGrid->GetPointData()->RemoveArray("desired_lengths");

        double distribution;
        double range = 0.0;

        while (range < 1.05)
        {
            distribution = adaptive_mesh.GetEdgeLengthDistribution(range);
            std::cout << "range = " << range << ", proportion of bad edges: " << distribution << std::endl;
            range = range + 0.1;
        }

//        adaptive_mesh.SetAdaptCriterion(0.8, 0.45);
//        adaptive_mesh.Adapt();
//
//        std::cout << "Adapted mesh has " <<  adaptive_mesh.GetNumNodes() << " nodes and " <<
//                     adaptive_mesh.GetNumElements() << " elements." << std::endl;
//
//        TS_ASSERT( !(adaptive_mesh.GetAdaptSuccess()) );

#endif
    }

    void TestAdaptingTheCoarseSlabMesh(void) throw(Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        EXIT_IF_PARALLEL;

        AdaptiveTetrahedralMesh adaptive_mesh;
        adaptive_mesh.ConstructFromVtuFile("heart/test/data/adaptivity/small_bidomain_slab.vtu");

        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveTetrahedralMesh");
        OutputFileHandler file_handler(HeartConfig::Instance()->GetOutputDirectory(), false);
        std::string output_directory = file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS( adaptive_mesh.GetNumNodes(), 1331U );
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumElements(), 6000U );

        double error = 100.0;
        double sigma = 0.01;
        double max_length = 0.04;
        double min_length = 0.0025;
        double gradation = 1.3;
        int max_nodes = 1000;
        int num_adaptive_sweeps = 5;

        HeartConfig::Instance()->SetAdaptivityParameters( error, sigma, max_length, min_length,
                                                          gradation, max_nodes, num_adaptive_sweeps );

        // Test no adapt happens when adaptivity criterion set so that any edge is "good".
        adaptive_mesh.SetAdaptCriterion( 1.0, 1.0 );
        adaptive_mesh.CalculateSENListAndSids();
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumSurfaceElements(), 1200U );
        adaptive_mesh.GetGeometryConstraints();
        adaptive_mesh.CalculateErrorMetric();
        adaptive_mesh.mpVtkUnstructuredGrid->GetPointData()->RemoveArray("mean_desired_lengths");
        adaptive_mesh.mpVtkUnstructuredGrid->GetPointData()->RemoveArray("desired_lengths");
        adaptive_mesh.SetAdaptCriterion( 1.0, 1.0 );
        adaptive_mesh.Adapt();

        TS_ASSERT_EQUALS( adaptive_mesh.GetNumNodes(), 1331U );
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumElements(), 6000U );
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumSurfaceElements(), 1200U );
        TS_ASSERT( !(adaptive_mesh.GetAdaptSuccess()) );

        // Test that we do get adaptivity when adaptivity criterion set back to default values
        adaptive_mesh.SetAdaptCriterion( 0.0, 0.0 );
        adaptive_mesh.GetGeometryConstraints();
//        adaptive_mesh.mpDiscreteGeometryConstraints->write_vtk(std::string("sids.vtu"));
        adaptive_mesh.CalculateErrorMetric();
        adaptive_mesh.WriteMeshToFile(output_directory, std::string("metric.vtu"));
        adaptive_mesh.mpVtkUnstructuredGrid->GetPointData()->RemoveArray("mean_desired_lengths");
        adaptive_mesh.mpVtkUnstructuredGrid->GetPointData()->RemoveArray("desired_lengths");
        adaptive_mesh.Adapt();

        adaptive_mesh.WriteMeshToFile(output_directory, std::string("adapted.vtu"));

        TS_ASSERT_DELTA( adaptive_mesh.GetNumNodes(), 700U, 50 );
        TS_ASSERT_DELTA( adaptive_mesh.GetNumElements(), 3000U, 150 );
        TS_ASSERT_DELTA( adaptive_mesh.GetNumSurfaceElements(), 662U, 50 );
        TS_ASSERT( adaptive_mesh.GetAdaptSuccess() );

#endif
    }

    void TestAdaptingTheCoarseSlabMeshInOneStep(void) throw(Exception)
    {
#ifdef CHASTE_ADAPTIVITY
        EXIT_IF_PARALLEL;

        AdaptiveTetrahedralMesh adaptive_mesh;
        adaptive_mesh.ConstructFromVtuFile("heart/test/data/adaptivity/small_bidomain_slab.vtu");

        HeartConfig::Instance()->SetOutputDirectory("TestAdaptiveTetrahedralMesh");

        TS_ASSERT_EQUALS( adaptive_mesh.GetNumNodes(), 1331U );
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumElements(), 6000U );

        double error = 100.0;
        double sigma = 0.01;
        double max_length = 0.04;
        double min_length = 0.0025;
        double gradation = 1.3;
        int max_nodes = 1000;
        int num_adaptive_sweeps = 5;

        HeartConfig::Instance()->SetAdaptivityParameters( error, sigma, max_length, min_length,
                                                          gradation, max_nodes, num_adaptive_sweeps );

        // Test no adapt happens when adaptivity criterion set so that any edge is "good".
        adaptive_mesh.SetAdaptCriterion( 1.0, 1.0 );
        adaptive_mesh.CalculateSENListAndSids();
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumSurfaceElements(), 1200U );
        adaptive_mesh.SetAdaptCriterion( 1.0, 1.0 );
        adaptive_mesh.AdaptMesh();

        TS_ASSERT_EQUALS( adaptive_mesh.GetNumNodes(), 1331U );
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumElements(), 6000U );
        TS_ASSERT_EQUALS( adaptive_mesh.GetNumSurfaceElements(), 1200U );
        TS_ASSERT( !(adaptive_mesh.GetAdaptSuccess()) );

        // Test that we do get adaptivity when adaptivity criterion set back to default values
        adaptive_mesh.SetAdaptCriterion( 0.0, 0.0 );
        adaptive_mesh.AdaptMesh();

        TS_ASSERT_DELTA( adaptive_mesh.GetNumNodes(), 700U, 50U );
        TS_ASSERT_DELTA( adaptive_mesh.GetNumElements(), 3000U, 150U );
        TS_ASSERT_DELTA( adaptive_mesh.GetNumSurfaceElements(), 662U, 50U );
        TS_ASSERT( adaptive_mesh.GetAdaptSuccess() );
#endif
    }
};

#endif /*TESTADAPTIVETETRAHEDRALMESH_HPP_*/
