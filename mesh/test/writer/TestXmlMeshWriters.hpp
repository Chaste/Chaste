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


#ifndef _TESTXMLMESHWRITERS_HPP_
#define _TESTXMLMESHWRITERS_HPP_

#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "TetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "XdmfMeshWriter.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "MixedDimensionMesh.hpp"
#include "QuadraticMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "FileComparison.hpp"
#include <iostream>

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

class TestXmlMeshWriters : public CxxTest::TestSuite
{
public:
    void TestBasicVtkMeshWriter()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "cube_2mm_12_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        //1.6K uncompressed, 1.3K compressed
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/";

        {
            //Check that the reader can see it
            VtkMeshReader<3,3> vtk_reader(results_dir+"cube_2mm_12_elements.vtu");
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumElements());
            TS_ASSERT_EQUALS(vtk_reader.GetNumFaces(), mesh.GetNumBoundaryElements());

            // Check we have the right number of nodes & elements when we re-construct it
            TetrahedralMesh<3,3> vtk_mesh;
            vtk_mesh.ConstructFromMeshReader(vtk_reader);
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), vtk_mesh.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumElements(), vtk_mesh.GetNumElements());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), vtk_mesh.GetNumBoundaryElements());
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestSequentialMeshCannotWriteParallelFiles()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "cube_2mm_12_elements");

        TS_ASSERT_THROWS_THIS( writer.SetParallelFiles(mesh),
                               "Cannot write parallel files using a sequential mesh");
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestParallelVtkMeshWriter()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "cube_2mm_12_elements_with_data");

        writer.SetParallelFiles(mesh);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (DistributedTetrahedralMesh<3,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
               node_iter != mesh.GetNodeIteratorEnd();
               ++node_iter)
        {
            distance.push_back(norm_2(node_iter->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        // Add location (vector) to "point" data
        std::vector< c_vector<double, 3> > location;
        for (DistributedTetrahedralMesh<3,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
               node_iter != mesh.GetNodeIteratorEnd();
               ++node_iter)
        {
            location.push_back(node_iter->rGetLocation());
        }
        writer.AddPointData("Location", location);

        // Add element quality into the element "cell" data
        std::vector<double> quality;
        for (DistributedTetrahedralMesh<3,3>::ElementIterator ele_iter = mesh.GetElementIteratorBegin();
               ele_iter != mesh.GetElementIteratorEnd();
               ++ele_iter)
        {
            quality.push_back(ele_iter->CalculateQuality());
        }
        writer.AddCellData("Quality", quality);

        // Add fibre type to "cell" data
        std::vector< c_vector<double, 3> > centroid;
        for (DistributedTetrahedralMesh<3,3>::ElementIterator ele_iter = mesh.GetElementIteratorBegin();
               ele_iter != mesh.GetElementIteratorEnd();
               ++ele_iter)
        {
            centroid.push_back(ele_iter->CalculateCentroid());
        }
        writer.AddCellData("Centroid", centroid);

        writer.WriteFilesUsingMesh(mesh);

        PetscTools::Barrier("Wait for files to be written");

        std::stringstream filepath;
        filepath << OutputFileHandler::GetChasteTestOutputDirectory() << "TestVtkMeshWriter/cube_2mm_12_elements_with_data";
        // Add suffix to VTK vtu file
        if (PetscTools::IsSequential())
        {
            filepath <<  ".vtu";
        }
        else
        {
            /*
             * Check that the pvtu file exists. Note that checking its content is hard
             * because the number of processes (.vtu file references) will vary.
             */
            FileFinder vtk_file(filepath.str() + ".pvtu", RelativeTo::Absolute);
            TS_ASSERT(vtk_file.Exists());
            // Add suffix to VTK vtu file
            filepath << "_" << PetscTools::GetMyRank() << ".vtu";
        }
        {
            // Check that the reader can see it
            VtkMeshReader<3,3> vtk_reader(filepath.str());
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumLocalNodes() + mesh.GetNumHaloNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumLocalElements());

            // Check that it has the correct data
            std::vector<double> distance_read;
            vtk_reader.GetPointData("Distance from origin", distance_read);
            TS_ASSERT_EQUALS(distance.size(), mesh.GetNumLocalNodes() );
            TS_ASSERT_EQUALS(distance_read.size(), mesh.GetNumLocalNodes()  + mesh.GetNumHaloNodes());

            for (unsigned i=0; i<distance.size(); i++)
            {
                TS_ASSERT_EQUALS(distance[i], distance_read[i]);
            }

            std::vector<c_vector<double,3> > location_read;
            vtk_reader.GetPointData("Location", location_read);
            TS_ASSERT_EQUALS(location.size(), mesh.GetNumLocalNodes() );
            TS_ASSERT_EQUALS(location_read.size(), mesh.GetNumLocalNodes()  + mesh.GetNumHaloNodes());
            for (unsigned i=0; i<location.size(); i++)
            {
                for (unsigned j=0; j<3; j++)
                {
                    TS_ASSERT_EQUALS(location[i][j], location_read[i][j]);
                }
            }

            std::vector<double> quality_read;
            vtk_reader.GetCellData("Quality", quality_read);
            TS_ASSERT_EQUALS(quality.size(), mesh.GetNumLocalElements() );
            TS_ASSERT_EQUALS(quality_read.size(), mesh.GetNumLocalElements());
            for (unsigned i=0; i<quality_read.size(); i++)
            {
                TS_ASSERT_EQUALS(quality[i], quality_read[i]);
            }

            std::vector<c_vector<double,3> > centroid_read;
            vtk_reader.GetCellData("Centroid", centroid_read);
            TS_ASSERT_EQUALS(centroid.size(), mesh.GetNumLocalElements() );
            TS_ASSERT_EQUALS(centroid_read.size(), mesh.GetNumLocalElements());
            for (unsigned i=0; i<centroid_read.size(); i++)
            {
                for (unsigned j=0; j<3; j++)
                {
                    TS_ASSERT_EQUALS(centroid[i][j], centroid_read[i][j]);
                }
            }
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestParallelVtkMeshWriter2d()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_200_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<2,2> writer("TestVtkMeshWriter", "2D_0_to_1mm_200_elements_parallel_data", false);
        writer.SetParallelFiles(mesh);
        // Add distance from origin into the node "point" data
        std::vector<double> rank;
        // Real rank for the owned nodes
        for (unsigned i=0; i<mesh.GetNumLocalNodes(); i++)
        {
            rank.push_back(PetscTools::GetMyRank());
        }
        writer.AddPointData("Process rank", rank);

        writer.WriteFilesUsingMesh(mesh);

        PetscTools::Barrier("Wait for files to be written");

        std::stringstream filepath;
        filepath << OutputFileHandler::GetChasteTestOutputDirectory() << "TestVtkMeshWriter/2D_0_to_1mm_200_elements_parallel_data";
        // Add suffix to VTK vtu file.
        if (PetscTools::IsSequential())
        {
            filepath <<  ".vtu";
        }
        else
        {
            /*
             * Check that the pvtu file exists. Note that checking its content is hard
             * because the number of processes (.vtu file references) will vary.
             */
            FileFinder vtk_file(filepath.str() + ".pvtu", RelativeTo::Absolute);
            TS_ASSERT(vtk_file.Exists());
            // Add suffix to VTK vtu file
            filepath << "_" << PetscTools::GetMyRank() << ".vtu";
        }

        {
            // Check that the reader can see it
            VtkMeshReader<2,2> vtk_reader(filepath.str());
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumLocalNodes() + mesh.GetNumHaloNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumLocalElements());

            // Check that it has the correct data
            std::vector<double> rank_read;
            vtk_reader.GetPointData("Process rank", rank_read);
            TS_ASSERT_EQUALS(rank.size(), mesh.GetNumLocalNodes() );
            TS_ASSERT_EQUALS(rank_read.size(), mesh.GetNumLocalNodes()  + mesh.GetNumHaloNodes());
            for (unsigned i=0; i<rank.size(); i++)
            {
                TS_ASSERT_EQUALS(rank[i], rank_read[i]);
                TS_ASSERT_EQUALS(rank_read[i], PetscTools::GetMyRank());
            }
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestVtkMeshWriter2D()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_200_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<2,2> writer("TestVtkMeshWriter", "2D_0_to_1mm_200_elements", false);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        // Add fibre type to "point" data
        std::vector< c_vector<double, 2> > location;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            location.push_back(mesh.GetNode(i)->rGetLocation());
        }
        writer.AddPointData("Location", location);

        // Add element quality into the element "cell" data
        std::vector<double> quality;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            quality.push_back(mesh.GetElement(i)->CalculateQuality());
        }
        writer.AddCellData("Quality", quality);

        // Add fibre type to "cell" data
        std::vector< c_vector<double, 2> > centroid;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            centroid.push_back(mesh.GetElement(i)->CalculateCentroid());
        }
        writer.AddCellData("Centroid", centroid);

        // Add Jacobian tensor cell data, covering both the symmetric and non-symmetric tensor methods
        std::vector< c_matrix<double, 2, 2> > jacobian;
        std::vector< c_vector<double, 3> > squared_jacobian;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            c_matrix<double, 2, 2> element_jacobian;
            double element_jacobian_determinant;
            mesh.GetElement(i)->CalculateJacobian(element_jacobian, element_jacobian_determinant);
            jacobian.push_back(element_jacobian);

            c_matrix<double, 2, 2> squared_element_jacobian;
            squared_element_jacobian = prod(element_jacobian, element_jacobian);

            c_vector<double, 3> tri_squared_element_jacobian;
            //We store [T00 T01 T02 T11 T12 T22]
            tri_squared_element_jacobian(0) = squared_element_jacobian(0, 0);
            tri_squared_element_jacobian(1) = squared_element_jacobian(0, 1);
            tri_squared_element_jacobian(2) = squared_element_jacobian(1, 1);

            squared_jacobian.push_back(tri_squared_element_jacobian);
        }
        writer.AddTensorCellData("Jacobian", jacobian);
        writer.AddTensorCellData("SquaredJacobian", squared_jacobian);

        //Add tensor point data, we use the outer product of the node's location
        std::vector< c_matrix<double, 2, 2> > location_outer_product;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_matrix<double, 2, 2> element_location_outer_product;
            element_location_outer_product = outer_prod(trans(mesh.GetNode(i)->rGetLocation()), mesh.GetNode(i)->rGetLocation());
            location_outer_product.push_back(element_location_outer_product);
        }
        writer.AddTensorPointData("LocationProduct", location_outer_product);

        writer.WriteFilesUsingMesh(mesh);
        //13K uncompressed, 3.7K compressed

        {
            // Check that the reader can see it
            VtkMeshReader<2,2> vtk_reader(OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/2D_0_to_1mm_200_elements.vtu");
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumElements());

            // Check that it has the correct data
            std::vector<double> distance_read;
            vtk_reader.GetPointData("Distance from origin", distance_read);
            for (unsigned i=0; i<distance_read.size(); i++)
            {
                TS_ASSERT_EQUALS(distance[i], distance_read[i]);
            }
            std::vector<c_vector<double,2> > location_read;
            vtk_reader.GetPointData("Location", location_read);
            for (unsigned i=0; i<location_read.size(); i++)
            {
                for (unsigned j=0; j<2; j++)
                {
                    TS_ASSERT_EQUALS(location[i][j], location_read[i][j]);
                }
            }
            std::vector<double> quality_read;
            vtk_reader.GetCellData("Quality", quality_read);
            for (unsigned i=0; i<quality_read.size(); i++)
            {
                TS_ASSERT_EQUALS(quality[i], quality_read[i]);
            }
            std::vector<c_vector<double,2> > centroid_read;
            vtk_reader.GetCellData("Centroid", centroid_read);
            for (unsigned i=0; i<centroid_read.size(); i++)
            {
                for (unsigned j=0; j<2; j++)
                {
                    TS_ASSERT_EQUALS(centroid[i][j], centroid_read[i][j]);
                }
            }

            ///\todo #2254  Implement reading of tensor data & test.
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestParallelVtkMeshWriter1d()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_10_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<1,1> writer("TestVtkMeshWriter", "1D_0_to_1_10_elements_parallel_data", false);
        writer.SetParallelFiles(mesh);
        // Add distance from origin into the node "point" data
        std::vector<double> rank;
        // Real rank for the owned nodes
        for (unsigned i=0; i<mesh.GetNumLocalNodes(); i++)
        {
            rank.push_back(PetscTools::GetMyRank());
        }
        writer.AddPointData("Process rank", rank);

        writer.WriteFilesUsingMesh(mesh);

        PetscTools::Barrier("Wait for files to be written");

        std::stringstream filepath;
        filepath << OutputFileHandler::GetChasteTestOutputDirectory() << "TestVtkMeshWriter/1D_0_to_1_10_elements_parallel_data";
        // Add suffix to VTK vtu file.
        if (PetscTools::IsSequential())
        {
            filepath <<  ".vtu";
        }
        else
        {
            /*
             * Check that the pvtu file exists. Note that checking its content is hard
             * because the number of processes (.vtu file references) will vary.
             */
            FileFinder vtk_file(filepath.str() + ".pvtu", RelativeTo::Absolute);
            TS_ASSERT(vtk_file.Exists());
            // Add suffix to VTK vtu file
            filepath << "_" << PetscTools::GetMyRank() << ".vtu";
        }
        {
            // Check that the reader can see it
            VtkMeshReader<1,1> vtk_reader(filepath.str());
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumLocalNodes() + mesh.GetNumHaloNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumLocalElements());

            // Check that it has the correct data
            std::vector<double> rank_read;
            vtk_reader.GetPointData("Process rank", rank_read);
            TS_ASSERT_EQUALS(rank.size(), mesh.GetNumLocalNodes() );
            TS_ASSERT_EQUALS(rank_read.size(), mesh.GetNumLocalNodes()  + mesh.GetNumHaloNodes());
            for (unsigned i=0; i<rank.size(); i++)
            {
                TS_ASSERT_EQUALS(rank[i], rank_read[i]);
                TS_ASSERT_EQUALS(rank_read[i], PetscTools::GetMyRank());
            }
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestVtkMeshWriter1D()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_10_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<1,1> writer("TestVtkMeshWriter", "1D_0_to_1_10_elements", false);

        writer.WriteFilesUsingMesh(mesh);

        if (PetscTools::AmMaster())
        {
            // Check that the reader can see it
            VtkMeshReader<1,1> vtk_reader(OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/1D_0_to_1_10_elements.vtu");
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumElements());

            //Construct a mesh
            TetrahedralMesh<1,1> read_mesh;
            read_mesh.ConstructFromMeshReader(vtk_reader);
            TS_ASSERT_DELTA(read_mesh.GetNode(5)->rGetLocation()[0], 0.5, 1e-8);
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }


    void TestVtkMeshWriterWithData()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("heart/test/data/UCSD_heart_decimated_173nodes");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "heart_decimation", false);

        // Add element quality into the element "cell" data
        std::vector<double> quality;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            quality.push_back(mesh.GetElement(i)->CalculateQuality());
        }
        writer.AddCellData("Quality", quality);

        // Add fibre type to "cell" data
        std::vector< c_vector<double, 3> > centroid;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            centroid.push_back(mesh.GetElement(i)->CalculateCentroid());
        }
        writer.AddCellData("Centroid", centroid);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        // Add fibre type to "point" data
        std::vector< c_vector<double, 3> > location;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            location.push_back(mesh.GetNode(i)->rGetLocation());
        }
        writer.AddPointData("Location", location);

        // Add Jacobian tensor cell data, covering both the symmetric and non-symmetric tensor methods
        std::vector< c_matrix<double, 3, 3> > jacobian;
        std::vector< c_vector<double, 6> > squared_jacobian;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            c_matrix<double, 3, 3> element_jacobian;
            double element_jacobian_determinant;
            mesh.GetElement(i)->CalculateJacobian(element_jacobian, element_jacobian_determinant);
            jacobian.push_back(element_jacobian);

            c_matrix<double, 3, 3> squared_element_jacobian = prod(element_jacobian, element_jacobian);
            c_vector<double, 6> tri_squared_element_jacobian;
            //We store [T00 T01 T02 T11 T12 T22]
            tri_squared_element_jacobian(0) = squared_element_jacobian(0, 0);
            tri_squared_element_jacobian(1) = squared_element_jacobian(0, 1);
            tri_squared_element_jacobian(2) = squared_element_jacobian(0, 2);
            tri_squared_element_jacobian(3) = squared_element_jacobian(1, 1);
            tri_squared_element_jacobian(4) = squared_element_jacobian(1, 2);
            tri_squared_element_jacobian(5) = squared_element_jacobian(2, 2);

            squared_jacobian.push_back(tri_squared_element_jacobian);
        }
        writer.AddTensorCellData("Jacobian", jacobian);
        writer.AddTensorCellData("SquaredJacobian", squared_jacobian);

        //Add tensor point data, we use the outer product of the node's location
        std::vector< c_matrix<double, 3, 3> > location_outer_product;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_matrix<double, 3, 3> element_location_outer_product;
            element_location_outer_product = outer_prod(trans(mesh.GetNode(i)->rGetLocation()), mesh.GetNode(i)->rGetLocation());
            location_outer_product.push_back(element_location_outer_product);
        }
        writer.AddTensorPointData("LocationProduct", location_outer_product);

        writer.WriteFilesUsingMesh(mesh);
        //32K uncompressed, 19K compressed


        {
            // Check that the reader can see it
            VtkMeshReader<3,3> vtk_reader(OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/heart_decimation.vtu");
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumElements());

            // Check that it has the correct data
            std::vector<double> distance_read;
            vtk_reader.GetPointData("Distance from origin", distance_read);
            for (unsigned i=0; i<distance_read.size(); i++)
            {
                TS_ASSERT_EQUALS(distance[i], distance_read[i]);
            }
            std::vector<double> quality_read;
            vtk_reader.GetCellData("Quality", quality_read);
            for (unsigned i=0; i<quality_read.size(); i++)
            {
                TS_ASSERT_EQUALS(quality[i], quality_read[i]);
            }

            std::vector<c_vector<double,3> > centroid_read;
            vtk_reader.GetCellData("Centroid", centroid_read);
            TS_ASSERT_EQUALS(centroid_read.size(),centroid.size());
            ///\todo #1731 - need to read the tensors too.
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestVtkMeshWriterForCables()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/cylinder");
        TrianglesMeshReader<3,3> reader(mesh_base);
        MixedDimensionMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "mixed_mesh_3d", false);
        writer.SetParallelFiles(mesh);

        // Add element quality into the element "cell" data
         std::vector<double> quality;
         for (DistributedTetrahedralMesh<3,3>::ElementIterator ele_iter = mesh.GetElementIteratorBegin();
                ele_iter != mesh.GetElementIteratorEnd();
                ++ele_iter)
         {
             quality.push_back(ele_iter->CalculateQuality());
         }
         writer.AddCellData("Quality", quality);

        // Add fibre type to "cell" data
        std::vector< c_vector<double, 3> > centroid;
        for (DistributedTetrahedralMesh<3,3>::ElementIterator ele_iter = mesh.GetElementIteratorBegin();
               ele_iter != mesh.GetElementIteratorEnd();
               ++ele_iter)
        {
            centroid.push_back(ele_iter->CalculateCentroid());
        }
        writer.AddCellData("Centroid", centroid);

        writer.WriteFilesUsingMesh(mesh);

        ///\todo #2052 We can't yet test if the cables are written correctly, because we don't have the reader part.
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestVtkMeshWriterForQuadraticMesh2D()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_200_elements");
        QuadraticMesh<2> mesh;
        mesh.ConstructFromLinearMeshReader(reader);

        VtkMeshWriter<2,2> writer("TestVtkMeshWriter", "2D_0_to_1mm_200_elements_quadratic", false);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        // Add fibre type to "point" data
        std::vector< c_vector<double, 2> > location;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            location.push_back(mesh.GetNode(i)->rGetLocation());
        }
        writer.AddPointData("Location", location);

        // Add element quality into the element "cell" data
        std::vector<double> quality;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            quality.push_back(mesh.GetElement(i)->CalculateQuality());
        }
        writer.AddCellData("Quality", quality);

        // Add fibre type to "cell" data
        std::vector< c_vector<double, 2> > centroid;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            centroid.push_back(mesh.GetElement(i)->CalculateCentroid());
        }
        writer.AddCellData("Centroid", centroid);


        writer.WriteFilesUsingMesh(mesh);

        {
            // Check that the reader can see it
            VtkMeshReader<2,2> vtk_reader(OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/2D_0_to_1mm_200_elements_quadratic.vtu");
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumElements());

            // Check that it has the correct data
            std::vector<double> distance_read;
            vtk_reader.GetPointData("Distance from origin", distance_read);
            for (unsigned i=0; i<distance_read.size(); i++)
            {
                TS_ASSERT_EQUALS(distance[i], distance_read[i]);
            }
            std::vector<c_vector<double,2> > location_read;
            vtk_reader.GetPointData("Location", location_read);
            for (unsigned i=0; i<location_read.size(); i++)
            {
                for (unsigned j=0; j<2; j++)
                {
                    TS_ASSERT_EQUALS(location[i][j], location_read[i][j]);
                }
            }
            std::vector<double> quality_read;
            vtk_reader.GetCellData("Quality", quality_read);
            for (unsigned i=0; i<quality_read.size(); i++)
            {
                TS_ASSERT_EQUALS(quality[i], quality_read[i]);
            }
            std::vector<c_vector<double,2> > centroid_read;
            vtk_reader.GetCellData("Centroid", centroid_read);
            for (unsigned i=0; i<centroid_read.size(); i++)
            {
                for (unsigned j=0; j<2; j++)
                {
                    TS_ASSERT_EQUALS(centroid[i][j], centroid_read[i][j]);
                }
            }
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestBasicQuadraticVtkMeshWriter()
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        QuadraticMesh<3> mesh;
        mesh.ConstructFromLinearMeshReader(reader);

        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "cube_2mm_12_elements_quadratic", false);

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        //1.6K uncompressed, 1.3K compressed
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/";

        {
            //Check that the reader can see it
            VtkMeshReader<3,3> vtk_reader(results_dir+"cube_2mm_12_elements_quadratic.vtu");
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumElements());

            ///\todo: The reader can open a quadratic vtu file, but not construct a QuadraticMesh
            //further tests of the written output should be made once this is supported.
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    //Test that the vtk mesh writer can output a 1D mesh embedded in 3D space
    void TestVtkMeshWriter1Din3D()
    {
#ifdef CHASTE_VTK
        TrianglesMeshReader<1,3> reader("mesh/test/data/branched_1d_in_3d_mesh");
        TetrahedralMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<1,3> writer("TestVtkMeshWriter", "branched_1d_in_3d_mesh", false);

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        {
            // Check that the reader can see it
            VtkMeshReader<1,3> vtk_reader(OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/branched_1d_in_3d_mesh.vtu");
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumElements());

            //There should be 3 terminals (boundary nodes/elements)
            TS_ASSERT_EQUALS(vtk_reader.GetNumEdges(), reader.GetNumEdges());
            TS_ASSERT_EQUALS(vtk_reader.GetNumEdges(), mesh.GetNumBoundaryElements());
        }

#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    //Test that the vtk mesh writer can output a 2D mesh embedded in 3D space
    void TestVtkMeshWriterWithSurfaceMesh()
    {
#ifdef CHASTE_VTK
        VtkMeshReader<2,3> mesh_reader("mesh/test/data/cylinder.vtu");
        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        VtkMeshWriter<2,3> writer("TestVtkMeshWriter", "cylinder", false);

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        {
            // Check that the reader can see it
            VtkMeshReader<2,3> vtk_reader(OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/cylinder.vtu");
            TS_ASSERT_EQUALS(vtk_reader.GetNumNodes(), mesh.GetNumNodes());
            TS_ASSERT_EQUALS(vtk_reader.GetNumElements(), mesh.GetNumElements());

            TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 1632u);
            TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 32u);
        }

#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestXdmfWriter()
    {
        /*Read as ascii*/
        TrianglesMeshReader<3,3> reader("mesh/test/data/simple_cube");

        XdmfMeshWriter<3,3> writer_from_reader("TestXdmfMeshWriter", "simple_cube", false);
#ifdef _MSC_VER
        if (PetscTools::AmMaster()) // Only the master does anything when using a mesh reader
        {
            TS_ASSERT_THROWS_THIS(writer_from_reader.WriteFilesUsingMeshReader(reader), "XDMF is not supported under Windows at present.");
        }
#else
        writer_from_reader.WriteFilesUsingMeshReader(reader);

        //Check that the files are the same as previously
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube.xdmf",
                       "mesh/test/data/TestXdmfMeshWriter/simple_cube.xdmf").CompareFiles();
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube_geometry_0.xml",
                       "mesh/test/data/TestXdmfMeshWriter/simple_cube_geometry_0.xml").CompareFiles();
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube_topology_0.xml",
                       "mesh/test/data/TestXdmfMeshWriter/simple_cube_topology_0.xml").CompareFiles();
#endif // _MSC_VER

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);
        XdmfMeshWriter<3,3> writer_from_mesh("TestXdmfMeshWriter", "simple_cube_from_mesh", false);
#ifdef _MSC_VER
        TS_ASSERT_THROWS_THIS(writer_from_mesh.WriteFilesUsingMesh(mesh), "XDMF is not supported under Windows at present.");
#else
        writer_from_mesh.WriteFilesUsingMesh(mesh);

        //Just check that the files are there
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube_from_mesh.xdmf",
                       "mesh/test/data/TestXdmfMeshWriter/simple_cube_from_mesh.xdmf").CompareFiles();
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube_from_mesh_geometry_0.xml",
                       "mesh/test/data/TestXdmfMeshWriter/simple_cube_from_mesh_geometry_0.xml").CompareFiles();
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube_from_mesh_topology_0.xml",
                       "mesh/test/data/TestXdmfMeshWriter/simple_cube_from_mesh_topology_0.xml").CompareFiles();
#endif // _MSC_VER
    }

    void TestXdmfWriter2D()
    {
#ifndef _MSC_VER
        /*Read as ascii*/
        TrianglesMeshReader<2,2> reader("mesh/test/data/disk_522_elements");

        XdmfMeshWriter<2,2> writer_from_reader("TestXdmfMeshWriter", "disk", false);
        writer_from_reader.WriteFilesUsingMeshReader(reader);

        //Check that the files are the same as previously
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/disk.xdmf",
                       "mesh/test/data/TestXdmfMeshWriter/disk.xdmf").CompareFiles();
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/disk_geometry_0.xml",
                       "mesh/test/data/TestXdmfMeshWriter/disk_geometry_0.xml").CompareFiles();
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/disk_topology_0.xml",
                       "mesh/test/data/TestXdmfMeshWriter/disk_topology_0.xml").CompareFiles();

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        XdmfMeshWriter<2,2> writer_from_mesh("TestXdmfMeshWriter", "disk_from_mesh", false);
        writer_from_mesh.WriteFilesUsingMesh(mesh);

        //Just check that the files are there
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/disk_from_mesh.xdmf",
                       "mesh/test/data/TestXdmfMeshWriter/disk_from_mesh.xdmf").CompareFiles();
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/disk_from_mesh_geometry_0.xml",
                       "mesh/test/data/TestXdmfMeshWriter/disk_from_mesh_geometry_0.xml").CompareFiles();
        FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/disk_from_mesh_topology_0.xml",
                       "mesh/test/data/TestXdmfMeshWriter/disk_from_mesh_topology_0.xml").CompareFiles();
#endif // _MSC_VER
    }

    void TestXdmfWriterDistributed()
    {
#ifndef _MSC_VER
        TrianglesMeshReader<3,3> reader("mesh/test/data/simple_cube");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        XdmfMeshWriter<3,3> writer_from_mesh("TestXdmfMeshWriter", "simple_cube_dist", false);
        writer_from_mesh.WriteFilesUsingMesh(mesh);

        if (PetscTools::AmMaster())
        {
            for (unsigned i=0; i<PetscTools::GetNumProcs(); i++)
            {
                std::stringstream chunk_name;
                chunk_name << OutputFileHandler::GetChasteTestOutputDirectory();
                chunk_name << "TestXdmfMeshWriter/simple_cube_dist_topology_" << i <<".xml";
                TS_ASSERT( FileFinder(chunk_name.str(), RelativeTo::Absolute).Exists());
            }
            // If we are running with exactly 2 processes, then we can check for the exact output
            if (PetscTools::GetNumProcs() == 2u)
            {
                FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube_dist.xdmf",
                               "mesh/test/data/TestXdmfMeshWriter/simple_cube_dist.xdmf", false /*not collective*/).CompareFiles();
                FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube_dist_geometry_0.xml",
                               "mesh/test/data/TestXdmfMeshWriter/simple_cube_dist_geometry_0.xml", false /*not collective*/).CompareFiles();
                FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube_dist_topology_0.xml",
                               "mesh/test/data/TestXdmfMeshWriter/simple_cube_dist_topology_0.xml", false /*not collective*/).CompareFiles();
                FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube_dist_geometry_1.xml",
                                               "mesh/test/data/TestXdmfMeshWriter/simple_cube_dist_geometry_1.xml", false /*not collective*/).CompareFiles();
                FileComparison(OutputFileHandler::GetChasteTestOutputDirectory() + "TestXdmfMeshWriter/simple_cube_dist_topology_1.xml",
                                               "mesh/test/data/TestXdmfMeshWriter/simple_cube_dist_topology_1.xml", false /*not collective*/).CompareFiles();
            }
        }
#endif // _MSC_VER
     }
};

#endif //_TESTXMLMESHWRITERS_HPP_
