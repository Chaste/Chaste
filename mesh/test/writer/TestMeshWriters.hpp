/*

Copyright (c) 2005-2012, University of Oxford.
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


#ifndef _TESTMESHWRITERS_HPP_
#define _TESTMESHWRITERS_HPP_

#include <cxxtest/TestSuite.h>
#include "MemfemMeshReader.hpp"
#include "FemlabMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "MutableMesh.hpp"
#include "CmguiMeshWriter.hpp"
#include "CmguiDeformedSolutionsWriter.hpp"
#include "QuadraticMesh.hpp"
#include <iostream>

class TestMeshWriters : public CxxTest::TestSuite
{

public:

    void TestMemfemToTetgen()
    {
        TrianglesMeshWriter<3,3> mesh_writer("", "MeshFromMemfem");
        MemfemMeshReader<3,3> import_mesh_reader("mesh/test/data/Memfem_slab");

        mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader);
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<3,3>* p_new_mesh_reader;
        p_new_mesh_reader = new TrianglesMeshReader<3,3>(output_dir + "MeshFromMemfem");

        delete p_new_mesh_reader;
    }

    void TestFemlabToTriangles()
    {
        TrianglesMeshWriter<2,2> mesh_writer("", "MeshFromFemlab");

        FemlabMeshReader<2,2> import_mesh_reader("mesh/test/data/",
                                                 "femlab_lshape_nodes.dat",
                                                 "femlab_lshape_elements.dat",
                                                 "femlab_lshape_edges.dat");

        //TS_ASSERT_EQUALS(import_mesh_reader.GetNumFaces(), 54U); //Has internal faces
        TS_ASSERT_EQUALS(import_mesh_reader.GetNumFaces(), 40U); //Has no internal faces
        mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader);
        std::string output_dir = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<2,2>* p_new_mesh_reader;
        p_new_mesh_reader = new TrianglesMeshReader<2,2>(output_dir + "MeshFromFemlab");
        //TS_ASSERT_EQUALS(p_new_mesh_reader->GetNumFaces(), 54U); //No faces have been culled
        TS_ASSERT_EQUALS(p_new_mesh_reader->GetNumFaces(), 40U); //No faces have been culled

        delete p_new_mesh_reader;
    }

    void TestTrianglesToMeshalyzer1d()
    {
        TrianglesMeshReader<1,1> import_mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MeshalyzerMeshWriter<1,1> mesh_writer("", "MeshFromTetgen");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }

    void TestTrianglesToMeshalyzer2d()
    {
        TrianglesMeshReader<2,2> import_mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        MeshalyzerMeshWriter<2,2> mesh_writer("", "MeshFromTetgen");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }

    void TestTrianglesToMeshalyzer3d()
    {
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        MeshalyzerMeshWriter<3,3> mesh_writer("", "MeshFromTetgen");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }

    void TestTrianglesToMeshalyzer1dIn3d()
    {
        TrianglesMeshReader<1,3> import_mesh_reader("mesh/test/data/trivial_1d_in_3d_mesh");
        MeshalyzerMeshWriter<1,3> mesh_writer("", "Mesh");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }

    void TestTrianglesToCoolGraphics()
    {
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        bool set_CG_format = true;

        MeshalyzerMeshWriter<3,3> mesh_writer("CGFromTetgen", "CGFromTetgen", true, set_CG_format);

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }

    void TestFemlabtoTrianglesViaMesh()
    {
        TrianglesMeshWriter<2,2> mesh_writer("", "MeshFromFemlabViaMesh");

        FemlabMeshReader<2,2> import_mesh_reader("mesh/test/data/",
                                                 "femlab_lshape_nodes.dat",
                                                 "femlab_lshape_elements.dat",
                                                 "femlab_lshape_edges.dat");

        TS_ASSERT_EQUALS(import_mesh_reader.GetNumFaces(), 40U);//Has no internal faces
        TetrahedralMesh<2,2> mesh;

        mesh.ConstructFromMeshReader(import_mesh_reader);

        mesh_writer.WriteFilesUsingMesh(mesh);
        std::string output_dir = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<2,2>* p_new_mesh_reader;
        p_new_mesh_reader = new TrianglesMeshReader<2,2>(output_dir + "MeshFromFemlabViaMesh");
        TS_ASSERT_EQUALS(p_new_mesh_reader->GetNumFaces(), 40U); //Internal faces have been culled

        delete p_new_mesh_reader;
    }

    void TestTrianglesToMeshalyzerViaMesh1d()
    {
        TrianglesMeshReader<1,1> import_mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MeshalyzerMeshWriter<1,1> mesh_writer("", "MeshFromTetgenViaMesh");

        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }

    void TestTrianglesToMeshalyzerViaMesh2d()
    {
        TrianglesMeshReader<2,2> import_mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        MeshalyzerMeshWriter<2,2> mesh_writer("", "MeshFromTetgenViaMesh");

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }

    void TestTrianglesToMeshalyzerViaMesh3d()
    {
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        MeshalyzerMeshWriter<3,3> mesh_writer("", "MeshFromTetgenViaMesh");

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        unsigned num_nodes = mesh.GetNumNodes();
        unsigned num_nodes_from_file;

        std::string filename = OutputFileHandler::GetChasteTestOutputDirectory() + "MeshFromTetgenViaMesh.pts";
        std::ifstream meshalyzer_output_file;
        meshalyzer_output_file.open(filename.c_str());
        assert(meshalyzer_output_file.is_open());
        meshalyzer_output_file >> num_nodes_from_file;

        TS_ASSERT_EQUALS(num_nodes,num_nodes_from_file);

        // Nothing appears to be tested here, so commenting it out to remove a system call.
//        std::string command = "wc -l " + filename;
//        int result = system(command.c_str());
//        std::cout<<result<<std::endl;
    }

    void TestTrianglesToCoolGraphicsViaMesh()
    {
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        bool set_CG_format = true;
        MeshalyzerMeshWriter<3,3> mesh_writer("CGFromTetgenViaMesh", "CGFromTetgenViaMesh", true, set_CG_format);

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }

    void TestTriangles1DClosedMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/circle_outline");
        TetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<1,2> mesh_writer("", "1dClosedMeshIn2dSpace");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();

        std::cout << output_dir + "1dClosedMeshIn2dSpace" << std::endl;

        TrianglesMeshReader<1,2> mesh_reader2(output_dir + "1dClosedMeshIn2dSpace");

        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 100u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 100u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 0u);
    }

    void TestTriangles1DMeshIn2DSpaceBinary()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        TetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<1,2> mesh_writer("", "1dMeshIn2dSpace");
        mesh_writer.SetWriteFilesAsBinary();

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<1,2> mesh_reader2(output_dir + "1dMeshIn2dSpace");

        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 51u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 50u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 2u);

        // Test for connectivity
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + output_dir + "/1dMeshIn2dSpace.ncl mesh/test/data/1dMeshIn2dSpace.ncl").c_str()), 0);
    }

    void TestTriangles1DMeshIn2DSpaceWithDeletedNode() throw (Exception)
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        MutableMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.DeleteBoundaryNodeAt(0);

        TrianglesMeshWriter<1,2> mesh_writer("", "1dMeshIn2dSpaceWithDeletedNode");
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        TrianglesMeshWriter<1,2> mesh_writer2("", "1dMeshIn2dSpaceWithDeletedNodeConst");
        //mesh_writer2.WriteFilesUsingMesh(static_cast<const MutableMesh<1,2> &>(mesh));
        mesh_writer2.WriteFilesUsingMesh(mesh);

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<1,2> mesh_reader2(output_dir + "1dMeshIn2dSpaceWithDeletedNode");
        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 50u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 49u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 2u);

        TrianglesMeshReader<1,2> mesh_reader3(output_dir + "1dMeshIn2dSpaceWithDeletedNodeConst");
        TS_ASSERT_EQUALS(mesh_reader3.GetNumNodes(), 50u);
        TS_ASSERT_EQUALS(mesh_reader3.GetNumElements(), 49u);
        TS_ASSERT_EQUALS(mesh_reader3.GetNumFaces(), 2u);
    }

    void Test2DClosedMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/slab_395_elements");

        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<2,3> mesh_writer("", "2dClosedMeshIn3dSpace");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,3> mesh_reader2(output_dir+"2dClosedMeshIn3dSpace");

        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 132u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 224u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 0u);
    }

    void Test2DMeshIn3DSpaceBinary()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");

        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<2,3> mesh_writer("", "2dMeshIn3dSpace");
        mesh_writer.SetWriteFilesAsBinary();

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,3> mesh_reader2(output_dir + "2dMeshIn3dSpace");

        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 312u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 522u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), mesh_reader.GetNumElements());

        // Note that this not a straight conversion, since we have culled internal data
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 100u); // culling now occurs in the reader
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 100u);
    }

    void TestQuadratic1D() throw (Exception)
    {
        QuadraticMesh<1> mesh;
        TrianglesMeshReader<1,1> mesh_reader1("mesh/test/data/1D_0_to_1_10_elements_quadratic", 2, 1, false);
        mesh.ConstructFromMeshReader(mesh_reader1);
        TrianglesMeshWriter<1,1> mesh_writer("", "1d_quadratic");
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<1,1> mesh_reader(output_dir + "1d_quadratic", 2, 2);

        // Test that reader is reading correctly
        TS_ASSERT_EQUALS(mesh_reader.GetNextNode().size(), 1u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextElementData().NodeIndices.size(), 3u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().NodeIndices.size(), 1u);

        // Test that mesh can be reconstructed
        QuadraticMesh<1> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNumNodes(), 1U);
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNodeGlobalIndex(0), mesh.GetBoundaryElement(0)->GetNodeGlobalIndex(0));
    }

    void TestQuadratic2D() throw (Exception)
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader1("mesh/test/data/square_128_elements_fully_quadratic", 2, 2, false);
        mesh.ConstructFromMeshReader(mesh_reader1);
        TrianglesMeshWriter<2,2> mesh_writer("", "2d_quadratic");
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,2> mesh_reader(output_dir + "2d_quadratic", 2, 2);

        // Test that reader is reading correctly
        TS_ASSERT_EQUALS(mesh_reader.GetNextNode().size(), 2u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextElementData().NodeIndices.size(), 6u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().NodeIndices.size(), 3u);

        // Test that mesh can be reconstructed
        QuadraticMesh<2> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNumNodes(), 3U);
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNodeGlobalIndex(2), mesh.GetBoundaryElement(0)->GetNodeGlobalIndex(2));
    }

    void TestQuadratic3D() throw (Exception)
    {
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader1("mesh/test/data/3D_Single_tetrahedron_element_quadratic", 2, 1, false);
        mesh.ConstructFromMeshReader(mesh_reader1);
        TrianglesMeshWriter<3,3> mesh_writer("", "3d_quadratic");
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<3,3> mesh_reader(output_dir + "3d_quadratic", 2, 2);

        // Test that reader is reading correctly
        TS_ASSERT_EQUALS(mesh_reader.GetNextNode().size(), 3u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextElementData().NodeIndices.size(), 10u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().NodeIndices.size(), 6u);

        // Test that mesh can be reconstructed
        QuadraticMesh<3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNumNodes(), 6U);
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNodeGlobalIndex(5), mesh.GetBoundaryElement(0)->GetNodeGlobalIndex(5));
    }

    void TestCmguiMeshWriter3D() throw(Exception)
    {
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        CmguiMeshWriter<3,3> writer("TestCmguiMeshWriter3D", "cube_2mm_12_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriter3D/";

        std::string node_file1 = results_dir + "/cube_2mm_12_elements.exnode";
        std::string node_file2 = "mesh/test/data/TestCmguiMeshWriter/cube_2mm_12_elements.exnode";
        std::string elem_file1 = results_dir + "/cube_2mm_12_elements.exelem";
        std::string elem_file2 = "mesh/test/data/TestCmguiMeshWriter/cube_2mm_12_elements.exelem";

        bool comparison_result = CmguiMeshWriter<3,3>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);
        comparison_result = CmguiMeshWriter<3,3>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);

        //now test the set method for additional fields. We set two fields.
        CmguiMeshWriter<3,3> writer2("TestCmguiMeshWriterAdditionalHeaders3D", "cube_2mm_12_elements");

        std::vector<std::string> field_names;
        field_names.push_back("V");
        field_names.push_back("Phi_e");
        std::string results_dir2 = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriterAdditionalHeaders3D";
        writer2.SetAdditionalFieldNames(field_names);
        TS_ASSERT_THROWS_NOTHING(writer2.WriteFilesUsingMesh(mesh));

        elem_file1 = results_dir2 + "/cube_2mm_12_elements.exelem";
        elem_file2 = "mesh/test/data/TestCmguiMeshWriter/cube_2mm_12_elements_additional_fields.exelem";

        comparison_result = CmguiMeshWriter<3,3>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);
    }

    void TestCmguiMeshWriter2D() throw(Exception)
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        CmguiMeshWriter<2,2> writer("TestCmguiMeshWriter2D", "square_128_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriter2D/";

        std::string node_file1 = results_dir + "/square_128_elements.exnode";
        std::string node_file2 = "mesh/test/data/TestCmguiMeshWriter/square_128_elements.exnode";
        std::string elem_file1 = results_dir + "/square_128_elements.exelem";
        std::string elem_file2 = "mesh/test/data/TestCmguiMeshWriter/square_128_elements.exelem";

        bool comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);

        comparison_result =  CmguiMeshWriter<2,2>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);

        //now test the set method for additional fields. We set two fields.
        CmguiMeshWriter<2,2> writer2("TestCmguiMeshWriterAdditionalHeaders2D", "square_128_elements");

        std::vector<std::string> field_names;
        field_names.push_back("V");
        field_names.push_back("Phi_e");
        std::string results_dir2 = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriterAdditionalHeaders2D";
        writer2.SetAdditionalFieldNames(field_names);
        TS_ASSERT_THROWS_NOTHING(writer2.WriteFilesUsingMesh(mesh));

        elem_file1 = results_dir2 + "/square_128_elements.exelem";
        elem_file2 = "mesh/test/data/TestCmguiMeshWriter/square_128_elements_additional_fields.exelem";

        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(elem_file1,elem_file2);

        TS_ASSERT(comparison_result);
    }

    void TestCmguiMeshWriter1D() throw(Exception)
    {
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        CmguiMeshWriter<1,1> writer("TestCmguiMeshWriter1D", "1D_0_to_1_100_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriter1D/";

        std::string node_file1 = results_dir + "/1D_0_to_1_100_elements.exnode";
        std::string node_file2 = "mesh/test/data/TestCmguiMeshWriter/1D_0_to_1_100_elements.exnode";
        std::string elem_file1 = results_dir + "/1D_0_to_1_100_elements.exelem";
        std::string elem_file2 = "mesh/test/data/TestCmguiMeshWriter/1D_0_to_1_100_elements.exelem";

        bool comparison_result =  CmguiMeshWriter<1,1>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);
        comparison_result =  CmguiMeshWriter<1,1>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);

        // Now test the set method for additional fields. We set two fields.
        CmguiMeshWriter<1,1> writer2("TestCmguiMeshWriterAdditionalHeaders1D", "1D_0_to_1_100_elements");

        std::vector<std::string> field_names;
        field_names.push_back("V");
        field_names.push_back("Phi_e");
        std::string results_dir2 = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriterAdditionalHeaders1D";
        writer2.SetAdditionalFieldNames(field_names);
        TS_ASSERT_THROWS_NOTHING(writer2.WriteFilesUsingMesh(mesh));

        elem_file1 = results_dir2 + "/1D_0_to_1_100_elements.exelem";
        elem_file2 = "mesh/test/data/TestCmguiMeshWriter/1D_0_to_1_100_elements_additional_fields.exelem";

        comparison_result =  CmguiMeshWriter<1,1>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);
    }

    void TestCmguiDeformedSolutionsWriter2dLinearViz() throw(Exception)
    {
        QuadraticMesh<2> mesh(0.5, 1.0, 1.0);
        mesh.Scale(1.0, 2.0); // historical reasons

        CmguiDeformedSolutionsWriter<2> writer("TestCmguiDeformedSolutionsWriter",
                                                "solution",
                                                mesh,
                                                WRITE_LINEAR_MESH);

        // writes solution_0.exnode and solution_0.exelem using the mesh
        writer.WriteInitialMesh();
        writer.WriteInitialMesh("some_other_name");

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiDeformedSolutionsWriter";

        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/solution_0.exnode " + results_dir + "/some_other_name.exnode").c_str()), 0);
        // we shouldn't really have called WriteInitialMesh() twice on the same writer. If we had used different writers
        // we could also check the exelem files are the same here (the second one isn't created due to WriteInitialMesh
        // being called twice

        // set up a deformed positions vector
        std::vector<c_vector<double,2> > deformed_positions(mesh.GetNumNodes(), zero_vector<double>(2));
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            deformed_positions[i](0) = 2*mesh.GetNode(i)->rGetLocation()[0];
            deformed_positions[i](1) = 3*mesh.GetNode(i)->rGetLocation()[1];
        }
        // write solution_1.exnode
        writer.WriteDeformationPositions(deformed_positions, 1);

        // ..and again
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            deformed_positions[i](0) = 4*mesh.GetNode(i)->rGetLocation()[0];
            deformed_positions[i](1) = 6*mesh.GetNode(i)->rGetLocation()[1];
        }
        // write solution_1.exnode
        writer.WriteDeformationPositions(deformed_positions, 2);
        // write LoadSolutions.com
        writer.WriteCmguiScript();

        deformed_positions.push_back(zero_vector<double>(2));

        TS_ASSERT_THROWS_CONTAINS(writer.WriteDeformationPositions(deformed_positions, 3), "The size of rDeformedPositions does not match the number of nodes in the mesh");

        std::string node_file1 = results_dir + "/solution_0.exnode";
        std::string node_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_0.exnode";
        std::string elem_file1 = results_dir + "/solution_0.exelem";
        std::string elem_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_0.exelem";

        bool comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);
        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);

        node_file1 = results_dir + "/solution_1.exnode";
        node_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_1.exnode";
        elem_file1 = results_dir + "/solution_2.exnode"; // Not actually an element file - a second node file to compare!
        elem_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_2.exnode";

        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);

        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);

        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/LoadSolutions.com mesh/test/data/TestCmguiDeformedSolutionsWriter/LoadSolutions.com").c_str()), 0);

        // Coverage
        writer.WriteCmguiScript("myfield");
    }

    void TestCmguiDeformedSolutionsWriter2dQuadraticViz() throw(Exception)
    {
        QuadraticMesh<2> mesh(0.5, 1.0, 2.0);

        // Displace the internal nodes so you can see if quadratically visualised or not
        for (unsigned i=mesh.GetNumVertices(); i<mesh.GetNumNodes(); i++)
        {
            mesh.GetNode(i)->rGetModifiableLocation()[0] += 0.02;
            mesh.GetNode(i)->rGetModifiableLocation()[1] += 0.02;
        }

        CmguiDeformedSolutionsWriter<2> writer("TestCmguiDeformedSolutionsWriterQuadViz",
                                                "solution_2d_quadviz",
                                                mesh,
                                                WRITE_QUADRATIC_MESH);

        // Writes solution_0.exnode and solution_0.exelem using the mesh
        writer.WriteInitialMesh();

        // Set up a deformed positions vector
        std::vector<c_vector<double,2> > deformed_positions(mesh.GetNumNodes(), zero_vector<double>(2));
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            deformed_positions[i](0) = 2*mesh.GetNode(i)->rGetLocation()[0];
            deformed_positions[i](1) = 3*mesh.GetNode(i)->rGetLocation()[1];
        }
        // Write solution_1.exnode
        writer.WriteDeformationPositions(deformed_positions, 1);

        // ..and again
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            deformed_positions[i](0) = 4*mesh.GetNode(i)->rGetLocation()[0];
            deformed_positions[i](1) = 6*mesh.GetNode(i)->rGetLocation()[1];
        }
        // Write solution_1.exnode
        writer.WriteDeformationPositions(deformed_positions, 2);
        // Write LoadSolutions.com
        writer.WriteCmguiScript();

        deformed_positions.push_back(zero_vector<double>(2));

        TS_ASSERT_THROWS_CONTAINS(writer.WriteDeformationPositions(deformed_positions, 3), "The size of rDeformedPositions does not match the number of nodes in the mesh");

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiDeformedSolutionsWriterQuadViz";

        std::string node_file1 = results_dir + "/solution_2d_quadviz_0.exnode";
        std::string node_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_2d_quadviz_0.exnode";
        std::string elem_file1 = results_dir + "/solution_2d_quadviz_0.exelem";
        std::string elem_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_2d_quadviz_0.exelem";

        bool comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);
        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);

        node_file1 = results_dir + "/solution_2d_quadviz_1.exnode";
        node_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_2d_quadviz_1.exnode";
        elem_file1 = results_dir + "/solution_2d_quadviz_2.exnode"; // Not actually an element file - a second node file to compare!
        elem_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_2d_quadviz_2.exnode";

        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);

        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);
    }

    void TestCmguiDeformedSolutionsWriter3dQuadraticViz() throw(Exception)
    {
        QuadraticMesh<3> mesh(0.5, 1.0, 2.0, 3.0);

        // Displace the internal nodes so you can see if quadratically visualised or not
        for (unsigned i=mesh.GetNumVertices(); i<mesh.GetNumNodes(); i++)
        {
            mesh.GetNode(i)->rGetModifiableLocation()[0] += 0.02;
            mesh.GetNode(i)->rGetModifiableLocation()[1] += 0.02;
            mesh.GetNode(i)->rGetModifiableLocation()[1] += 0.03;
        }

        CmguiDeformedSolutionsWriter<3> writer("TestCmguiDeformedSolutionsWriterQuadViz3d",
                                                "solution_3d_quadviz",
                                                mesh,
                                                WRITE_QUADRATIC_MESH);

        // Writes solution_0.exnode and solution_0.exelem using the mesh
        writer.WriteInitialMesh();

        // Set up a deformed positions vector
        std::vector<c_vector<double,3> > deformed_positions(mesh.GetNumNodes(), zero_vector<double>(3));
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            deformed_positions[i](0) = 2*mesh.GetNode(i)->rGetLocation()[0];
            deformed_positions[i](1) = 3*mesh.GetNode(i)->rGetLocation()[1];
            deformed_positions[i](2) = 0.5*mesh.GetNode(i)->rGetLocation()[2];
        }
        // Write solution_1.exnode
        writer.WriteDeformationPositions(deformed_positions, 1);

        // Write LoadSolutions.com
        writer.WriteCmguiScript();

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiDeformedSolutionsWriterQuadViz3d";

        std::string node_file1 = results_dir + "/solution_3d_quadviz_0.exnode";
        std::string node_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_3d_quadviz_0.exnode";
        std::string elem_file1 = results_dir + "/solution_3d_quadviz_0.exelem";
        std::string elem_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_3d_quadviz_0.exelem";

        bool comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);
        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);

        node_file1 = results_dir + "/solution_3d_quadviz_1.exnode";
        node_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_3d_quadviz_1.exnode";

        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);
    }

    void TestCmguiDeformedSolutionsWriterConvertOutput() throw(Exception)
    {
        QuadraticMesh<2> mesh(0.5, 1.0, 1.0);
        mesh.Scale(1.0, 2.0); // historical reasons

        CmguiDeformedSolutionsWriter<2> writer("TestCmguiDeformedSolutionsWriter_ConvertOutput",
                                               "solution",
                                               mesh,
                                               WRITE_LINEAR_MESH);

        // throws as file doesn't exist
        TS_ASSERT_THROWS_CONTAINS(writer.ConvertOutput("blahblahblah", "blah", 1), "Could not open file:");

        // convert some files
        writer.ConvertOutput("mesh/test/data/TestCmguiDeformedSolutionsWriter", "myoldsolution", 2);

        // myoldsolution_1.nodes and myoldsolution_2.nodes are in chaste format but match the two deformed_positions written in the
        // test above, so the new output should match the good output from the test above.

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiDeformedSolutionsWriter_ConvertOutput";
        std::string node_file1 = results_dir + "/solution_0.exnode";
        std::string node_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_0.exnode";
        std::string elem_file1 = results_dir + "/solution_0.exelem";
        std::string elem_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_0.exelem";

        bool comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);
        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);

        node_file1 = results_dir + "/solution_1.exnode";
        node_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_1.exnode";
        elem_file1 = results_dir + "/solution_2.exnode"; // Not actually an element file - a second node file to compare!
        elem_file2 = "mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_2.exnode";

        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);
        comparison_result = CmguiMeshWriter<2,2>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);

        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/LoadSolutions.com mesh/test/data/TestCmguiDeformedSolutionsWriter/LoadSolutions.com").c_str()), 0);
        // throws as is incomplete
        TS_ASSERT_THROWS_CONTAINS(writer.ConvertOutput("mesh/test/data/TestCmguiDeformedSolutionsWriter", "bad_myoldsolution", 1), "Error occurred when reading file");
    }

    /**
     * This test is based on TestTrianglesMeshReader.hpp TestReadingElementAttributes.
     */
    void TestWritingElementAttributesInTrianglesFormat() throw (Exception)
    {
        std::string source_mesh = "mesh/test/data/1D_0_to_1_10_elements_with_attributes";
        std::string output_dir = "element_attrs";
        std::string file_from_reader = "from_reader";
        std::string file_from_mesh = "from_mesh";
        std::string file_from_mesh_bin = "from_mesh_binary";

        // Firstly, write directly from a mesh reader
        {
            TrianglesMeshReader<1,1> mesh_reader(source_mesh);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 10u);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

            TrianglesMeshWriter<1,1> mesh_writer(output_dir, file_from_reader, true);
            mesh_writer.WriteFilesUsingMeshReader(mesh_reader);
        }


        // Next, write using a mesh object
        {
            TrianglesMeshReader<1,1> mesh_reader(source_mesh);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 10u);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

            TetrahedralMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            TrianglesMeshWriter<1,1> mesh_writer(output_dir, file_from_mesh, false);
            mesh_writer.WriteFilesUsingMesh(mesh);
        }

        // Next, write using a mesh object in binary
        {
            TrianglesMeshReader<1,1> mesh_reader(source_mesh);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 10u);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

            TetrahedralMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            TrianglesMeshWriter<1,1> mesh_writer(output_dir, file_from_mesh_bin, false);
            mesh_writer.SetWriteFilesAsBinary();
            mesh_writer.WriteFilesUsingMesh(mesh);
        }

        // Now check the written meshes
        OutputFileHandler handler(output_dir, false);
        TrianglesMeshReader<1,1> reader1(handler.GetOutputDirectoryFullPath() + file_from_reader);
        TS_ASSERT_EQUALS(reader1.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(reader1.GetNumElementAttributes(), 1u);
        TrianglesMeshReader<1,1> reader2(handler.GetOutputDirectoryFullPath() + file_from_mesh);
        TS_ASSERT_EQUALS(reader2.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(reader2.GetNumElementAttributes(), 1u);
        TrianglesMeshReader<1,1> reader3(handler.GetOutputDirectoryFullPath() + file_from_mesh_bin);
        TS_ASSERT_EQUALS(reader3.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(reader3.GetNumElementAttributes(), 1u);
        for (unsigned i=0; i<10; i++)
        {
            ElementData next_element_info = reader1.GetNextElementData();
            std::vector<unsigned> nodes = next_element_info.NodeIndices;
            TS_ASSERT_EQUALS(nodes.size(), 2u);
            TS_ASSERT_EQUALS(next_element_info.AttributeValue, i%5 + 1);

            next_element_info = reader2.GetNextElementData();
            nodes = next_element_info.NodeIndices;
            TS_ASSERT_EQUALS(nodes.size(), 2u);
            TS_ASSERT_EQUALS(next_element_info.AttributeValue, i%5 + 1);

            next_element_info = reader3.GetNextElementData();
            nodes = next_element_info.NodeIndices;
            TS_ASSERT_EQUALS(nodes.size(), 2u);
            TS_ASSERT_EQUALS(next_element_info.AttributeValue, i%5 + 1);
        }
    }
    void TestWritingBinaryFormat()
    {
        /*Read as ascii*/
        TrianglesMeshReader<3,3> reader("mesh/test/data/simple_cube");

        TrianglesMeshWriter<3,3> writer_from_reader("TestMeshWriter", "simple_cube_binary_from_reader", false);
        writer_from_reader.SetWriteFilesAsBinary();
        writer_from_reader.WriteFilesUsingMeshReader(reader);

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);
        TrianglesMeshWriter<3,3> writer_from_mesh("TestMeshWriter", "simple_cube_binary_from_mesh", false);
        writer_from_mesh.SetWriteFilesAsBinary();
        writer_from_mesh.WriteFilesUsingMesh(mesh);
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestMeshWriter/";

        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/simple_cube_binary_from_reader.node mesh/test/data/simple_cube_binary.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/simple_cube_binary_from_reader.ele mesh/test/data/simple_cube_binary.ele").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/simple_cube_binary_from_reader.face mesh/test/data/simple_cube_binary.face").c_str()), 0);

        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/simple_cube_binary_from_mesh.node mesh/test/data/simple_cube_binary.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/simple_cube_binary_from_mesh.ele mesh/test/data/simple_cube_binary.ele").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/simple_cube_binary_from_mesh.face mesh/test/data/simple_cube_binary.face").c_str()), 0);

        // Test for connectivity
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/simple_cube_binary_from_mesh.ncl mesh/test/data/simple_cube_binary.ncl").c_str()), 0);

        /* Looking for beginning of provenance line: "#Created by Chaste"
         *          Ascii   Binary
         *  .node   182     267
         *  .ele    196     290
         *  .face   170     192
         *
         * Files are bigger!  That's because "\t1" is 2 bytes, but an unsigned "1" takes up
         * 4 bytes in an element file and a double "1" takes 8 bytes in a node file.
         * This won't be an issue on larger meshes with lots of floats.
            TS_ASSERT_EQUALS(sizeof(char), 1u);
            TS_ASSERT_EQUALS(sizeof(unsigned), 4u);
            TS_ASSERT_EQUALS(sizeof(double), 8u);
         */
    }
};

#endif //_TESTMESHWRITERS_HPP_
