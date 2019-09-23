/*

Copyright (C) Fujitsu Laboratories of Europe, 2009

*/

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

#ifndef TESTVTKMESHREADER_
#define TESTVTKMESHREADER_

#include <cxxtest/TestSuite.h>
#include <fstream>

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#endif //CHASTE_VTK

#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "MixedDimensionMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "VtkMeshReader.hpp"
#include "GenericMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UblasVectorInclude.hpp"

#ifdef CHASTE_VTK
typedef VtkMeshReader<3,3> MESH_READER3;
#endif //CHASTE_VTK

class TestVtkMeshReader : public CxxTest::TestSuite
{
//Requires  "sudo aptitude install libvtk5-dev" or similar
public:

    /**
     * Check that input files are opened correctly and non-existent input files throw an Exception.
     */
    void TestFilesOpen(void)
    {
#ifdef CHASTE_VTK
        TS_ASSERT_THROWS_NOTHING(MESH_READER3 mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu"));
        TS_ASSERT_THROWS_ANYTHING(MESH_READER3 mesh_reader("mesh/test/data/nofile.vtu"));
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    /**
     * Check outputting as a in VTKUnstructuredGrid format.
     */
    void TestOutputVtkUnstructuredGrid(void)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu");

        vtkUnstructuredGrid* vtk_unstructed_grid = mesh_reader.OutputMeshAsVtkUnstructuredGrid();

        TS_ASSERT_EQUALS(vtk_unstructed_grid->GetNumberOfPoints(), 12);
        TS_ASSERT_EQUALS(vtk_unstructed_grid->GetNumberOfCells(), 12);
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    /**
     * Check that the nodes are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing nodes) then an exception is thrown.
     */
    void TestGetNextNode(void)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu");

        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 12u);

        std::vector<double> first_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA(first_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[2], 0.0, 1e-6);

        std::vector<double> next_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA(next_node[0], 0.2, 1e-6);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(next_node[2], 0.0, 1e-6);

        for (unsigned node=2; node < mesh_reader.GetNumNodes(); node++)
        {
            next_node = mesh_reader.GetNextNode();
        }

        TS_ASSERT_THROWS_THIS( next_node = mesh_reader.GetNextNode(),
                               "Trying to read data for a node that doesn't exist" );
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestGetNextElementData(void)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu");

        // Coverage of GetOrderOfBoundaryElements()
        TS_ASSERT_EQUALS(mesh_reader.GetOrderOfBoundaryElements(), 1u);

        // Coverage of GetReadContainingElementOfBoundaryElement()
        TS_ASSERT_EQUALS(mesh_reader.GetReadContainingElementOfBoundaryElement(), false);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 12u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 0u);

        ElementData first_element_data = mesh_reader.GetNextElementData();
        TS_ASSERT_EQUALS(first_element_data.NodeIndices[0], 11u);
        TS_ASSERT_EQUALS(first_element_data.NodeIndices[1], 3u);
        TS_ASSERT_EQUALS(first_element_data.NodeIndices[2], 8u);
        TS_ASSERT_EQUALS(first_element_data.NodeIndices[3], 0u);
        TS_ASSERT_EQUALS(first_element_data.AttributeValue, 0u);

        ElementData next_element_data = mesh_reader.GetNextElementData();
        TS_ASSERT_EQUALS(next_element_data.NodeIndices[0], 10u);
        TS_ASSERT_EQUALS(next_element_data.NodeIndices[1], 8u);
        TS_ASSERT_EQUALS(next_element_data.NodeIndices[2], 5u);
        TS_ASSERT_EQUALS(next_element_data.NodeIndices[3], 11u);
        TS_ASSERT_EQUALS(next_element_data.AttributeValue, 0u);

        for (unsigned i=2; i<mesh_reader.GetNumElements(); i++)
        {
            next_element_data = mesh_reader.GetNextElementData();
            TS_ASSERT_EQUALS(next_element_data.AttributeValue, 0u);
        }

        TS_ASSERT_THROWS_THIS( mesh_reader.GetNextElementData(),
                               "Trying to read data for an element that doesn't exist" );

        // Test on a .vtu file where the elements are triangles, rather than tetrahedra
        VtkMeshReader<3,3> invalid_mesh_reader("mesh/test/data/sids.vtu");
        TS_ASSERT_EQUALS(invalid_mesh_reader.GetNumElements(), 736U);
        TS_ASSERT_EQUALS(invalid_mesh_reader.GetNumElementAttributes(), 0u);

        TS_ASSERT_THROWS_THIS( first_element_data = invalid_mesh_reader.GetNextElementData(),
                               "Element is not of expected type (vtkTetra/vtkTriangle)" );

#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestGetNextFaceData(void)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu");

        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 20u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumEdges(), 20u);

        ElementData first_face_data = mesh_reader.GetNextFaceData();
        TS_ASSERT_EQUALS(first_face_data.NodeIndices[0], 11u);
        TS_ASSERT_EQUALS(first_face_data.NodeIndices[1], 3u);
        TS_ASSERT_EQUALS(first_face_data.NodeIndices[2], 0u);

        ElementData next_face_data = mesh_reader.GetNextFaceData();
        TS_ASSERT_EQUALS(next_face_data.NodeIndices[0], 3u);
        TS_ASSERT_EQUALS(next_face_data.NodeIndices[1], 8u);
        TS_ASSERT_EQUALS(next_face_data.NodeIndices[2], 0u);

        mesh_reader.Reset();

        first_face_data = mesh_reader.GetNextEdgeData();
        TS_ASSERT_EQUALS(first_face_data.NodeIndices[0], 11u);
        TS_ASSERT_EQUALS(first_face_data.NodeIndices[1], 3u);
        TS_ASSERT_EQUALS(first_face_data.NodeIndices[2], 0u);

        next_face_data = mesh_reader.GetNextEdgeData();
        TS_ASSERT_EQUALS(next_face_data.NodeIndices[0], 3u);
        TS_ASSERT_EQUALS(next_face_data.NodeIndices[1], 8u);
        TS_ASSERT_EQUALS(next_face_data.NodeIndices[2], 0u);

        for (unsigned face=2; face<mesh_reader.GetNumFaces(); face++)
        {
            next_face_data = mesh_reader.GetNextEdgeData();
        }

        TS_ASSERT_THROWS_THIS( next_face_data = mesh_reader.GetNextEdgeData(),
                               "Trying to read data for a boundary element that doesn't exist" );
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestConstructFromVtkUnstructuredGridObject()
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader_1("mesh/test/data/cube_2mm_12_elements.vtu");
        vtkUnstructuredGrid* vtk_unstructed_grid = mesh_reader_1.OutputMeshAsVtkUnstructuredGrid();

        VtkMeshReader<3,3> mesh_reader(vtk_unstructed_grid);

        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 12u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 20u);

        std::vector<double> first_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA(first_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[2], 0.0, 1e-6);

        std::vector<double> next_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA(next_node[0], 0.2, 1e-6);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(next_node[2], 0.0, 1e-6);
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestGenericReader()
    {
#ifdef CHASTE_VTK
        std::shared_ptr<AbstractMeshReader<3,3> > p_mesh_reader = GenericMeshReader<3,3>("mesh/test/data/cube_2mm_12_elements.vtu");

        TS_ASSERT_EQUALS(p_mesh_reader->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_mesh_reader->GetNumElements(), 12u);
        TS_ASSERT_EQUALS(p_mesh_reader->GetNumFaces(), 20u);

        std::vector<double> first_node = p_mesh_reader->GetNextNode();
        TS_ASSERT_DELTA(first_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[2], 0.0, 1e-6);

        std::vector<double> next_node = p_mesh_reader->GetNextNode();
        TS_ASSERT_DELTA(next_node[0], 0.2, 1e-6);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(next_node[2], 0.0, 1e-6);

        // Exception coverage
        TS_ASSERT_THROWS_THIS((GenericMeshReader<3,3>("mesh/test/data/cube_2mm_12_elements.vtu", 2, 2)),
                              "Quadratic meshes are only supported in Triangles format.");
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    /**
     * Check that we can build a TetrahedralMesh using the mesh reader.
     */
    void TestBuildTetrahedralMeshFromMeshReader(void)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/heart_decimation.vtu");

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 173u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 610u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 312u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0963, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.3593, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[2], 0.9925, 1e-4);

        TS_ASSERT_DELTA(mesh.GetNode(8)->GetPoint()[0], 1.0969, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(8)->GetPoint()[1], 0.6678, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(8)->GetPoint()[2], 0.7250, 1e-4);

        // Check first element has the right nodes
        TetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(0), 47u);
        TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(1), 31u);
        TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(2), 131u);
        TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(3), 51u);
        TS_ASSERT_EQUALS((it)->GetNode(1), mesh.GetNode(31));

        /**
         * Check that point and cell data attributes work properly.
         */
        // Check element quality cell attribute is read properly
        std::vector<double> quality;
        TS_ASSERT_THROWS_THIS(mesh_reader.GetCellData( "Centroid", quality), "The cell data \'Centroid\' is not scalar data.");
        mesh_reader.GetCellData( "Quality", quality);
        for (unsigned i = 0; i < 610; i+=60)
        {

            TS_ASSERT_DELTA(quality[i], mesh.GetElement(i)->CalculateQuality(), 1e-4 );

        }

        // Check centroid attribute is read properly
        std::vector<c_vector<double, 3> > centroid;
        TS_ASSERT_THROWS_THIS(mesh_reader.GetCellData( "Quality", centroid), "The cell data \'Quality\' is not 3-vector data.");
        mesh_reader.GetCellData("Centroid", centroid);
        for (unsigned i = 0; i < 610; i+=60)
        {

            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(centroid[i][j], mesh.GetElement(i)->CalculateCentroid()[j], 1e-4);
            }
        }


        // Check distance from origin point attribute is read properly
        std::vector<double> distance;
        TS_ASSERT_THROWS_THIS(mesh_reader.GetPointData( "Location", distance), "The point data \'Location\' is not scalar data.");
        mesh_reader.GetPointData( "Distance from origin", distance);
        for (unsigned i = 0; i < 173; i+=17)
        {
            TS_ASSERT_DELTA(distance[i], norm_2(mesh.GetNode(i)->rGetLocation()), 1e-4);
        }

        // Check location attribute is read properly
        std::vector<c_vector<double, 3> > location;
        TS_ASSERT_THROWS_THIS(mesh_reader.GetPointData( "Distance from origin", location), "The point data \'Distance from origin\' is not 3-vector data.");
        mesh_reader.GetPointData("Location", location);
        for (unsigned i = 0; i < 173; i+=17)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(location[i][j], mesh.GetNode(i)->rGetLocation()[j], 1e-4);
            }
        }

        // Check that we can't ask for cell or point data that doesn't exist
        std::vector<double> not_there;
        std::vector<c_vector<double, 3> > vectors_not_there;
        TS_ASSERT_THROWS_ANYTHING( mesh_reader.GetCellData( "Non-existent data", not_there) );
        TS_ASSERT_THROWS_ANYTHING( mesh_reader.GetPointData( "Non-existent data", not_there) );
        TS_ASSERT_THROWS_ANYTHING( mesh_reader.GetCellData( "Non-existent data", vectors_not_there) );
        TS_ASSERT_THROWS_ANYTHING( mesh_reader.GetPointData( "Non-existent data", vectors_not_there) );
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    /**
     * Check that we can build a DistributedTetrahedralMesh using the VTK mesh reader.
     */
    void TestBuildDistributedTetrahedralMeshFromVtkMeshReader(void)
    {
 #ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/heart_decimation.vtu");

        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 173u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 610u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 312u);

        // Check some node co-ordinates
        try
        {
            Node<3> *node = mesh.GetNode(0);
            TS_ASSERT_DELTA(node->GetPoint()[0], 0.0963, 1e-4);
            TS_ASSERT_DELTA(node->GetPoint()[1], 0.3593, 1e-4);
            TS_ASSERT_DELTA(node->GetPoint()[2], 0.9925, 1e-4);
        }
        catch (Exception&)
        {
            // Don't own this node
        }

        try
        {
            Node<3> *node = mesh.GetNode(8);
            TS_ASSERT_DELTA(node->GetPoint()[0], 1.0969, 1e-4);
            TS_ASSERT_DELTA(node->GetPoint()[1], 0.6678, 1e-4);
            TS_ASSERT_DELTA(node->GetPoint()[2], 0.7250, 1e-4);
        }
        catch (Exception&)
        {
            // Don't own this node
        }

        // Check first element has the right nodes
        DistributedTetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(47))
        {
            //Owner of node 47 has to own (or part own) element 0
            TS_ASSERT_EQUALS((it)->GetIndex(), 0U);
            TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(0), 47u);
            TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(1), 31u);
            TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(2), 131u);
            TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(3), 51u);

            Node<3> *mesh_node = mesh.GetNodeOrHaloNode(31);
            Node<3> *iterator_node = (it)->GetNode(1);
            TS_ASSERT_EQUALS(iterator_node, mesh_node);
        }

#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    /**
     * Check that we can build a MixedDimensionMesh using the VTK mesh reader.
     */
    void TestBuild2DFromVtkMeshReader(void)
    {
#ifdef CHASTE_VTK
    VtkMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements.vtu");

    TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 40u);

    TetrahedralMesh<2,2> mesh;
    mesh.ConstructFromMeshReader(mesh_reader);

    TS_ASSERT_EQUALS(mesh.GetNumNodes(), 121u);
    TS_ASSERT_EQUALS(mesh.GetNumElements(), 200u);
    TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 40u);

    std::vector<c_vector<double, 2> > centroid;
    mesh_reader.GetCellData("Centroid", centroid);

    TS_ASSERT_EQUALS(centroid.size(), 200u);
    TS_ASSERT_DELTA(centroid[0](0), 0.0033, 1e-4);
    TS_ASSERT_DELTA(centroid[0](1), 0.0033, 1e-4);
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }
    /**
     * Check that we can build a 2D MixedDimensionMesh using the VTK mesh reader.
     */
    void TestBuildMixedMesh2DFromVtkMeshReader(void)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<2,2> mesh_reader("mesh/test/data/mixed_dimension_meshes/mixed_mesh_2d.vtu");
        TS_ASSERT_EQUALS(mesh_reader.GetNumCableElements(), 10u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumCableElementAttributes(), 1u);

        MixedDimensionMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 121u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 200u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 40u);
        TS_ASSERT_EQUALS(mesh.GetNumCableElements(), 10u);

        // Cover exception - note that ConstructFromMeshReader has reset the reader
        for (unsigned i=0; i<10; i++)
        {
            ElementData element_data = mesh_reader.GetNextCableElementData();
            TS_ASSERT_EQUALS(element_data.AttributeValue, i + 1.5);
        }
        TS_ASSERT_THROWS_THIS(mesh_reader.GetNextCableElementData(), "Trying to read data for a cable element that doesn't exist");

#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }
    /**
     * Check that we can build a 3D MixedDimensionMesh using the VTK mesh reader.
     */
    void TestBuildMixedMeshTetrahedralMeshFromVtkMeshReader(void)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/mixed_dimension_meshes/mixed_mesh_3d.vtu");
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 616u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumCableElements(), 5u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumCableElementAttributes(), 1u);

        MixedDimensionMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes and elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 387u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1298u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 616u);
        TS_ASSERT_EQUALS(mesh.GetNumCableElements(), 5u);

        for (MixedDimensionMesh<3,3>::CableElementIterator iter = mesh.GetCableElementIteratorBegin();
             iter != mesh.GetCableElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT_DELTA((*iter)->GetAttribute(), 2.0, 1e-9);
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    /**
     * Check that we can build a 2D in 3D mesh using the VTK mesh reader.
     */
    void TestLoadingSurfaceMeshFromVtkMeshReader(void)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<2,3> mesh_reader("mesh/test/data/cylinder.vtu");
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 1632u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 32u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumCableElements(), 0u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumCableElementAttributes(), 0u);

        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 832u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1632u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 32u);
        TS_ASSERT_EQUALS(mesh.GetNumCableElements(), 0u);
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }


    /**
     * Check that we can build a 1D in 3D Mesh using the VTK mesh reader.
     */
    void TestLoading1Din3DMeshFromVtkMeshReader(void)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<1,3> mesh_reader("mesh/test/data/branched_1d_in_3d_mesh.vtu");
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 30u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 31u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumEdges(), 3u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumCableElements(), 0u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumCableElementAttributes(), 0u);

        TetrahedralMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumEdges());
        TS_ASSERT_EQUALS(mesh.GetNumCableElements(), 0u);
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }
};

#endif /*TESTVTKMESHREADER_*/
