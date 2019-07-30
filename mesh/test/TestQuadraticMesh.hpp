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
#ifndef _TESTQUADRATICMESH_HPP_
#define _TESTQUADRATICMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "QuadraticMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"
#include "ArchiveOpener.hpp"
#include "Warnings.hpp"
#include "PetscMatTools.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestQuadraticMesh : public CxxTest::TestSuite
{
public:

    void TestQuadraticMesh1d()
    {
        QuadraticMesh<1> mesh;
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10u);

        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 11u);

        // Node 2 (ie middle) of element 0
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(2), 11u);
        TS_ASSERT_DELTA(mesh.GetNode(11)->rGetLocation()[0], 0.05, 1e-12);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            bool is_boundary_node = mesh.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node, ((x==0)||(x==1)));
        }
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            // Check internal nodes have corrent element associated with them
            std::set<unsigned> internal_node_elems;
            internal_node_elems.insert(mesh.GetElement(i)->GetIndex());
            TS_ASSERT_EQUALS(internal_node_elems,mesh.GetElement(i)->GetNode(2)->rGetContainingElementIndices());
        }
    }

    // Identical to above, except mesh is generated not read
    void TestQuadraticMesh1dAutomaticallyGenerated()
    {
        QuadraticMesh<1> mesh(0.1, 1.0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10u);

        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 11u);

        // Node 2 (ie middle) of element 0
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(2), 11u);
        TS_ASSERT_DELTA(mesh.GetNode(11)->rGetLocation()[0], 0.05, 1e-12);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            bool is_boundary_node = mesh.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node, ((x==0)||(x==1)));
        }
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            // Check internal nodes have correct element associated with them
            std::set<unsigned> internal_node_elems;
            internal_node_elems.insert(mesh.GetElement(i)->GetIndex());
            TS_ASSERT_EQUALS(internal_node_elems,mesh.GetElement(i)->GetNode(2)->rGetContainingElementIndices());
        }
    }

    void TestQuadraticMesh2d()
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1, false);
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 289u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 128u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 81u);

        // Each element should have 6 nodes
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetNumNodes(), 6u);

            for (unsigned j=0; j<2; j++)
            {
                // Check internal nodes have corrent element associated with them
                TS_ASSERT(mesh.GetElement(i)->GetNode(j+3)->GetNumContainingElements() <= 2u);
                TS_ASSERT(mesh.GetElement(i)->GetNode(j+3)->GetNumContainingElements() > 0u);

                std::set<unsigned> current_node_indices = mesh.GetElement(i)->GetNode(j)->rGetContainingElementIndices();
                TS_ASSERT_EQUALS(current_node_indices.count(mesh.GetElement(i)->GetIndex()), 1u);

                current_node_indices = mesh.GetElement(i)->GetNode(j+3)->rGetContainingElementIndices();
                TS_ASSERT_EQUALS(current_node_indices.count(mesh.GetElement(i)->GetIndex()), 1u);
            }
        }

        // Node 3 (ie fourth) of element 0
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(3), 82u);
        // Node 4 (ie fifth) of element 0
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(4), 83u);
        // Node 5 (ie last) of element 0
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(5), 81u);

        // Each boundary element should have three nodes
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            TS_ASSERT_EQUALS((*iter)->GetNumNodes(), 3u);
        }

        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();

        // The first edge has nodes 53 and 0, according to the edge file...
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(0), 53u);
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(1), 0u);
        // ...the midnode has to be computed (found) by the QuadraticMesh class
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(2), 81u);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            bool is_boundary_node = mesh.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node, ((x==0)||(x==1)||(y==0)||(y==1)));
        }
    }

    void TestQuadraticMesh3d()
    {
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element_quadratic",2,1, false);
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);

        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 4u);

        // Check getting global numbers of nodes 4-9 (in non-vertices)
        for (unsigned i=4; i<10; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(i), i);
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetNode(i)->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(mesh.GetNode(i)->GetNumContainingElements(), 1u);
        }

        // Lots of internal and boundary nodes in this mesh..
        QuadraticMesh<3> mesh2;
        TrianglesMeshReader<3,3> mesh_reader2("mesh/test/data/cube_1626_elements_quadratic",2,1, false);
        mesh2.ConstructFromMeshReader(mesh_reader2);

        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), 2570u);
        TS_ASSERT_EQUALS(mesh2.GetNumElements(), 1626u);
        TS_ASSERT_EQUALS(mesh2.GetNumVertices(), 375u);

        // Each element should have 10 nodes
        for (unsigned i=0; i<mesh2.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh2.GetElement(i)->GetNumNodes(), 10u);

            for (unsigned j=3; j<9; j++)
            {
                // Check internal nodes have corrent element associated with them
                TS_ASSERT(mesh2.GetElement(i)->GetNode(j)->GetNumContainingElements() > 0u);
                std::set<unsigned> current_node_indices = mesh2.GetElement(i)->GetNode(j)->rGetContainingElementIndices();
                TS_ASSERT_EQUALS(current_node_indices.count(mesh2.GetElement(i)->GetIndex()),1u);
            }
        }

        // Each boundary element should have 6 nodes
        for (TetrahedralMesh<3,3>::BoundaryElementIterator iter= mesh2.GetBoundaryElementIteratorBegin();
             iter != mesh2.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT_EQUALS((*iter)->GetNumNodes(), 6u);
        }

        TetrahedralMesh<3,3>::BoundaryElementIterator iter = mesh2.GetBoundaryElementIteratorBegin();

        // The first boundary elem has these nodes, according to the edge file..
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(0), 177u);
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(1), 43u);
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(2), 85u);
        // .. the internal nodes have to be computed (found) by the QuadraticMesh.
        // The nodes 177,43,85 are all in the third element in the ele file, and
        // they are nodes 1,3,2 respectively. Therefore, the internals are the local
        // nodes 9,5,8 respectively (look the the ordering picture), so..
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(3), 392u);
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(4), 388u);
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(5), 391u);

        for (unsigned i=0; i<mesh2.GetNumNodes(); i++)
        {
            double x = mesh2.GetNode(i)->rGetLocation()[0];
            double y = mesh2.GetNode(i)->rGetLocation()[1];
            double z = mesh2.GetNode(i)->rGetLocation()[2];
            bool is_boundary_node = mesh2.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node,  ((x==0)||(x==1)||(y==0)||(y==1)||(z==0)||(z==1)));
        }
    }

    void TestAutomaticallyGenerated2dMesh1()
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);

        TS_ASSERT_THROWS_CONTAINS(QuadraticMesh<2> bad_mesh(0.645, 1.0, 1.0), "does not divide");

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 8u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 4u);

        // Each element should have 6 nodes and a valid Jacobian
        double det;
        c_matrix<double, 2, 2> jacob;
        c_matrix<double, 2, 2> inv;

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetNumNodes(), 6u);
            mesh.GetInverseJacobianForElement(i, jacob, det, inv);
            TS_ASSERT_EQUALS(det, 1.0);
        }

        // Test vertex containment
        TS_ASSERT_EQUALS(mesh.GetNode(0)->GetNumContainingElements(), 1u); //(0,0)
        TS_ASSERT_EQUALS(mesh.GetNode(1)->GetNumContainingElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->GetNumContainingElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNode(3)->GetNumContainingElements(), 1u); //(1,1)

        // Test internal node containment
        TS_ASSERT_EQUALS(mesh.GetNode(4)->GetNumContainingElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNode(5)->GetNumContainingElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNode(6)->GetNumContainingElements(), 2u); //(.5,.5)
        TS_ASSERT_EQUALS(mesh.GetNode(7)->GetNumContainingElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNode(8)->GetNumContainingElements(), 1u);

        TS_ASSERT_DELTA( mesh.GetNode(3)->rGetLocation()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA( mesh.GetNode(3)->rGetLocation()[1], 1.0, 1e-6);

        // Test boundary elements
        unsigned num_boundary_elements=0;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT_EQUALS((*iter)->GetNumNodes(), 3u);

            bool all_x_zero =     (fabs((*iter)->GetNode(0)->rGetLocation()[0])<1e-6)
                               && (fabs((*iter)->GetNode(1)->rGetLocation()[0])<1e-6)
                               && (fabs((*iter)->GetNode(2)->rGetLocation()[0])<1e-6);

            bool all_x_one  =     (fabs((*iter)->GetNode(0)->rGetLocation()[0] - 1.0)<1e-6)
                               && (fabs((*iter)->GetNode(1)->rGetLocation()[0] - 1.0)<1e-6)
                               && (fabs((*iter)->GetNode(2)->rGetLocation()[0] - 1.0)<1e-6);

            bool all_y_zero =     (fabs((*iter)->GetNode(0)->rGetLocation()[1])<1e-6)
                               && (fabs((*iter)->GetNode(1)->rGetLocation()[1])<1e-6)
                               && (fabs((*iter)->GetNode(2)->rGetLocation()[1])<1e-6);

            bool all_y_one  =     (fabs((*iter)->GetNode(0)->rGetLocation()[1] - 1.0)<1e-6)
                               && (fabs((*iter)->GetNode(1)->rGetLocation()[1] - 1.0)<1e-6)
                               && (fabs((*iter)->GetNode(2)->rGetLocation()[1] - 1.0)<1e-6);

            TS_ASSERT_EQUALS(true, all_x_zero || all_x_one || all_y_zero || all_y_one);
            num_boundary_elements++;
        }
        TS_ASSERT_EQUALS(num_boundary_elements, 4u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4u);

    }

    void TestAutomaticallyGenerated2dMesh2()
    {
        QuadraticMesh<2> mesh(3.14159/10,  3.14159, 3.14159/2);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21*11u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 60u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 100u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 11*6u);

        // Each element should have 6 nodes
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetNumNodes(), 6u);
        }
        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetBoundaryElement(i)->GetNumNodes(), 3u);
        }

        TS_ASSERT_DELTA( mesh.GetNode(65)->rGetLocation()[0], 3.14159, 1e-4);
        TS_ASSERT_DELTA( mesh.GetNode(65)->rGetLocation()[1], 3.14159/2, 1e-5);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 30u);
    }

    void TestAutomaticallyGenerated3dMeshSimple()
    {
        double h = 3.14159;
        double width = h;

        QuadraticMesh<3> mesh(h, width, 2*width, 3*width);

        TS_ASSERT_THROWS_CONTAINS(QuadraticMesh<3> bad_mesh(0.645, 1.0, 1.0, 1.0), "does not divide");

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 3*5*7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6*1*2*3u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 2*3*4u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 90u);

        for (unsigned i=1; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> x = mesh.GetNode(i)->rGetLocation();

            // Check the extra nodes aren't (0,0,0).
            // This fails with 32bit outdated binary.
            TS_ASSERT_LESS_THAN(1e-12, norm_2(x)); // assert x not equal to 0
        }

        TS_ASSERT_DELTA( mesh.GetNode(23)->rGetLocation()[0], width, 1e-8);
        TS_ASSERT_DELTA( mesh.GetNode(23)->rGetLocation()[1], 2*width, 1e-8);
        TS_ASSERT_DELTA( mesh.GetNode(23)->rGetLocation()[2], 3*width, 1e-8);

        // Second 1 by 1 by 1 mesh
        QuadraticMesh<3> mesh2(h, width, width, width);

        for (unsigned i=1; i<mesh2.GetNumNodes(); i++)
        {
            //Check that all nodes have containg elements
            TS_ASSERT_LESS_THAN(0u, mesh2.GetNode(i)->GetNumContainingElements());
            //Mid-point of cube will have access to all 6 elements
            TS_ASSERT_LESS_THAN_EQUALS(mesh2.GetNode(i)->GetNumContainingElements(), 6u);
        }

        TS_ASSERT_EQUALS(mesh2.CalculateMaximumContainingElementsPerProcess(), 6U); //The midpoint, as given above
        // There are  8 vertex nodes in the cube
        // There are 12 internal nodes on the cube edges
        // There are  6 internal nodes on the diagonals to the cube faces
        // There is   1 interal node on the separating diagonal
        // 8V + 19I = 27 nodes
        TS_ASSERT_EQUALS(mesh2.CalculateMaximumNodeConnectivityPerProcess(), 27U); //The midpoint, as given above
    }

    void TestAutomaticallyGenerated3dMesh()
    {
        QuadraticMesh<3> mesh(0.5,  2.5, 2.5, 2.5);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 1331u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 750u); // 5 cubes in each direction = 125 cubes => 125 x 6 tetrahedra per cube = 750
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 216u); // 6^3 = 216
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 602u);

        // Each element should have 10 nodes
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetNumNodes(), 10u);
        }

        TS_ASSERT_DELTA(mesh.GetNode(215)->rGetLocation()[0], 2.5, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(215)->rGetLocation()[1], 2.5, 1e-5);
        TS_ASSERT_DELTA(mesh.GetNode(215)->rGetLocation()[2], 2.5, 1e-4);
        TS_ASSERT_EQUALS(mesh.CalculateMaximumContainingElementsPerProcess(), 24U); // Four surrounding cubes may have all 6 tetrahedra meeting at a node
        TS_ASSERT_EQUALS(mesh.CalculateMaximumNodeConnectivityPerProcess(), 65U);
    }

    void TestWritingReadingBoundaryElementsWithContainingElementInfo()
    {
        // This mesh has quadratic node and ele files, a linear face file that has containing element info
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_152_elements_v3",2,1,true);
        mesh.ConstructFromMeshReader(mesh_reader);

        for (QuadraticMesh<3>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT_EQUALS( (*iter)->GetNumNodes(), 6u );
        }
    }

    void TestArchiving()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "quadratic_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("quadratic_mesh");

        AbstractTetrahedralMesh<3,3>* const p_mesh = new QuadraticMesh<3>;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements_fully_quadratic", 2, 2, false);
        static_cast<QuadraticMesh<3>*>(p_mesh)->ConstructFromMeshReader(mesh_reader);

        {
            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) << p_mesh;
        }

        {
            // Should archive the most abstract class you can to check boost knows what individual classes are.
            // (but here AbstractMesh doesn't have the methods below).
            AbstractTetrahedralMesh<3,3>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_mesh2;

            // Compare the boundary elements of both meshes, should be identical (as one was created from the other)
            QuadraticMesh<3>::BoundaryElementIterator iter1 = p_mesh->GetBoundaryElementIteratorBegin();

            for (QuadraticMesh<3>::BoundaryElementIterator iter2 = p_mesh2->GetBoundaryElementIteratorBegin();
                 iter2 != p_mesh2->GetBoundaryElementIteratorEnd();
                 ++iter2)
            {
                TS_ASSERT_EQUALS( (*iter1)->GetNumNodes(), 6u );
                TS_ASSERT_EQUALS( (*iter2)->GetNumNodes(), 6u );

                for (unsigned i=0; i<6; i++)
                {
                   TS_ASSERT_EQUALS( (*iter1)->GetNodeGlobalIndex(i), (*iter2)->GetNodeGlobalIndex(i));
                }
                iter1++;
            }

            delete p_mesh2;
        }
        delete p_mesh;
    }

    void TestConstructRegularSlabMesh_Directly_1d()
    {
        QuadraticMesh<1> mesh;
        TS_ASSERT_THROWS_THIS(mesh.ConstructRegularSlabMesh(0.75, 1.0), "Space step does not divide the size of the mesh");
        mesh.ConstructRegularSlabMesh(0.1, 1.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 11u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2u);

        for (unsigned i=0; i<mesh.GetNumVertices(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], (i+0.0)/10, 1e-8);

            bool is_boundary = (i==0 || i+1==mesh.GetNumVertices());
            TS_ASSERT_EQUALS(mesh.GetNode(i)->IsBoundaryNode(), is_boundary);

            std::set<unsigned> containing_elems = mesh.GetNode(i)->rGetContainingElementIndices();
            TS_ASSERT_EQUALS(containing_elems.size(), (is_boundary ? 1u : 2u));
        }

        for (unsigned i=mesh.GetNumVertices(); i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], (i-11.0)/10 + 0.05, 1e-8);

            std::set<unsigned> containing_elems = mesh.GetNode(i)->rGetContainingElementIndices();
            TS_ASSERT_EQUALS(containing_elems.size(), 1u);
        }

        TrianglesMeshWriter<1,1> mesh_writer("TestQuadraticMesh", "QuadraticSlab1D", false);
        mesh_writer.WriteFilesUsingMesh(mesh);

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<1,1> quadratic_mesh_reader(output_dir + "QuadraticSlab1D", 2, 2);

        QuadraticMesh<1> quad_mesh_read_back;
        quad_mesh_read_back.ConstructFromMeshReader(quadratic_mesh_reader);

        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumNodes(), 21u);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumVertices(), 11u);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumBoundaryNodes(), 2u);
    }


    void TestConstructRegularSlabMesh_Directly_2d()
    {
        QuadraticMesh<2> mesh;
        mesh.ConstructRegularSlabMesh(0.1, 1.0, 2.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21*41u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 11*21u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2*10*20u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 120u);


        TrianglesMeshWriter<2,2> mesh_writer("TestQuadraticMesh", "QuadraticSlab2D", false);
        mesh_writer.WriteFilesUsingMesh(mesh);

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,2> quadratic_mesh_reader(output_dir + "QuadraticSlab2D", 2, 2);

        QuadraticMesh<2> quad_mesh_read_back;
        quad_mesh_read_back.ConstructFromMeshReader(quadratic_mesh_reader);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumNodes(), 21*41u);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumVertices(), 11*21u);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumElements(), 2*10*20u);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumBoundaryNodes(), 120u);
    }

    void TestConstructRegularSlabMesh_Directly_3d()
    {
        QuadraticMesh<3> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0, 2.0, 3.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 3*5*7u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 2*3*4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6*1*2*3u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 90u);

        TrianglesMeshWriter<3,3> mesh_writer("TestQuadraticMesh", "QuadraticSlab3D", false);
        mesh_writer.WriteFilesUsingMesh(mesh);

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<3,3> quadratic_mesh_reader(output_dir + "QuadraticSlab3D", 2, 2);

        QuadraticMesh<3> quad_mesh_read_back;
        quad_mesh_read_back.ConstructFromMeshReader(quadratic_mesh_reader);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumNodes(), 3*5*7u);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumVertices(), 2*3*4u);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumElements(), 6*1*2*3u);
        TS_ASSERT_EQUALS(quad_mesh_read_back.GetNumBoundaryNodes(), 90u);

    }

    void TestConstructionConversionVersusConstruction2dNoStagger()
    {
        QuadraticMesh<2> quad_mesh_read_back;
        QuadraticMesh<2> quad_mesh_constructed;
        unsigned width  = 1.0;
        unsigned height = 2.0;
        {
            //Two-dimensional two squares
            TetrahedralMesh<2,2> mesh;
            bool stagger=false;
            mesh.ConstructRectangularMesh(width, height, stagger);
            TrianglesMeshWriter<2,2> mesh_writer("TestQuadraticMesh", "TempGrid2d", false);
            mesh_writer.WriteFilesUsingMesh(mesh);

            //Convert to quadratic
            std::string output_dir = mesh_writer.GetOutputDirectory();
            TrianglesMeshReader<2,2> quadratic_mesh_reader(output_dir + "TempGrid2d");
            quad_mesh_read_back.ConstructFromLinearMeshReader(quadratic_mesh_reader);

            quad_mesh_constructed.ConstructRectangularMesh(width, height, stagger);
        }
        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumNodes(), quad_mesh_read_back.GetNumNodes());
        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumBoundaryNodes(), quad_mesh_read_back.GetNumBoundaryNodes());
        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumVertices(), quad_mesh_read_back.GetNumVertices());

        for (unsigned elem=0; elem<quad_mesh_constructed.GetNumElements(); elem++)
        {
            Element<2,2>* p_elem_constructed =  quad_mesh_constructed.GetElement(elem);
            Element<2,2>* p_elem_read_back =  quad_mesh_read_back.GetElement(elem);
            TS_ASSERT_EQUALS(p_elem_constructed->GetNumNodes(), p_elem_read_back->GetNumNodes());
            for (unsigned i = 0; i < p_elem_constructed->GetNumNodes(); i++)
            {
                c_vector<double, 2> loc_read_back;
                loc_read_back = p_elem_read_back->GetNode(i)->rGetLocation();

                c_vector<double, 2> loc_constructed;
                loc_constructed = p_elem_constructed->GetNode(i)->rGetLocation();

                TS_ASSERT_DELTA(loc_read_back[0], loc_constructed[0], 1e-10);
                TS_ASSERT_DELTA(loc_read_back[1], loc_constructed[1], 1e-10);

                TS_ASSERT_EQUALS(p_elem_read_back->GetNode(i)->IsBoundaryNode(), p_elem_constructed->GetNode(i)->IsBoundaryNode());
            }
        }
        // Can't check edges exactly because the linear to quadratic converter doesn't send any boundary information to the
        // external mesher.  (So the edges come back in a different order.)
        for (unsigned b_elem=0; b_elem<quad_mesh_constructed.GetNumBoundaryElements(); b_elem++)
        {
            BoundaryElement<1,2>* p_b_elem_constructed =  quad_mesh_constructed.GetBoundaryElement(b_elem);
            BoundaryElement<1,2>* p_b_elem_read_back =  quad_mesh_read_back.GetBoundaryElement(b_elem);
            TS_ASSERT_EQUALS(p_b_elem_constructed->GetNumNodes(), p_b_elem_read_back->GetNumNodes());
        }
    }

    void TestConstructionConversionVersusConstruction2dWithStagger()
    {
        QuadraticMesh<2> quad_mesh_read_back;
        QuadraticMesh<2> quad_mesh_constructed;
        double width  = 1.0;
        double height = 2.0;
        {
            //Two-dimensional two squares
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(1.0, width, height);
            TrianglesMeshWriter<2,2> mesh_writer("TestQuadraticMesh", "TempGrid2d", false);
            mesh_writer.WriteFilesUsingMesh(mesh);

            //Convert to quadratic
            std::string output_dir = mesh_writer.GetOutputDirectory();
            TrianglesMeshReader<2,2> quadratic_mesh_reader(output_dir + "TempGrid2d");
            quad_mesh_read_back.ConstructFromLinearMeshReader(quadratic_mesh_reader);

            quad_mesh_constructed.ConstructRegularSlabMesh(1.0, width, height);
        }
        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumNodes(), quad_mesh_read_back.GetNumNodes());
        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumBoundaryNodes(), quad_mesh_read_back.GetNumBoundaryNodes());
        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumVertices(), quad_mesh_read_back.GetNumVertices());

        for (unsigned elem=0; elem<quad_mesh_constructed.GetNumElements(); elem++)
        {
            Element<2,2>* p_elem_constructed =  quad_mesh_constructed.GetElement(elem);
            Element<2,2>* p_elem_read_back =  quad_mesh_read_back.GetElement(elem);
            TS_ASSERT_EQUALS(p_elem_constructed->GetNumNodes(), p_elem_read_back->GetNumNodes());
            for (unsigned i = 0; i < p_elem_constructed->GetNumNodes(); i++)
            {
                c_vector<double, 2> loc_read_back;
                loc_read_back = p_elem_read_back->GetNode(i)->rGetLocation();
                c_vector<double, 2> loc_constructed;
                loc_constructed = p_elem_constructed->GetNode(i)->rGetLocation();
                TS_ASSERT_DELTA(loc_read_back[0], loc_constructed[0], 1e-10);
                TS_ASSERT_DELTA(loc_read_back[1], loc_constructed[1], 1e-10);

                TS_ASSERT_EQUALS(p_elem_read_back->GetNode(i)->IsBoundaryNode(), p_elem_constructed->GetNode(i)->IsBoundaryNode());
            }
        }
        // Can't check edges exactly because the linear to quadratic converter doesn't send any boundary information to the
        // external mesher.  (So the edges come back in a different order.)
        for (unsigned b_elem=0; b_elem<quad_mesh_constructed.GetNumBoundaryElements(); b_elem++)
        {
            BoundaryElement<1,2>* p_b_elem_constructed =  quad_mesh_constructed.GetBoundaryElement(b_elem);
            BoundaryElement<1,2>* p_b_elem_read_back =  quad_mesh_read_back.GetBoundaryElement(b_elem);
            TS_ASSERT_EQUALS(p_b_elem_constructed->GetNumNodes(), p_b_elem_read_back->GetNumNodes());
        }
    }

    void TestConstructionConversionVersusConstruction3d()
    {
        QuadraticMesh<3> quad_mesh_read_back;
        QuadraticMesh<3> quad_mesh_constructed;
        double width  = 1.0;
        double height = 2.0;
        double depth  = 3.0;
        {
            //Three-dimensional cubes
            TetrahedralMesh<3,3> mesh;
            mesh.ConstructRegularSlabMesh(1.0, width, height, depth);
            TrianglesMeshWriter<3,3> mesh_writer("TestQuadraticMesh", "TempGrid3d", false);
            mesh_writer.WriteFilesUsingMesh(mesh);

            //Convert to quadratic
            std::string output_dir = mesh_writer.GetOutputDirectory();
            TrianglesMeshReader<3,3> quadratic_mesh_reader(output_dir + "TempGrid3d");
            quad_mesh_read_back.ConstructFromLinearMeshReader(quadratic_mesh_reader);

            quad_mesh_constructed.ConstructRegularSlabMesh(1.0, width, height, depth);
        }

        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumNodes(), quad_mesh_read_back.GetNumNodes());
        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumBoundaryNodes(), quad_mesh_read_back.GetNumBoundaryNodes());
        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumVertices(), quad_mesh_read_back.GetNumVertices());
        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumElements(), quad_mesh_read_back.GetNumElements());
        TS_ASSERT_EQUALS(quad_mesh_constructed.GetNumBoundaryElements(), quad_mesh_read_back.GetNumBoundaryElements());

        for (unsigned elem=0; elem<quad_mesh_constructed.GetNumElements(); elem++)
        {
            Element<3,3>* p_elem_constructed =  quad_mesh_constructed.GetElement(elem);
            Element<3,3>* p_elem_read_back =  quad_mesh_read_back.GetElement(elem);
            TS_ASSERT_EQUALS(p_elem_constructed->GetNumNodes(), p_elem_read_back->GetNumNodes());
            for (unsigned i = 0; i < p_elem_constructed->GetNumNodes(); i++)
            {
                c_vector<double, 3> loc_read_back;
                loc_read_back = p_elem_read_back->GetNode(i)->rGetLocation();

                c_vector<double, 3> loc_constructed;
                loc_constructed = p_elem_constructed->GetNode(i)->rGetLocation();

                TS_ASSERT_DELTA(loc_read_back[0], loc_constructed[0], 1e-10);
                TS_ASSERT_DELTA(loc_read_back[1], loc_constructed[1], 1e-10);
                TS_ASSERT_DELTA(loc_read_back[2], loc_constructed[2], 1e-10);

                TS_ASSERT_EQUALS(p_elem_read_back->GetNode(i)->IsBoundaryNode(), p_elem_constructed->GetNode(i)->IsBoundaryNode());
            }
        }

        for (unsigned b_elem=0; b_elem<quad_mesh_constructed.GetNumBoundaryElements(); b_elem++)
        {
            BoundaryElement<2,3>* p_b_elem_constructed =  quad_mesh_constructed.GetBoundaryElement(b_elem);
            BoundaryElement<2,3>* p_b_elem_read_back =  quad_mesh_read_back.GetBoundaryElement(b_elem);
            TS_ASSERT_EQUALS(p_b_elem_constructed->GetNumNodes(), p_b_elem_read_back->GetNumNodes());
        }
    }

    void TestLinearToQuadraticMeshConversion2d()
    {
        QuadraticMesh<2> quad_mesh;
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements");
        quad_mesh.ConstructFromLinearMeshReader(reader);

        TS_ASSERT_EQUALS(quad_mesh.GetNumNodes(), 289u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumVertices(), 81u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumElements(), 128u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumBoundaryNodes(), 64u);

        // Output
        TrianglesMeshWriter<2,2> mesh_writer("TestQuadraticMesh", "converted_square", false);
        mesh_writer.WriteFilesUsingMesh(quad_mesh);  //Compare with

        // Read in the new mesh to check it worked
        OutputFileHandler handler("TestQuadraticMesh",false);
        std::string full_new_mesh = handler.GetOutputDirectoryFullPath() + "converted_square";

        // Input again
        QuadraticMesh<2> quad_mesh_after_conversion;
        TrianglesMeshReader<2,2> quad_reader(full_new_mesh, 2, 2);
        quad_mesh_after_conversion.ConstructFromMeshReader(quad_reader);

        TS_ASSERT_EQUALS(quad_mesh_after_conversion.GetNumNodes(), 17*17u);
        TS_ASSERT_EQUALS(quad_mesh_after_conversion.GetNumVertices(), 81u);
        TS_ASSERT_EQUALS(quad_mesh_after_conversion.GetNumElements(), 128u);
        TS_ASSERT_EQUALS(quad_mesh_after_conversion.GetNumBoundaryNodes(), 64u);
    }

    void TestLinearToQuadraticMeshConversion2dNonconvex()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/l_shape");

        TetrahedralMesh<2,2> linear_mesh;
        linear_mesh.ConstructFromMeshReader(reader);

        TS_ASSERT_EQUALS(linear_mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(linear_mesh.GetNumElements(), 6u);
        TS_ASSERT_EQUALS(linear_mesh.GetNumBoundaryNodes(), 8u);
        TS_ASSERT_EQUALS(linear_mesh.GetVolume(), 3.0);
        TS_ASSERT_EQUALS(linear_mesh.GetSurfaceArea(), 8.0);

        reader.Reset();
        QuadraticMesh<2> quad_mesh;
        quad_mesh.ConstructFromLinearMeshReader(reader);
        TS_ASSERT_EQUALS(quad_mesh.GetNumVertices(), 8u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumElements(), 6u);
        TS_ASSERT_EQUALS(quad_mesh.GetVolume(), 3.0);
        TS_ASSERT_EQUALS(quad_mesh.GetSurfaceArea(), 8.0);
    }

    /* HOW_TO_TAG Mesh
     * Convert a linear tetrahedral mesh to quadratic and write back to file.
     */
    void TestLinearToQuadraticMeshConversion3d()
    {
        QuadraticMesh<3> quad_mesh;
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_136_elements");
        quad_mesh.ConstructFromLinearMeshReader(reader);


        TS_ASSERT_EQUALS(quad_mesh.GetNumNodes(), 285u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumVertices(), 51u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumElements(), 136u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumBoundaryNodes(), 194u);

        //Output
        TrianglesMeshWriter<3,3> mesh_writer("TestQuadraticMesh", "converted_cube", false);
        mesh_writer.WriteFilesUsingMesh(quad_mesh);

        // read in the new mesh to check it worked
        OutputFileHandler handler("TestQuadraticMesh", false);
        std::string full_new_mesh = handler.GetOutputDirectoryFullPath() + "converted_cube";

        //Input again
        QuadraticMesh<3> quad_mesh_after_conversion;
        TrianglesMeshReader<3,3> quad_reader(full_new_mesh, 2, 2);
        quad_mesh_after_conversion.ConstructFromMeshReader(quad_reader);

        TS_ASSERT_EQUALS(quad_mesh_after_conversion.GetNumNodes(), 285u);
        TS_ASSERT_EQUALS(quad_mesh_after_conversion.GetNumVertices(), 51u);
        TS_ASSERT_EQUALS(quad_mesh_after_conversion.GetNumElements(), 136u);
        TS_ASSERT_EQUALS(quad_mesh_after_conversion.GetNumBoundaryNodes(), 194u);
    }

    void TestLinearToQuadraticMeshConversion3dNonconvex()
    {
        TrianglesMeshReader<3,3> reader("mesh/test/data/l_shape3d");

        TetrahedralMesh<3,3> linear_mesh;
        linear_mesh.ConstructFromMeshReader(reader);

        TS_ASSERT_EQUALS(linear_mesh.GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(linear_mesh.GetNumElements(), 18u); // 3 cubes
        TS_ASSERT_EQUALS(linear_mesh.GetNumBoundaryNodes(), 16u); // All on boundary
        TS_ASSERT_DELTA(linear_mesh.GetVolume(), 3.0, 1e-15);  // 3 cubes
        TS_ASSERT_DELTA(linear_mesh.GetSurfaceArea(), 14.0, 1e-15);

        reader.Reset();
        QuadraticMesh<3> quad_mesh;
        quad_mesh.ConstructFromLinearMeshReader(reader);
        TS_ASSERT_EQUALS(quad_mesh.GetNumVertices(), 16u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumElements(), 18u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumNodes(), 63u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumBoundaryNodes(), 58u); // All on boundary
        TS_ASSERT_DELTA(quad_mesh.GetVolume(), 3.0, 1e-15);
        TS_ASSERT_DELTA(quad_mesh.GetSurfaceArea(), 14.0, 1e-15);
    }

    void TestQuadraticMesh2dReordered()
    {
        // Quadratics mesh - with different ordering
        QuadraticMesh<2> quad_mesh;
        TrianglesMeshReader<2,2> mesh_reader1("mesh/test/data/square_128_elements_quadratic_reordered",2,1,false);
        quad_mesh.ConstructFromMeshReader(mesh_reader1);


        // Linear mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), quad_mesh.GetNumVertices());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned quad_index=i;
            //Quad mesh has a minor permutation: vertex node 4 now appears at index 81
            if (i==4)
            {
                quad_index=81;
            }

            double lin_x = mesh.GetNode(i)->rGetLocation()[0];
            double lin_y = mesh.GetNode(i)->rGetLocation()[1];
            double quad_x = quad_mesh.GetNode(quad_index)->rGetLocation()[0];
            double quad_y = quad_mesh.GetNode(quad_index)->rGetLocation()[1];
            TS_ASSERT_DELTA(lin_x, quad_x, 1e-8);
            TS_ASSERT_DELTA(lin_y, quad_y, 1e-8);
        }
    }


    /**
     * Check that we can build a QuadraticMesh using the VTK mesh reader.
     */
    void TestBuildQuadraticMeshFromVtkMeshReader(void)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/heart_decimation.vtu");

        TetrahedralMesh<3,3> tet_mesh;
        tet_mesh.ConstructFromMeshReader(mesh_reader);
        mesh_reader.Reset();

        QuadraticMesh<3> quad_mesh;
        quad_mesh.ConstructFromLinearMeshReader(mesh_reader);


        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(tet_mesh.GetNumNodes(), 173u);
        TS_ASSERT_EQUALS(tet_mesh.GetNumElements(), 610u);
        TS_ASSERT_EQUALS(tet_mesh.GetNumBoundaryElements(), 312u);
        TS_ASSERT_EQUALS(tet_mesh.GetNumBoundaryNodes(), 158u);

        TS_ASSERT_EQUALS(quad_mesh.GetNumNodes(), 1110u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumVertices(), 173u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumElements(), 610u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumBoundaryElements(), 312u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumBoundaryNodes(), 626u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(tet_mesh.GetNode(0)->GetPoint()[0], 0.0963, 1e-4);
        TS_ASSERT_DELTA(tet_mesh.GetNode(0)->GetPoint()[1], 0.3593, 1e-4);
        TS_ASSERT_DELTA(tet_mesh.GetNode(0)->GetPoint()[2], 0.9925, 1e-4);

        TS_ASSERT_DELTA(tet_mesh.GetNode(8)->GetPoint()[0], 1.0969, 1e-4);
        TS_ASSERT_DELTA(tet_mesh.GetNode(8)->GetPoint()[1], 0.6678, 1e-4);
        TS_ASSERT_DELTA(tet_mesh.GetNode(8)->GetPoint()[2], 0.7250, 1e-4);

        TS_ASSERT_DELTA(quad_mesh.GetNode(0)->GetPoint()[0], 0.0963, 1e-4);
        TS_ASSERT_DELTA(quad_mesh.GetNode(0)->GetPoint()[1], 0.3593, 1e-4);
        TS_ASSERT_DELTA(quad_mesh.GetNode(0)->GetPoint()[2], 0.9925, 1e-4);

        TS_ASSERT_DELTA(quad_mesh.GetNode(8)->GetPoint()[0], 1.0969, 1e-4);
        TS_ASSERT_DELTA(quad_mesh.GetNode(8)->GetPoint()[1], 0.6678, 1e-4);
        TS_ASSERT_DELTA(quad_mesh.GetNode(8)->GetPoint()[2], 0.7250, 1e-4);

        //Use ordinary functionality - using the "wrong" method gives back a warning
        quad_mesh.Clear();
        mesh_reader.Reset();
        TS_ASSERT_EQUALS(quad_mesh.GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        quad_mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Reading a (linear) tetrahedral mesh and converting it to a QuadraticMesh.  This involves making an external library call to Triangle/Tetgen in order to compute internal nodes");
        TS_ASSERT_EQUALS(quad_mesh.GetNumNodes(), 1110u);

#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }


    void CalculateConnectivityMatrix(Mat& matrix, AbstractTetrahedralMesh<2,2>& rMesh)
    {
        PetscTools::SetupMat(matrix, rMesh.GetNumNodes(), rMesh.GetNumNodes(), 17);
        for (TetrahedralMesh<2,2>::ElementIterator iter
                    = rMesh.GetElementIteratorBegin();
                    iter != rMesh.GetElementIteratorEnd();
                    ++iter)
        {
            for (unsigned i=0; i<iter->GetNumNodes(); i++)
            {
                unsigned global_index1 = iter->GetNodeGlobalIndex(i);
                for (unsigned j=i+1; j<iter->GetNumNodes(); j++)
                {
                    unsigned global_index2 = iter->GetNodeGlobalIndex(j);
                    PetscMatTools::SetElement(matrix, global_index1, global_index2, 1.0);
                    PetscMatTools::SetElement(matrix, global_index2, global_index1, 1.0);
                }
            }
        }
        PetscMatTools::Finalise(matrix);
    }

    std::vector<unsigned> CalculateMatrixFill(AbstractTetrahedralMesh<2,2>& rMesh)
    {
        //Get some statistics about matrix fill
        Mat matrix;

        CalculateConnectivityMatrix(matrix, rMesh);

        std::vector<unsigned> upper_hist(rMesh.GetNumNodes(), 0.0);
        double error_sum = 0;

        PetscInt lo, hi;
        PetscMatTools::GetOwnershipRange(matrix, lo, hi);
        for (PetscInt row=lo; row<hi; row++)
        {
            PetscInt num_entries;
            const PetscInt* column_indices;
            const PetscScalar* values;
            MatGetRow(matrix, row, &num_entries, &column_indices, &values);
            for (PetscInt col=0; col<num_entries; col++)
            {
                if (column_indices[col] >= row)
                {
                    error_sum += (column_indices[col] - row)*(column_indices[col] - row);
                    upper_hist[ column_indices[col] - row]++;
                }
            }
            MatRestoreRow(matrix, row, &num_entries, &column_indices, &values);
        }
        PetscTools::Destroy(matrix);

        std::vector<unsigned> global_hist(rMesh.GetNumNodes());
        MPI_Allreduce( &upper_hist[0], &global_hist[0], rMesh.GetNumNodes(), MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);

        double global_error_sum = 0;
        MPI_Allreduce( &error_sum, &global_error_sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

        return global_hist;
    }

    void TestElementsContainedByNodes3d()
    {
        QuadraticMesh<3> mesh;
        double h = 1.0;
        mesh.ConstructRegularSlabMesh(h, 2.0, 1.0, 1.0);

        for (unsigned node_index = 0; node_index<mesh.GetNumNodes(); node_index++)
        {
            std::set<unsigned> elements = mesh.GetNode(node_index)->rGetContainingElementIndices();
            for (std::set<unsigned>::iterator iter = elements.begin(); iter != elements.end(); iter++)
            {
                Element<3,3>* p_element = mesh.GetElement(*iter);
                bool found_node = false;
                for (unsigned i=0; i<p_element->GetNumNodes(); i++)
                {
                    unsigned this_node = p_element->GetNodeGlobalIndex(i);
                    if (this_node == node_index)
                    {
                        found_node = true;
                    }
                }
                TS_ASSERT(found_node);
            }

            std::set<unsigned> boundary_elements = mesh.GetNode(node_index)->rGetContainingBoundaryElementIndices();
            for (std::set<unsigned>::iterator iter = boundary_elements.begin(); iter != boundary_elements.end(); iter++)
            {
                BoundaryElement<2,3>* p_element = mesh.GetBoundaryElement(*iter);
                bool found_node = false;
                for (unsigned i=0; i<p_element->GetNumNodes(); i++)
                {
                    unsigned this_node = p_element->GetNodeGlobalIndex(i);
                    if (this_node == node_index)
                    {
                        found_node = true;
                    }
                }
                TS_ASSERT(found_node);
            }
        }
    }
};

#endif // _TESTQUADRATICMESH_HPP_
