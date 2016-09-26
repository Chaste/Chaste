/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTVORONOIPRISM3DVERTEXMESHGENERATOR_HPP_
#define TESTVORONOIPRISM3DVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "VoronoiPrism3dVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"

class TestVoronoiPrism3dVertexMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestSimplestMesh() throw(Exception)
    {
        VoronoiPrism3dVertexMeshGenerator generator(10, 10, 1, 0);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        VertexMeshWriter<3, 3> vertex_mesh_writer("TestVoronoiPrism3dVertexMesh/Test0Relax", "10x10 0relax");

        std::vector<double> cell_ids;
        for (unsigned id=0 ; id < p_mesh->GetNumElements() ; id++)
        {
            cell_ids.push_back(double(id));
        }

        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);
        vertex_mesh_writer.WriteVtkUsingMesh(*p_mesh);

    }

    void TestSimplestMesh2Relaxation() throw(Exception)
    {
        VoronoiPrism3dVertexMeshGenerator generator(10, 10, 1, 2);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        VertexMeshWriter<3, 3> vertex_mesh_writer("TestVoronoiPrism3dVertexMesh/Test2Relax", "10x10 2relax");

        std::vector<double> cell_ids;
        for (unsigned id=0 ; id < p_mesh->GetNumElements() ; id++)
        {
            cell_ids.push_back(double(id));
        }

        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);
        vertex_mesh_writer.WriteVtkUsingMesh(*p_mesh);

    }

    void TestSimplestMesh6Relaxation() throw(Exception)
    {
        VoronoiPrism3dVertexMeshGenerator generator(10, 10, 1, 6);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        VertexMeshWriter<3, 3> vertex_mesh_writer("TestVoronoiPrism3dVertexMesh/Test6Relax", "10x10 6relax");

        std::vector<double> cell_ids;
        for (unsigned id=0 ; id < p_mesh->GetNumElements() ; id++)
        {
            cell_ids.push_back(double(id));
        }

        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);
        vertex_mesh_writer.WriteVtkUsingMesh(*p_mesh);

    }

    void TestSimplestMesh15Relaxation() throw(Exception)
    {
        VoronoiPrism3dVertexMeshGenerator generator(10, 10, 1, 15);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        VertexMeshWriter<3, 3> vertex_mesh_writer("TestVoronoiPrism3dVertexMesh/Test15Relax", "10x10 15relax");

        std::vector<double> cell_ids;
        for (unsigned id=0 ; id < p_mesh->GetNumElements() ; id++)
        {
            cell_ids.push_back(double(id));
        }

        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);
        vertex_mesh_writer.WriteVtkUsingMesh(*p_mesh);

    }

    void TestGenerateAnotherOne() throw(Exception)
    {
        VoronoiPrism3dVertexMeshGenerator generator(20, 15, 3, 1, 1);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        VertexMeshWriter<3, 3> vertex_mesh_writer("TestVoronoiPrism3dVertexMesh/BiggerTest", "20x15x3");

        std::vector<double> cell_ids;
        for (unsigned id=0 ; id < p_mesh->GetNumElements() ; id++)
        {
            cell_ids.push_back(double(id));
        }

        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);
        vertex_mesh_writer.WriteVtkUsingMesh(*p_mesh);
    }



    void TestSimpleMesh() throw(Exception)
    {

        // Generate a mesh that is 20 cells wide in x, 12 cells wide in y, 1 unit high in Z,
        // with 4 Lloyd's relaxation steps and target average element apical area 1.23
        VoronoiPrism3dVertexMeshGenerator generator(20, 12, 1, 4, 1.23);
        MutableVertexMesh<3, 3>* p_mesh_a = generator.GetMesh();
        TS_ASSERT_THROWS_THIS(generator.GetMeshAfterReMesh(),
                              "Remeshing has not been implemented in 3D (see #827 and #860)\n");

        TS_ASSERT_EQUALS(p_mesh_a->GetNumNodes(), 1096u);
        TS_ASSERT_EQUALS(p_mesh_a->GetNumElements(), 240u);

        // Check average cell area is correct
        double average_voloume = 0.0;
        for (unsigned elem_idx = 0 ; elem_idx < p_mesh_a->GetNumElements() ; elem_idx++)
        {
            average_voloume += p_mesh_a->GetVolumeOfElement(elem_idx);
        }
        average_voloume /= double(p_mesh_a->GetNumElements());

        TS_ASSERT_DELTA(average_voloume, 1.23, 1e-6);

        VertexMeshWriter<3, 3> vertex_mesh_writer("TestVoronoiPrism3dVertexMesh/SimpleMesh",
                                                  "20x12x1 4relax ApicalArea1.23");

        std::vector<double> cell_ids;
        for (unsigned id=0 ; id < p_mesh_a->GetNumElements() ; id++)
        {
            cell_ids.push_back(double(id));
        }

        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);
        vertex_mesh_writer.WriteVtkUsingMesh(*p_mesh_a);

    }

    void TestBoundaryNodes() throw(Exception)
    {
        // Generate a mesh that is 3 cells wide in x, 2 cells wide in y, 1 unit high in z, with 3 Lloyd's
        // relaxation steps and target average element area 100.0
        VoronoiPrism3dVertexMeshGenerator generator(3, 2, 5, 3, 100.0);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        VertexMeshWriter<3, 3> vertex_mesh_writer("TestVoronoiPrism3dVertexMesh/RowOfThree",
                                                  "3x2x5 3relax ApicalArea100");

        std::vector<double> cell_ids;
        for (unsigned id=0 ; id < p_mesh->GetNumElements() ; id++)
        {
            cell_ids.push_back(double(id));
        }

        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);
        vertex_mesh_writer.WriteVtkUsingMesh(*p_mesh);

        // Check basic mesh properties are correct
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 44u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 6u);

        // Check if all nodes are boundary nodes
        unsigned num_boundary_nodes = 0;
        for (unsigned node_idx = 0 ; node_idx < p_mesh->GetNumNodes() ; node_idx++ )
        {
            if (p_mesh->GetNode(node_idx)->IsBoundaryNode())
            {
                num_boundary_nodes++;
            }
        }
        TS_ASSERT_EQUALS(num_boundary_nodes, p_mesh->GetNumNodes());
    }

    void TestConstructorExceptions() throw(Exception)
    {
        // Throws because first parameter < 2
        TS_ASSERT_THROWS_THIS(VoronoiPrism3dVertexMeshGenerator generator(1, 9, 2, 1.23),
                "Need at least 2 by 2 cells");

        // Throws because second parameter < 2
        TS_ASSERT_THROWS_THIS(VoronoiPrism3dVertexMeshGenerator generator(9, 1, 2, 1.23),
                "Need at least 2 by 2 cells");

        //Throws because third parameter <= 0.0
        TS_ASSERT_THROWS_THIS(VoronoiPrism3dVertexMeshGenerator generator(9, 9, -2, 1.23),
                              "Specified element height must be strictly positive");

        // Throws because fifth parameter <= 0.0
        TS_ASSERT_THROWS_THIS(VoronoiPrism3dVertexMeshGenerator generator(9, 9, 2, 4, -1.23),
                              "Specified target apical area must be strictly positive");

    }

    void TestValidateSeedLocations() throw(Exception) //todo maybe check again the tests
    {
        // The instance of the generator class
        VoronoiPrism3dVertexMeshGenerator generator(3, 3, 1.0, 0, 1.0);

        // A helper vector
        c_vector<double, 2> one_one;
        one_one[0] = 1.0;
        one_one[1] = 1.0;

        // The sampling multiplier that should be used by the generator
        double sampling_multiplier = 0.5 * double(INT_MAX);

        // Vector of points that will be passed to ValidateSeedLocations()
        std::vector<c_vector<double, 2> > points;

        // Test two points at 0,0 are moved as expected
        points.push_back(zero_vector<double>(2));
        points.push_back(zero_vector<double>(2));

        generator.ValidateSeedLocations(points);

        TS_ASSERT_DELTA(points[0][0], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[0][1], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[1][0], 1.5 / sampling_multiplier, 1e-10);
        TS_ASSERT_DELTA(points[1][1], 0.0, 1e-10);

        // Test three points at 0,0 are moved as expected
        points.clear();
        points.push_back(zero_vector<double>(2));
        points.push_back(zero_vector<double>(2));
        points.push_back(zero_vector<double>(2));

        generator.ValidateSeedLocations(points);

        TS_ASSERT_DELTA(points[0][0], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[0][1], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[1][0], 1.5 / sampling_multiplier, 1e-10);
        TS_ASSERT_DELTA(points[1][1], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[2][0], 3.0 / sampling_multiplier, 1e-10);
        TS_ASSERT_DELTA(points[2][1], 0.0, 1e-10);

        // Test that periodicity is working correctly
        points.clear();
        points.push_back(one_one);
        points.push_back(one_one);
        generator.ValidateSeedLocations(points);

        TS_ASSERT_DELTA(points[0][0], 1.0, 1e-10);
        TS_ASSERT_DELTA(points[0][1], 1.0, 1e-10);
        TS_ASSERT_DELTA(points[1][0], 1.5 / sampling_multiplier, 1e-10);
        TS_ASSERT_DELTA(points[1][1], 1.0, 1e-10);

        // Test that two close non-equal points are moved correctly
        points.clear();
        points.push_back(zero_vector<double>(2));
        points.push_back( (2.0 * DBL_EPSILON) * one_one );
        generator.ValidateSeedLocations(points);

        TS_ASSERT_DELTA(points[0][0], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[0][1], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[1][0], 1.5 / (sampling_multiplier * sqrt(2.0)), 1e-10);
        TS_ASSERT_DELTA(points[1][1], 1.5 / (sampling_multiplier * sqrt(2.0)), 1e-10);

        // Test that periodicity for two close non-equal points is working as expected
        points.clear();
        points.push_back(one_one);
        points.push_back((1.0 + 2.0 * DBL_EPSILON) * one_one);
        generator.ValidateSeedLocations(points);

        TS_ASSERT_DELTA(points[0][0], 1.0, 1e-10);
        TS_ASSERT_DELTA(points[0][1], 1.0, 1e-10);
        TS_ASSERT_DELTA(points[1][0], 1.5 / (sampling_multiplier * sqrt(2.0)), 1e-10);
        TS_ASSERT_DELTA(points[1][1], 1.5 / (sampling_multiplier * sqrt(2.0)), 1e-10);
    }

    void TestGetToroidalMesh() throw(Exception)
    {
/*#if BOOST_VERSION >= 105200

        // Generate and get a toroidal mesh
        VoronoiPrism3dVertexMeshGenerator generator(19, 11, 7, 1.0);
        Toroidal2dVertexMesh* p_tor_mesh = generator.GetToroidalMesh();

        // Check nothing happens when we ReMesh()
        TS_ASSERT_THROWS_NOTHING(p_tor_mesh->ReMesh());

#endif // BOOST_VERSION >= 105200
*/    }

    void TestNodesAreRepositionedInToroidalMesh() throw(Exception)
    {
/*#if BOOST_VERSION >= 105200

        // Generate a toroidal mesh, all of whose nodes should lie within a prescribed bounding box
        unsigned num_x = 9;
        unsigned num_y = 11;
        unsigned num_relaxation_steps = 1;
        double area = 1.23;

        double width = num_x * sqrt(area);
        double height = num_y * sqrt(area);

        VoronoiPrism3dVertexMeshGenerator generator(num_x, num_y, num_relaxation_steps, area);
        Toroidal2dVertexMesh* p_tor_mesh = generator.GetToroidalMesh();

        // Verify that all node locations have been moved to be >= 0 and <= width & height
        for (unsigned node_idx = 0 ; node_idx < p_tor_mesh->GetNumNodes() ; node_idx++ )
        {
            c_vector<double, 2> this_location = p_tor_mesh->GetNode(node_idx)->rGetLocation();

            TS_ASSERT(this_location[0] >= 0.0);
            TS_ASSERT(this_location[1] >= 0.0);
            TS_ASSERT(this_location[0] < width);
            TS_ASSERT(this_location[1] < height);
        }

#endif // BOOST_VERSION >= 105200
*/    }

    void TestGetPolygonDistributionAndAreaVariation() throw(Exception) //todo expand it
    {
        unsigned num_x = 3;
        unsigned num_y = 4;
        double height_z =5;
        unsigned num_relaxation_steps = 1;
        double apical_area = 1.23;

        VoronoiPrism3dVertexMeshGenerator generator(num_x, num_y, height_z, num_relaxation_steps, apical_area);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        VertexMeshWriter<3, 3> vertex_mesh_writer("TestVoronoiPrism3dVertexMesh/PolygonDistributionAndStuff",
                                                  "3x4x1 1relax");

        std::vector<double> cell_ids;
        for (unsigned id=0 ; id < p_mesh->GetNumElements() ; id++)
        {
            cell_ids.push_back(double(id));
        }

        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);
        vertex_mesh_writer.WriteVtkUsingMesh(*p_mesh);

        // Get the polgyon distribution and check it
        std::vector<double> polygon_dist = generator.GetPolygonDistribution();

        TS_ASSERT(polygon_dist.size() > 0);

        //for manual visual validation
        MARK; PRINT_VECTOR(polygon_dist);
        double cumulative_proportion = 0.0;
        for (unsigned poly_idx = 0 ; poly_idx < polygon_dist.size() ; poly_idx++)
        {
            cumulative_proportion += polygon_dist[poly_idx];
        }

        TS_ASSERT_DELTA(cumulative_proportion, 1.0, 1e-6);

        // Get the area variation coefficient and check it
        double area_variation = generator.GetApicalAreaCoefficientOfVariation();

        TS_ASSERT_DELTA(area_variation, 0.32886, 1e-6);
        VoronoiPrism3dVertexMeshGenerator generator2(num_x, num_y, height_z + 5, num_relaxation_steps, apical_area);    //todo get same seed
        TS_ASSERT_DELTA(area_variation, generator2.GetApicalAreaCoefficientOfVariation(), 1e-6);
    }

    void TestRefreshSeedsAndRegenerateMeshCoverage() throw(Exception)
    {
        unsigned num_x = 3;
        unsigned num_y = 4;
        unsigned height_z = 3;
        unsigned num_relaxation_steps = 1;
        double area = 1.23;

        VoronoiPrism3dVertexMeshGenerator generator(num_x, num_y, height_z, num_relaxation_steps, area);

        generator.RefreshSeedsAndRegenerateMesh();
    }

    void TestSetAndGetMethods() throw(Exception)
    {
        unsigned num_x = 3;
        unsigned num_y = 4;
        double height_z = 5;
        unsigned num_relaxation_steps = 1;
        double area = 1.23;

        VoronoiPrism3dVertexMeshGenerator generator(num_x, num_y, height_z, num_relaxation_steps, area);

        generator.SetMaxExpectedNumSidesPerPolygon(5);
        TS_ASSERT_EQUALS(generator.GetMaxExpectedNumSidesPerPolygon(), 5u);
    }

    void TestDummyClassCoverage()
    {
#if BOOST_VERSION < 105200

        TS_ASSERT_THROWS_THIS(VoronoiPrism3dVertexMeshGenerator generator,
                "This is a dummy class. Build with Boost version 1.52 or above for functionality.");
        WARNING("Build with Boost version 1.52 or above for functionality.");
#endif // BOOST_VERSION < 105200
    }
};

#endif /*TESTVORONOIPRISM3DVERTEXMESHGENERATOR_HPP_*/
