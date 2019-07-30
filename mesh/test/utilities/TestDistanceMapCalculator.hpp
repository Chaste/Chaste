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
#ifndef TESTDISTANCEMAPCALCULATOR_
#define TESTDISTANCEMAPCALCULATOR_

#include "TrianglesMeshReader.hpp"
#include "DistanceMapCalculator.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "RandomNumberGenerator.hpp"

class TestDistanceMapCalculator : public CxxTest::TestSuite
{
public:
    void TestDistances1D()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");

        std::vector<unsigned> map_origin;
        map_origin.push_back(0u);
        std::vector<double> distances_serial;
        {
            //This is in a block so that we can minimise to scope of the serial mesh (to avoid using it in error)
            TetrahedralMesh<1,1> serial_mesh;
            serial_mesh.ConstructFromMeshReader(mesh_reader);
            TS_ASSERT_EQUALS(serial_mesh.GetNumNodes(), 11u);
            TS_ASSERT_EQUALS(serial_mesh.GetNumElements(), 10u);
            TS_ASSERT_EQUALS(serial_mesh.GetNumBoundaryElements(), 2u);
            DistanceMapCalculator<1,1> distance_calculator_serial(serial_mesh);
            distance_calculator_serial.ComputeDistanceMap(map_origin, distances_serial);
        }

        DistributedTetrahedralMesh<1,1> parallel_mesh;
        parallel_mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(parallel_mesh.GetNumNodes(), 11u);
        TS_ASSERT_EQUALS(parallel_mesh.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(parallel_mesh.GetNumBoundaryElements(), 2u);

        DistanceMapCalculator<1,1> distance_calculator_parallel(parallel_mesh);
        std::vector<double> distances_parallel;
        distance_calculator_parallel.ComputeDistanceMap(map_origin, distances_parallel);

        TS_ASSERT_EQUALS(distances_serial.size(), distances_parallel.size());
        for (unsigned index=0; index<distances_parallel.size(); index++)
        {
            try
            {
                // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
                c_vector<double, 1> node;
                node = parallel_mesh.GetNode(index)->rGetLocation(); // throws if not owned

                TS_ASSERT_DELTA(distances_serial[index], node(0), 1e-12);
                TS_ASSERT_DELTA(distances_parallel[index], node(0), 1e-12);
            }
            catch (Exception&)
            {
            }
        }
    }

    void TestDistancesToCorner()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_21_nodes_side/Cube21"); // 5x5x5mm cube (internode distance = 0.25mm)

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        unsigned num_nodes=9261u;
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes); // 21x21x21 nodes
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 48000u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4800u);

        DistributedTetrahedralMesh<3,3> parallel_mesh(DistributedTetrahedralMeshPartitionType::DUMB); // No reordering;
        parallel_mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(parallel_mesh.GetNumNodes(), num_nodes); // 21x21x21 nodes
        TS_ASSERT_EQUALS(parallel_mesh.GetNumElements(), 48000u);
        TS_ASSERT_EQUALS(parallel_mesh.GetNumBoundaryElements(), 4800u);

        unsigned far_index=9260u;
        c_vector<double,3> far_corner=mesh.GetNode(far_index)->rGetLocation();
        TS_ASSERT_DELTA( far_corner[0], 0.25, 1e-11);
        TS_ASSERT_DELTA( far_corner[1], 0.25, 1e-11);
        TS_ASSERT_DELTA( far_corner[2], 0.25, 1e-11);
        try
        {
            c_vector<double,3> parallel_far_corner=parallel_mesh.GetNode(far_index)->rGetLocation();
            TS_ASSERT_DELTA( parallel_far_corner[0], 0.25, 1e-11);
            TS_ASSERT_DELTA( parallel_far_corner[1], 0.25, 1e-11);
            TS_ASSERT_DELTA( parallel_far_corner[2], 0.25, 1e-11);
        }
        catch (Exception&)
        {
        }

        std::vector<unsigned> map_far_corner;
        map_far_corner.push_back(far_index);

        DistanceMapCalculator<3,3> distance_calculator(mesh);
        std::vector<double> distances;
        distance_calculator.ComputeDistanceMap(map_far_corner, distances);

        DistanceMapCalculator<3,3> parallel_distance_calculator(parallel_mesh);
        std::vector<double> parallel_distances;
        parallel_distance_calculator.ComputeDistanceMap(map_far_corner, parallel_distances);

        TS_ASSERT_EQUALS(distance_calculator.mRoundCounter, 1u);
        //Nodes in mesh are order such that a dumb partitioning will give a sequential handover from proc0 to proc1...
        TS_ASSERT_EQUALS(parallel_distance_calculator.mRoundCounter, PetscTools::GetNumProcs());
        //Note unsigned division is okay here
        TS_ASSERT_DELTA(parallel_distance_calculator.mPopCounter, num_nodes/PetscTools::GetNumProcs(), 1u);
        TS_ASSERT_DELTA(distance_calculator.mPopCounter, num_nodes, 1u);
        for (unsigned index=0; index<distances.size(); index++)
        {
            c_vector<double, 3> node = mesh.GetNode(index)->rGetLocation();

            //Straightline distance
            double euclidean_distance = norm_2(far_corner - node);
            // x + y + z distance
            double manhattan_distance = norm_1(far_corner - node);
            //If they differ, then allow the in-mesh distance to be in between
            double error_bound = (manhattan_distance - euclidean_distance)/2.0;
            //If they don't differ, then we expect the in-mesh distance to be similar
            if (error_bound < 1e-15)
            {
                error_bound = 1e-15;
            }
            TS_ASSERT_LESS_THAN_EQUALS(distances[index], manhattan_distance+DBL_EPSILON);
            TS_ASSERT_LESS_THAN_EQUALS(euclidean_distance, distances[index]+DBL_EPSILON);
            TS_ASSERT_DELTA(distances[index], euclidean_distance, error_bound);

            TS_ASSERT_DELTA(distances[index], parallel_distances[index], 1e-15);
        }

        // Test some point-to-point distances
        RandomNumberGenerator::Instance()->Reseed(1);
        unsigned trials=25;
        unsigned pops=0;
        unsigned sequential_pops=0;
        for (unsigned i=0; i<trials; i++)
        {
            unsigned index=RandomNumberGenerator::Instance()->randMod(parallel_distances.size());
            TS_ASSERT_DELTA(parallel_distance_calculator.SingleDistance(9260u, index), parallel_distances[index], 1e-15);
            TS_ASSERT_DELTA(distance_calculator.SingleDistance(9260u, index), parallel_distances[index], 1e-15);
            pops += parallel_distance_calculator.mPopCounter;
            sequential_pops += distance_calculator.mPopCounter;
            TS_ASSERT_LESS_THAN_EQUALS(parallel_distance_calculator.mRoundCounter, PetscTools::GetNumProcs()+2);
        }

        // Without A*: TS_ASSERT_DELTA(sequential_pops/(double)trials, num_nodes/2, 300);
        TS_ASSERT_LESS_THAN(sequential_pops/(double)trials, num_nodes/20.0);
        if (PetscTools::IsSequential())
        {
            //Early termination
            TS_ASSERT_EQUALS(pops, sequential_pops);
        }
        else
        {
            //Early termination on remote processes is not yet possible
            //This may lead to multiple updates from remote
            //A* Leads to even more updates on average
            // Without A*: TS_ASSERT_DELTA(pops/(double)trials, num_nodes/PetscTools::GetNumProcs(), 700.0);
            TS_ASSERT_LESS_THAN(pops/(double)trials, num_nodes/10.0 );
         }

        //Reverse - to check that cached information is flushed.
        for (unsigned i=0; i<3; i++)
        {
            unsigned index=RandomNumberGenerator::Instance()->randMod(parallel_distances.size());
            TS_ASSERT_DELTA(parallel_distance_calculator.SingleDistance(index, 9260u), parallel_distances[index], 1e-15);
        }
    }

    void TestDistancesToFaceDumb()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_21_nodes_side/Cube21"); // 5x5x5mm cube (internode distance = 0.25mm)

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9261u); // 21x21x21 nodes
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 48000u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4800u);

        DistributedTetrahedralMesh<3,3> parallel_mesh(DistributedTetrahedralMeshPartitionType::DUMB); // No reordering
        parallel_mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(parallel_mesh.GetNumNodes(), 9261u); // 21x21x21 nodes
        TS_ASSERT_EQUALS(parallel_mesh.GetNumElements(), 48000u);
        TS_ASSERT_EQUALS(parallel_mesh.GetNumBoundaryElements(), 4800u);


        std::vector<unsigned> map_left;
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            // Get the nodes at the left face of the cube
            if (mesh.GetNode(index)->rGetLocation()[0] + 0.25 < 1e-6)
            {
                map_left.push_back(index);
            }
        }

        TS_ASSERT_EQUALS(map_left.size(), 21u*21u);

        DistanceMapCalculator<3,3> distance_calculator(mesh);
        std::vector<double> distances;
        distance_calculator.ComputeDistanceMap(map_left, distances);

        DistanceMapCalculator<3,3> parallel_distance_calculator(parallel_mesh);
        std::vector<double> parallel_distances;
        parallel_distance_calculator.ComputeDistanceMap(map_left, parallel_distances);

        TS_ASSERT_EQUALS(distance_calculator.mRoundCounter, 1u);
        TS_ASSERT_DELTA(parallel_distance_calculator.mRoundCounter, 2u, 1u);// 1 2 or 3

        for (unsigned index=0; index<distances.size(); index++)
        {
            // The distance should be equal to the x-coordinate of the point (minus the offset of the left face of the cube)
            c_vector<double, 3> node = mesh.GetNode(index)->rGetLocation();
            TS_ASSERT_DELTA(distances[index], node[0]+0.25,1e-11);
            TS_ASSERT_DELTA(parallel_distances[index], node[0]+0.25,1e-11);
        }
    }
    void TestDistancesToFace()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_21_nodes_side/Cube21"); // 5x5x5mm cube (internode distance = 0.25mm)

        DistributedTetrahedralMesh<3,3> parallel_mesh;
        parallel_mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(parallel_mesh.GetNumNodes(), 9261u); // 21x21x21 nodes
        TS_ASSERT_EQUALS(parallel_mesh.GetNumElements(), 48000u);
        TS_ASSERT_EQUALS(parallel_mesh.GetNumBoundaryElements(), 4800u);


        std::vector<unsigned> map_left;
        for (unsigned index=0; index<parallel_mesh.GetNumNodes(); index++)
        {
            // Get the *only local* nodes at the left face of the cube
            try
            {
                if (parallel_mesh.GetNode(index)->rGetLocation()[0] + 0.25 < 1e-6)
                {
                    map_left.push_back(index);
                }
            }
            catch (Exception&)
            {
            }
        }

        DistanceMapCalculator<3,3> parallel_distance_calculator(parallel_mesh);
        std::vector<double> parallel_distances;
        parallel_distance_calculator.ComputeDistanceMap(map_left, parallel_distances);

        for (unsigned index=0; index<parallel_distances.size(); index++)
        {
            try
            {
                // The distance should be equal to the x-coordinate of the point (minus the offset of the left face of the cube)
                c_vector<double, 3> node = parallel_mesh.GetNode(index)->rGetLocation();
                TS_ASSERT_DELTA(parallel_distances[index], node[0]+0.25, 1e-15);
            }
            catch (Exception&)
            {
                //If we don't know the geometry of this node, then we still know the distance, which ought to be in [0, 0.5 ] for left and  right faces at extremes
                TS_ASSERT_DELTA(parallel_distances[index], 0.25, 0.2500001);
            }
        }
    }

    void TestDistancesWithEmptySource()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_21_nodes_side/Cube21"); // 5x5x5mm cube (internode distance = 0.25mm)

        DistributedTetrahedralMesh<3,3> parallel_mesh;
        parallel_mesh.ConstructFromMeshReader(mesh_reader);


        std::vector<unsigned> empty_sources;
        TS_ASSERT_EQUALS(empty_sources.size(), 0U);

        DistanceMapCalculator<3,3> parallel_distance_calculator(parallel_mesh);
        std::vector<double> parallel_distances;
        parallel_distance_calculator.ComputeDistanceMap(empty_sources, parallel_distances);

        TS_ASSERT_EQUALS(parallel_distances.size(), 9261U);
        for (unsigned index=0; index<parallel_distances.size(); index++)
        {
            TS_ASSERT_EQUALS(parallel_distances[index], DBL_MAX);
        }
    }
};

#endif /*TESTDISTANCEMAPCALCULATOR_*/
