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

#ifndef TESTHONEYCOMBMESHGENERATOR_HPP_
#define TESTHONEYCOMBMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "HoneycombMeshGenerator.hpp"

#include "FakePetscSetup.hpp"

class TestHoneycombMeshGenerator : public CxxTest::TestSuite
{
private:

    void Output2DNodesToFile(MutableMesh<2,2>* pMesh, std::string fileName)
    {
        OutputFileHandler handler("");
        out_stream file = handler.OpenOutputFile(fileName);

        unsigned num_nodes = pMesh->GetNumNodes();

        for (unsigned i=0; i<num_nodes; i++)
        {
            c_vector<double, 2> location;
            location = pMesh->GetNode(i)->rGetLocation();
            (*file) << location[0] << "\t" << location[1] << "\n" << std::flush;
        }

        file->close();
    }

public:

    void TestSimpleMesh()
    {
        unsigned cells_across = 2;
        unsigned cells_up = 2;
        double crypt_width = 0.5;
        unsigned thickness_of_ghost_layer = 1;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);

        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 16u);

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        TS_ASSERT_EQUALS(location_indices.size(), 4u);

        TS_ASSERT_DELTA(generator.GetDomainDepth(), 0.4330, 1e-4);
        TS_ASSERT_DELTA(generator.GetDomainWidth(), 0.5000, 1e-4);
    }

    void TestMonolayerHoneycombMeshGeneratorRelaxed()
    {
        int num_cells_width = 8;
        int num_cells_depth = 22;
        double width = 8.0;
        unsigned ghosts = 2;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, ghosts);
        double length = (double)num_cells_depth*(sqrt(3.0)/2)*width/(double)num_cells_width;

        // Check the mesh
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS((unsigned)p_mesh->GetNumNodes(),(num_cells_width+2*ghosts)*(num_cells_depth+2*ghosts));

        // Scaling factor
        double spooky = (double) ghosts;

        // Zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0], -spooky, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1], -spooky*sqrt(3.0)/2, 1e-6);

        unsigned this_many_ghosts_at_start = ((2*ghosts+num_cells_width)*ghosts+ghosts);

        // First real node
        TS_ASSERT_DELTA(p_mesh->GetNode(this_many_ghosts_at_start)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(this_many_ghosts_at_start)->GetPoint()[1], 0.0, 1e-6);

        // Last real node
        int index = (2*ghosts+num_cells_width)*(ghosts+num_cells_depth)+ghosts+num_cells_width;
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], width, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], length, 1e-4);

        // Last node
        int last_node = p_mesh->GetNumNodes()-1;
        double last_node_y = length+(spooky-1)*(sqrt(3.0)/2);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[0], width+(spooky-0.5), 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[1], last_node_y, 1e-5);

        // Check the ghost nodes
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        TS_ASSERT_EQUALS(location_indices.size(), p_mesh->GetNumNodes() - 2*(ghosts*(num_cells_width + 2*ghosts + num_cells_depth)));

        std::set<unsigned> correct_ghost_node_indices;
        for (unsigned i=0; i<this_many_ghosts_at_start; i++)
        {
            correct_ghost_node_indices.insert(i);
        }
        correct_ghost_node_indices.insert( this_many_ghosts_at_start+num_cells_width+1 );

        // Create a set of node indices corresponding to ghost nodes
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices_set;
        std::set<unsigned> ghost_node_indices;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            node_indices.insert(p_mesh->GetNode(i)->GetIndex());
        }
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            location_indices_set.insert(location_indices[i]);
        }

        std::set_difference(node_indices.begin(), node_indices.end(),
                            location_indices_set.begin(), location_indices_set.end(),
                            std::inserter(ghost_node_indices, ghost_node_indices.begin()));

        bool all_included = includes(ghost_node_indices.begin(), ghost_node_indices.end(),
                                     correct_ghost_node_indices.begin(),correct_ghost_node_indices.end());

        TS_ASSERT_EQUALS(all_included, true);

        TS_ASSERT_DELTA(p_mesh->GetWidth(0u), (width-1.0) + 2.0 * (double) ghosts + 0.5, 1e-7); // 11.5 as mesh is stagered and includes ghosts
        TS_ASSERT_DELTA(p_mesh->GetWidth(1u), 21.6506, 1e-4); // includes ghosts
    }

    void TestMonolayerHoneycombMeshGeneratorCompressed()
    {
        int num_cells_width = 8;
        int num_cells_depth = 12;
        double width = 6.0;
        unsigned ghosts = 4;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, ghosts, width/num_cells_width);
        double length = (double)num_cells_depth*(sqrt(3.0)/2)*width/(double)num_cells_width;

        // Check the mesh
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS((unsigned)p_mesh->GetNumNodes(), (num_cells_width+2*ghosts)*(num_cells_depth+2*ghosts));

        // Scaling factor
        double factor = (width/(double)num_cells_width);
        double spooky = (double) ghosts;

        // Zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0], -spooky*factor, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1], -spooky*factor*sqrt(3.0)/2, 1e-6);

        unsigned this_many_ghosts_at_start = ((2*ghosts+num_cells_width)*ghosts+ghosts);

        // First real node
        TS_ASSERT_DELTA(p_mesh->GetNode(this_many_ghosts_at_start)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(this_many_ghosts_at_start)->GetPoint()[1], 0.0, 1e-6);

        // Last real node
        int index = (2*ghosts+num_cells_width)*(ghosts+num_cells_depth)+ghosts+num_cells_width;
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], width, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], length, 1e-4);

        // Last node
        int last_node = p_mesh->GetNumNodes()-1;
        double last_node_y = length+(spooky-1)*factor*(sqrt(3.0)/2);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[0], width+(spooky-0.5)*factor, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[1], last_node_y, 1e-6);

        // Check the ghost nodes
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        TS_ASSERT_EQUALS(location_indices.size(), p_mesh->GetNumNodes() - 2*(ghosts*(num_cells_width + 2*ghosts + num_cells_depth)));

        std::set<unsigned> correct_ghost_node_indices;
        for (unsigned i=0; i<this_many_ghosts_at_start; i++)
        {
            correct_ghost_node_indices.insert(i);
        }
        correct_ghost_node_indices.insert( this_many_ghosts_at_start+num_cells_width+1 );

        // Create a set of node indices corresponding to ghost nodes
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices_set;
        std::set<unsigned> ghost_node_indices;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            node_indices.insert(p_mesh->GetNode(i)->GetIndex());
        }
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            location_indices_set.insert(location_indices[i]);
        }

        std::set_difference(node_indices.begin(), node_indices.end(),
                            location_indices_set.begin(), location_indices_set.end(),
                            std::inserter(ghost_node_indices, ghost_node_indices.begin()));

        bool all_included = includes(ghost_node_indices.begin(), ghost_node_indices.end(),
                                     correct_ghost_node_indices.begin(),correct_ghost_node_indices.end());

        TS_ASSERT_EQUALS(all_included, true);

        TS_ASSERT_DELTA(p_mesh->GetWidth(0u),  11.6250 , 1e-4); // as mesh is staggered, compressed and includes ghosts
        TS_ASSERT_DELTA(p_mesh->GetWidth(1u), 12.3408, 1e-4); // includes ghosts
    }

    void TestBoundaryNodes()
    {
        unsigned cells_across = 4;
        unsigned cells_up = 4;
        double crypt_width = 4;
        unsigned thickness_of_ghost_layer = 0;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 16u);

        unsigned num_non_boundary_nodes = 0;
        for (unsigned node_index=0; node_index<16u; node_index++)
        {
            if (!p_mesh->GetNode(node_index)->IsBoundaryNode())
            {
                num_non_boundary_nodes++;
            }
        }
        TS_ASSERT_EQUALS(num_non_boundary_nodes, 4u);
    }

    void TestGetCircularMesh()
    {
        unsigned num_cells_depth = 10;
        unsigned num_cells_width = 10;
        double radius = 3.5;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(radius);

        double epsilon = 1e-5;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TS_ASSERT_LESS_THAN_EQUALS(norm_2(p_mesh->GetNode(i)->rGetLocation()), radius+epsilon);
        }
        Output2DNodesToFile(p_mesh, "circular_mesh.dat");
    }

    void TestGetCircularMeshWithGhostNodesThrowsException()
    {
        unsigned num_cells_depth = 10;
        unsigned num_cells_width = 10;
        double radius = 3.5;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 3);

        TS_ASSERT_THROWS_THIS(generator.GetCircularMesh(radius), "Cannot call GetCircularMesh on a HoneycombMesh with ghost nodes");
    }

    void TestCircularMeshIsJacobian()
    {
        HoneycombMeshGenerator generator(20, 20, 0);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(10);

        NodeMap map(p_mesh->GetNumAllNodes());
        p_mesh->ReMesh(map);
    }

    void TestLargeMesh()
    {
        unsigned cells_across = 100;
        unsigned cells_up = 100;
        unsigned thickness_of_ghost_layer = 1;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 10404u);

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        TS_ASSERT_EQUALS(location_indices.size(), 10000u);

        TS_ASSERT_DELTA(generator.GetDomainDepth(), 86.6025, 1e-4);
        TS_ASSERT_DELTA(generator.GetDomainWidth(), 100, 1e-4);
    }
};

#endif /*TESTHONEYCOMBMESHGENERATOR_HPP_*/
