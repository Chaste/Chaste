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

#ifndef TESTHONEYCOMBMESHGENERATOR_HPP_
#define TESTHONEYCOMBMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "HoneycombMeshGenerator.hpp"

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
            c_vector<double, 2> location = pMesh->GetNode(i)->rGetLocation();
            (*file) << location[0] << "\t" << location[1] << "\n" << std::flush;
        }

        file->close();
    }

public:

    void TestSimpleMesh() throw(Exception)
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

    void TestMonolayerHoneycombMeshGeneratorRelaxed() throw(Exception)
    {
        int num_cells_width = 8;
        int num_cells_depth = 22;
        double width = 8.0;
        unsigned ghosts = 2;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, ghosts);
        double length = (double)num_cells_depth*(sqrt(3)/2)*width/(double)num_cells_width;

        // Check the mesh
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS((unsigned)p_mesh->GetNumNodes(),(num_cells_width+2*ghosts)*(num_cells_depth+2*ghosts));

        // Scaling factor
        double spooky = (double) ghosts;

        // Zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0], -spooky, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1], -spooky*sqrt(3)/2, 1e-6);

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
        double last_node_y = length+(spooky-1)*(sqrt(3)/2);
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

    void TestMonolayerHoneycombMeshGeneratorCompressed() throw(Exception)
    {
        int num_cells_width = 8;
        int num_cells_depth = 12;
        double width = 6.0;
        unsigned ghosts = 4;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, ghosts, width/num_cells_width);
        double length = (double)num_cells_depth*(sqrt(3)/2)*width/(double)num_cells_width;

        // Check the mesh
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS((unsigned)p_mesh->GetNumNodes(), (num_cells_width+2*ghosts)*(num_cells_depth+2*ghosts));

        // Scaling factor
        double factor = (width/(double)num_cells_width);
        double spooky = (double) ghosts;

        // Zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0], -spooky*factor, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1], -spooky*factor*sqrt(3)/2, 1e-6);

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
        double last_node_y = length+(spooky-1)*factor*(sqrt(3)/2);
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

    void TestBoundaryNodes() throw(Exception)
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

    void TestGetCircularMesh() throw(Exception)
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

    void TestGetCircularMeshWithGhostNodesThrowsException() throw(Exception)
    {
        unsigned num_cells_depth = 10;
        unsigned num_cells_width = 10;
        double radius = 3.5;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 3);

        TS_ASSERT_THROWS_THIS(generator.GetCircularMesh(radius), "Cannot call GetCircularMesh on a HoneycombMesh with ghost nodes");
    }

    void TestCircularMeshIsJacobian() throw(Exception)
    {
        HoneycombMeshGenerator generator(20, 20, 0);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(10);

        NodeMap map(p_mesh->GetNumAllNodes());
        p_mesh->ReMesh(map);
    }

    void TestLargeMesh() throw(Exception)
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
