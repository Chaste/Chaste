/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTCYLINDRICALHONEYCOMBMESHGENERATOR_HPP_
#define TESTCYLINDRICALHONEYCOMBMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "CylindricalHoneycombMeshGenerator.hpp"

class TestCylindricalHoneycombMeshGenerator : public CxxTest::TestSuite
{
private:

    void Output2DNodesToFileCylindrical(Cylindrical2dMesh* pMesh, std::string fileName)
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

    void TestCylindricalHoneycombMeshGeneratorRelaxed() throw(Exception)
    {
        unsigned num_cells_width = 8;
        unsigned num_cells_depth = 22;
        unsigned ghosts = 2;

        CylindricalHoneycombMeshGenerator generator(num_cells_width, num_cells_depth, ghosts);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Check the mesh
        TS_ASSERT_THROWS_THIS(generator.GetMesh(),"A cylindrical mesh was created but a normal mesh is being requested.");

        Output2DNodesToFileCylindrical(p_mesh, "cylindrical_node_positions.dat");
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), (num_cells_width)*(num_cells_depth+2*ghosts));

        // Zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1], -(double)ghosts*sqrt(3.0)/2.0, 1e-5);

        // First real node
        int index = (num_cells_width)*ghosts; // 4 here is the number of ghost nodes in a row
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], 0.0, 1e-12);

        // Last real node
        index = (ghosts+num_cells_depth)*(num_cells_width)-1;
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], 7.5, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], 21.0*sqrt(3)/2.0, 1e-4);

        // Last node
        int last_node = p_mesh->GetNumNodes()-1;
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[0], 7.5, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[1], (ghosts+num_cells_depth-1)*sqrt(3.0)/2.0, 1e-4);

        // Check the ghost nodes
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        TS_ASSERT_EQUALS(location_indices.size(), p_mesh->GetNumNodes() - 2*(ghosts*(num_cells_width)));

        std::set<unsigned> correct_ghost_node_indices;

        for (unsigned i=0; i< num_cells_width*ghosts; i++)
        {
            correct_ghost_node_indices.insert(i);
        }
        correct_ghost_node_indices.insert( (ghosts+num_cells_depth)*num_cells_width+1 );

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
                                     correct_ghost_node_indices.begin(), correct_ghost_node_indices.end());

        TS_ASSERT_EQUALS(all_included, true);

        TS_ASSERT_DELTA(p_mesh->GetWidth(0u), (double)num_cells_width, 1e-7);
        TS_ASSERT_DELTA(p_mesh->GetWidth(1u), 21.6506, 1e-4);
    }

    void TestCylindricalHoneycombMeshGeneratorCompressed() throw(Exception)
    {
        unsigned num_cells_width = 8;
        unsigned num_cells_depth = 22;
        double width = 7.0;
        unsigned ghosts = 2;

        double x_factor = width/(double)num_cells_width;

        CylindricalHoneycombMeshGenerator generator(num_cells_width, num_cells_depth, ghosts, width/num_cells_width);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Check the mesh
        TS_ASSERT_THROWS_THIS(generator.GetMesh(),"A cylindrical mesh was created but a normal mesh is being requested.");

        Output2DNodesToFileCylindrical(p_mesh, "cylindrical_node_positions.dat");
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),(num_cells_width)*(num_cells_depth+2*ghosts));

        // Zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1], -x_factor*(double)ghosts*sqrt(3.0)/2.0, 1e-5);

        // First real node
        int index = (num_cells_width)*ghosts; // 4 here is the number of ghost nodes in a row
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], 0.0, 1e-12);

        // Last real node
        index = (ghosts+num_cells_depth)*(num_cells_width)-1;
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], x_factor*7.5, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], x_factor*21.0*sqrt(3)/2.0, 1e-4);

        // Last node
        int last_node = p_mesh->GetNumNodes()-1;
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[0], x_factor*7.5, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[1], x_factor*(ghosts+num_cells_depth-1)*sqrt(3.0)/2.0, 1e-4);

        // Check the ghost nodes
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        TS_ASSERT_EQUALS(location_indices.size(), p_mesh->GetNumNodes() - 2*(ghosts*(num_cells_width)));

        std::set<unsigned> correct_ghost_node_indices;

        for (unsigned i=0; i< num_cells_width*ghosts; i++)
        {
            correct_ghost_node_indices.insert(i);
        }
        correct_ghost_node_indices.insert( (ghosts+num_cells_depth)*num_cells_width+1 );

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

        TS_ASSERT_DELTA(p_mesh->GetWidth(0u), x_factor*(double)num_cells_width, 1e-7);
        TS_ASSERT_DELTA(p_mesh->GetWidth(1u), 18.9443, 1e-4);
    }
};

#endif /*TESTCYLINDRICALHONEYCOMBMESHGENERATOR_HPP_*/
