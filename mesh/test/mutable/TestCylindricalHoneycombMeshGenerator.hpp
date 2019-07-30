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

#ifndef TESTCYLINDRICALHONEYCOMBMESHGENERATOR_HPP_
#define TESTCYLINDRICALHONEYCOMBMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "CylindricalHoneycombMeshGenerator.hpp"

#include "FakePetscSetup.hpp"

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
            c_vector<double, 2> location;
            location = pMesh->GetNode(i)->rGetLocation();
            (*file) << location[0] << "\t" << location[1] << "\n" << std::flush;
        }

        file->close();
    }

public:

    void TestCylindricalHoneycombMeshGeneratorRelaxed()
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
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], 21.0*sqrt(3.0)/2.0, 1e-4);

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

    void TestCylindricalHoneycombMeshGeneratorCompressed()
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
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], x_factor*21.0*sqrt(3.0)/2.0, 1e-4);

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
