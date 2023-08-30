/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TESTTOROIDALHONEYCOMBMESHGENERATOR_HPP_
#define TESTTOROIDALHONEYCOMBMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>

#include "ToroidalHoneycombMeshGenerator.hpp"

#include "FakePetscSetup.hpp"

class TestToroidalHoneycombMeshGenerator : public CxxTest::TestSuite
{
private:

    void Output2DNodesToFileToroidal(boost::shared_ptr<Toroidal2dMesh> pMesh, std::string fileName)
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

    void TestToroidalHoneycombMeshGeneratorRelaxed()
    {
        unsigned num_cells_width = 8;
        unsigned num_cells_height = 22;
        double width = 8.0;
        double height = 22.0*sqrt(3)/2.0;

        double x_factor = width/(double)num_cells_width;
        double y_factor = height/(double)num_cells_height/sqrt(3.0)*2.0;

        ToroidalHoneycombMeshGenerator generator(num_cells_width, num_cells_height, x_factor, y_factor);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        // Check the mesh
        TS_ASSERT_THROWS_THIS(generator.GetMesh(),"A Toroidal mesh was created but a normal mesh is being requested.");

        Output2DNodesToFileToroidal(p_mesh, "Toroidal_node_positions.dat");
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), (num_cells_width)*(num_cells_height));

        // Zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1], 0.0, 1e-5);

        // Last node
        int last_node = p_mesh->GetNumNodes()-1;
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[0], 7.5, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[1], (num_cells_height-1)*sqrt(3.0)/2.0, 1e-4);

        TS_ASSERT_DELTA(p_mesh->GetWidth(0u), width, 1e-7);
        TS_ASSERT_DELTA(p_mesh->GetWidth(1u), height, 1e-4);
    }

    void TestToroidalHoneycombMeshGeneratorCompressed()
    {
        unsigned num_cells_width = 8;
        unsigned num_cells_height = 8;
        double width = 7.0;
        double height = 4.0;

        double x_factor = width/(double)num_cells_width;
        double y_factor = height/(double)num_cells_height/sqrt(3.0)*2.0;

        ToroidalHoneycombMeshGenerator generator(num_cells_width, num_cells_height, x_factor, y_factor);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        // Check the mesh
        TS_ASSERT_THROWS_THIS(generator.GetMesh(),"A Toroidal mesh was created but a normal mesh is being requested.");

        Output2DNodesToFileToroidal(p_mesh, "Toroidal_node_positions.dat");
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),(num_cells_width)*(num_cells_height));

        // First node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1], 0.0, 1e-12);

        // Last node
        int last_node = p_mesh->GetNumNodes()-1;
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[0], x_factor*7.5, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[1], y_factor*(num_cells_height-1)*sqrt(3.0)/2.0, 1e-4);

        TS_ASSERT_DELTA(p_mesh->GetWidth(0u), width, 1e-7);
        TS_ASSERT_DELTA(p_mesh->GetWidth(1u), height, 1e-4);
    }
};

#endif /*TESTTOROIDALHONEYCOMBMESHGENERATOR_HPP_*/
