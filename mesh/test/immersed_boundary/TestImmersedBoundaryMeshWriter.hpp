/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYMESHWRITER_HPP_
#define TESTIMMERSEDBOUNDARYMESHWRITER_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

// Includes from projects/ImmersedBoundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryMeshWriter.hpp"
#include "ImmersedBoundaryMeshReader.hpp"
#include "FileComparison.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

/*
  This file is just missing a test which handles multiple overlaps
*/

class TestImmersedBoundaryMeshWriter : public CxxTest::TestSuite
{
public:

    void TestImmersedBoundaryMeshWriterIn2d()
    {
        /*
         * In this test, we generate a mesh with multiple elements and lamina, and check that it is written correctly.
         */
        std::vector<Node<2>*> nodes;

        std::vector<Node<2>*> nodes_elem1;
        std::vector<Node<2>*> nodes_elem2;
        std::vector<Node<2>*> nodes_elem3;

        std::vector<Node<2>*> nodes_lam1;
        std::vector<Node<2>*> nodes_lam2;

        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.05, 0.05 * sqrt(3)));

        nodes.push_back(new Node<2>(3, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(4, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(5, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(6, true, 0.0, 0.1));

        nodes.push_back(new Node<2>(7, true, 0.5 + 0.1 * cos(0.0 * M_PI / 8.0), 0.5 + 0.1 * sin(0.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(8, true, 0.5 + 0.1 * cos(2.0 * M_PI / 8.0), 0.5 + 0.1 * sin(2.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(9, true, 0.5 + 0.1 * cos(4.0 * M_PI / 8.0), 0.5 + 0.1 * sin(4.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(10, true, 0.5 + 0.1 * cos(6.0 * M_PI / 8.0), 0.5 + 0.1 * sin(6.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(11, true, 0.5 + 0.1 * cos(8.0 * M_PI / 8.0), 0.5 + 0.1 * sin(8.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(12, true, 0.5 + 0.1 * cos(10.0 * M_PI / 8.0), 0.5 + 0.1 * sin(10.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(13, true, 0.5 + 0.1 * cos(12.0 * M_PI / 8.0), 0.5 + 0.1 * sin(12.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(14, true, 0.5 + 0.1 * cos(14.0 * M_PI / 8.0), 0.5 + 0.1 * sin(14.0 * M_PI / 8.0)));

        nodes.push_back(new Node<2>(15, true, 0.1, 0.3));
        nodes.push_back(new Node<2>(16, true, 0.3, 0.3));
        nodes.push_back(new Node<2>(17, true, 0.5, 0.3));
        nodes.push_back(new Node<2>(18, true, 0.7, 0.3));
        nodes.push_back(new Node<2>(19, true, 0.9, 0.3));

        nodes.push_back(new Node<2>(20, true, 0.2, 0.2));
        nodes.push_back(new Node<2>(21, true, 0.2, 0.7));

        // Triangle
        nodes_elem1.push_back(nodes[0]);
        nodes_elem1.push_back(nodes[1]);
        nodes_elem1.push_back(nodes[2]);

        // Square
        nodes_elem2.push_back(nodes[3]);
        nodes_elem2.push_back(nodes[4]);
        nodes_elem2.push_back(nodes[5]);
        nodes_elem2.push_back(nodes[6]);

        // Octagon
        nodes_elem3.push_back(nodes[7]);
        nodes_elem3.push_back(nodes[8]);
        nodes_elem3.push_back(nodes[9]);
        nodes_elem3.push_back(nodes[10]);
        nodes_elem3.push_back(nodes[11]);
        nodes_elem3.push_back(nodes[12]);
        nodes_elem3.push_back(nodes[13]);
        nodes_elem3.push_back(nodes[14]);

        // Lam 1 (x varies)
        nodes_lam1.push_back(nodes[15]);
        nodes_lam1.push_back(nodes[16]);
        nodes_lam1.push_back(nodes[17]);
        nodes_lam1.push_back(nodes[18]);
        nodes_lam1.push_back(nodes[19]);

        // Lam 2 (y varies)
        nodes_lam2.push_back(nodes[20]);
        nodes_lam2.push_back(nodes[21]);

        std::vector<ImmersedBoundaryElement<2, 2>*> elems;
        elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes_elem1));
        elems[0]->AddCornerNode(nodes_elem1[0]);
        elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes_elem2));
        elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes_elem3));

        std::vector<ImmersedBoundaryElement<1, 2>*> lams;
        lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes_lam1));
        lams.push_back(new ImmersedBoundaryElement<1, 2>(1, nodes_lam2));

        ImmersedBoundaryMesh<2, 2> mesh(nodes, elems, lams, 8, 8);

        // Create a vertex mesh writer
        ImmersedBoundaryMeshWriter<2,2> ib_mesh_writer("TestIbMeshWriterIn2d", "ib_mesh_2d");

        // Test files are written correctly
        ib_mesh_writer.WriteFilesUsingMesh(mesh);

        OutputFileHandler handler("TestIbMeshWriterIn2d", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "ib_mesh_2d.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "ib_mesh_2d.elem";
        std::string results_file3 = handler.GetOutputDirectoryFullPath() + "ib_mesh_2d.lam";

        // Test adding point data
        ib_mesh_writer.AddPointData("test", {1.0, 1.0});

        FileComparison comparer1(results_file1,"mesh/test/data/TestIbMeshWriter/ib_mesh_2d.node");
        TS_ASSERT(comparer1.CompareFiles());

        FileComparison comparer2(results_file2,"mesh/test/data/TestIbMeshWriter/ib_mesh_2d.elem");
        TS_ASSERT(comparer2.CompareFiles());

        FileComparison comparer3(results_file3,"mesh/test/data/TestIbMeshWriter/ib_mesh_2d.lam");
        TS_ASSERT(comparer3.CompareFiles());
    }

    void TestImmersedBoundaryMeshWriterVTKCornerOverlap()
    {
#ifdef CHASTE_VTK
        // Overlap top right
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.1));
            nodes.push_back(new Node<2>(2, true, 0.9, 0.9));
            nodes.push_back(new Node<2>(3, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems);

            ImmersedBoundaryMeshWriter<2, 2> ib_mesh_writer("TestIbMeshWriterVTKCornerOverlap", "ib_mesh_vtk_corner_overlap");
            ib_mesh_writer.WriteVtkUsingMesh(mesh);

            OutputFileHandler handler("TestIbMeshWriterVTKCornerOverlap", false);
            std::string results_file = handler.GetOutputDirectoryFullPath() + "ib_mesh_vtk_corner_overlap.vtu";
            FileFinder vtk_file(results_file, RelativeTo::Absolute);
            TS_ASSERT(vtk_file.Exists());
        }
        // Overlap everything
        {
            // Generate nodes
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(1, true, 0.9, 0.1));
            nodes.push_back(new Node<2>(2, true, 0.9, 0.9));
            nodes.push_back(new Node<2>(3, true, 0.1, 0.9));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems);

            ImmersedBoundaryMeshWriter<2, 2> ib_mesh_writer("TestIbMeshWriterVTKCornerOverlap", "ib_mesh_vtk_corner_overlap");
            ib_mesh_writer.WriteVtkUsingMesh(mesh);

            OutputFileHandler handler("TestIbMeshWriterVTKCornerOverlap", false);
            std::string results_file = handler.GetOutputDirectoryFullPath() + "ib_mesh_vtk_corner_overlap.vtu";
            FileFinder vtk_file(results_file, RelativeTo::Absolute);
            TS_ASSERT(vtk_file.Exists());
        }
#else
        std::cout << "This test ran, but did not test VTK-dependent functions as VTK visualization is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }
};

#endif /*TESTIMMERSEDBOUNDARYMESHWRITER_HPP_*/
