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

#ifndef TESTMAJORAIRWAYSCENTRELINESCLEANER_HPP_
#define TESTMAJORAIRWAYSCENTRELINESCLEANER_HPP_

#include <cxxtest/TestSuite.h>
#include "OutputFileHandler.hpp"
#include "TetrahedralMesh.hpp"
#include "MutableMesh.hpp"
#include "MajorAirwaysCentreLinesCleaner.hpp"
#include "TrianglesMeshWriter.hpp"
#include "VtkMeshWriter.hpp"


//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestMajorAirwaysCentreLinesCleaner : public CxxTest::TestSuite
{
public:
    void TestDeleteOrderSimpleMesh()
    {
        TrianglesMeshReader<1,3> mesh_reader("lung/test/airway_generation/data/test_major_airways_mesh");
        MutableMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        //Assign valid radii
        for (unsigned node_index = 0; node_index < mesh.GetNumNodes(); ++node_index)
        {
            mesh.GetNode(node_index)->rGetNodeAttributes()[0] = 1.0;
        }

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 5u);

        MajorAirwaysCentreLinesCleaner cleaner(mesh, 0u);
        cleaner.CleanUsingHorsfieldOrder(1u);

        NodeMap node_map(mesh.GetNumAllNodes());
        mesh.ReIndex(node_map);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetLocation()[2], -2.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetLocation()[2], 0.0, 1e-6);
    }

    void TestDeleteFirstOrder()
    {
#ifdef CHASTE_VTK
        VtkMeshReader<1,3> mesh_reader("lung/test/data/TestSubject002MajorAirways.vtu");

        MutableMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        //Assign valid radii
        for (unsigned node_index = 0; node_index < mesh.GetNumNodes(); ++node_index)
        {
            mesh.GetNode(node_index)->AddNodeAttribute(1.0);
        }

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 12065u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 12064u);

        MajorAirwaysCentreLinesCleaner cleaner(mesh, 0u);
        cleaner.CleanUsingHorsfieldOrder(1u);

        NodeMap node_map(mesh.GetNumAllNodes());
        mesh.ReIndex(node_map);

        //Confirmed visually to be correct
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 3683u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3682u);

// Uncomment to visualise
//        VtkMeshWriter<1,3> mesh_writer("TestMajorAirwaysCentreLinesCleaner", "Novartis002Trimmed");
//        mesh_writer.WriteFilesUsingMesh(mesh);
#endif //CHASTE_VTK
    }

    void TestHeuristicCleanSimpleMesh()
    {
        {
            TrianglesMeshReader<1,3> mesh_reader("lung/test/airway_generation/data/test_major_airways_mesh");
            MutableMesh<1,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            //Assign valid radii
            for (unsigned node_index = 0; node_index < mesh.GetNumNodes(); ++node_index)
            {
                 mesh.GetNode(node_index)->rGetNodeAttributes()[0] = 1.0;
            }

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);
            TS_ASSERT_EQUALS(mesh.GetNumElements(), 5u);

            MajorAirwaysCentreLinesCleaner cleaner(mesh, 0u);
            cleaner.CleanTerminalsHueristic();

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);
            TS_ASSERT_EQUALS(mesh.GetNumElements(), 5u);

            //Trips an added nodes exception. Not sure if this is needed!
            //NodeMap node_map(mesh.GetNumAllNodes());
            //mesh.ReIndex(node_map);

            TS_ASSERT_DELTA(mesh.GetNode(2)->rGetLocation()[0], 1.6, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(3)->rGetLocation()[0], -1.6, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[1], 1.6, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[1], -1.6, 1e-6);


            //Uncomment to visualise
//            VtkMeshWriter<1,3> mesh_writer("TestMajorAirwaysCentreLinesCleaner", "heuristic_shorten", false);
//            mesh_writer.WriteFilesUsingMesh(mesh);
        }

        {
            TrianglesMeshReader<1,3> mesh_reader("lung/test/airway_generation/data/test_major_airways_mesh_short_terminals");
            MutableMesh<1,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            //Assign valid radii
            for (unsigned node_index = 0; node_index < mesh.GetNumNodes(); ++node_index)
            {
                 mesh.GetNode(node_index)->rGetNodeAttributes()[0] = 1.0;
            }

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);
            TS_ASSERT_EQUALS(mesh.GetNumElements(), 5u);

            MajorAirwaysCentreLinesCleaner cleaner(mesh, 0u);
            cleaner.CleanTerminalsHueristic();

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);
            TS_ASSERT_EQUALS(mesh.GetNumElements(), 5u);

            //Trips an added nodes exception. Not sure if this is needed!
            //NodeMap node_map(mesh.GetNumAllNodes());
            //mesh.ReIndex(node_map);

            TS_ASSERT_DELTA(mesh.GetNode(2)->rGetLocation()[0], 1.4, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(3)->rGetLocation()[0], -1.4, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[1], 1.4, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[1], -1.4, 1e-6);

            //Uncomment to visualise
//            VtkMeshWriter<1,3> mesh_writer("TestMajorAirwaysCentreLinesCleaner", "heuristic_lengthen", false);
//            mesh_writer.WriteFilesUsingMesh(mesh);
        }
    }

    void TestRemoveIsolatedNodesSimpleMesh()
    {
        TrianglesMeshReader<1,3> mesh_reader("lung/test/airway_generation/data/test_isolated_nodes_major_airways_mesh");
        MutableMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Assign valid radii
        for (unsigned node_index = 0; node_index < mesh.GetNumNodes(); ++node_index)
        {
          mesh.GetNode(node_index)->rGetNodeAttributes()[0] = 1.0;
        }

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);

        MajorAirwaysCentreLinesCleaner cleaner(mesh, 0u);
        cleaner.CleanIsolatedNodes();

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);
    }
};

#endif /* TESTMAJORAIRWAYSCENTRELINESCLEANER_HPP_ */
