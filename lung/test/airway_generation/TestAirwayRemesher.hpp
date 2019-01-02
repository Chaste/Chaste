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

#ifndef TESTAIRWAYREMESHER_HPP_
#define TESTAIRWAYREMESHER_HPP_

#include <cxxtest/TestSuite.h>
#include "OutputFileHandler.hpp"
#include "TetrahedralMesh.hpp"
#include "MutableMesh.hpp"
#include "AirwayRemesher.hpp"
#include "TrianglesMeshWriter.hpp"
#include "VtkMeshWriter.hpp"


//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestAirwayRemesher : public CxxTest::TestSuite
{
public:
    void TestRemeshSingleBranch()
    {
        //Load a single branch mesh file
        TrianglesMeshReader<1,3> reader("mesh/test/data/1D_in_3D_0_to_1mm_10_elements");
        TetrahedralMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        //We need to add some attributes to the mesh
        for (TetrahedralMesh<1,3>::NodeIterator iter = mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            iter->AddNodeAttribute(0.05);
        }

        //Create remesher object
        AirwayRemesher remesher(mesh, 0u);

        //Check intermediate elements are removed. Poiseuille resistance of the branch is
        // C*0.1/((0.05^4) = C*1.6 * 10^4

        MutableMesh<1,3> output_mesh_one;
        remesher.Remesh(output_mesh_one, 1e5); //With this tolerance all intermediate nodes should be removed.

        TS_ASSERT_EQUALS(output_mesh_one.GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(output_mesh_one.GetNumElements(), 1u);
        TS_ASSERT_DELTA(output_mesh_one.GetElement(0)->GetAttribute(), 0.05, 1e-6);

        MutableMesh<1,3> output_mesh_two;
        remesher.Remesh(output_mesh_two, 0.8e4); //With this tolerance there should be one intermediate node.

        TS_ASSERT_EQUALS(output_mesh_two.GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(output_mesh_two.GetNumElements(), 2u);
        TS_ASSERT_DELTA(output_mesh_two.GetElement(0)->GetAttribute(), 0.05, 1e-6);

        MutableMesh<1,3> output_mesh_three;
        remesher.Remesh(output_mesh_three, 1.6e3); //With this tolerance there should be ten elements.

        TS_ASSERT_EQUALS(output_mesh_three.GetNumNodes(), 11u);
        TS_ASSERT_EQUALS(output_mesh_three.GetNumElements(), 10u);
        TS_ASSERT_DELTA(output_mesh_three.GetElement(0)->GetAttribute(), 0.05, 1e-6);
        TS_ASSERT_DELTA(output_mesh_three.GetElement(5)->GetAttribute(), 0.05, 1e-6);

        MutableMesh<1,3> output_mesh_four;
        remesher.Remesh(output_mesh_four); //Utility method to remove all intermediate nodes

        TS_ASSERT_EQUALS(output_mesh_one.GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(output_mesh_one.GetNumElements(), 1u);
        TS_ASSERT_DELTA(output_mesh_one.GetElement(0)->GetAttribute(), 0.05, 1e-6);

        //To visualise
        //VtkMeshWriter<1,3> writer("TestAirwayRemesher", "1D_remeshed");
        //writer.WriteFilesUsingMesh(output_mesh_three);
    }

    void TestRemeshFullTree()
    {
        TrianglesMeshReader<1,3> reader("lung/test/data/TestSubject002");

        TetrahedralMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        AirwayPropertiesCalculator calculator(mesh, 0);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 136625u);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 136624u);

        //Create remesher object
        AirwayRemesher remesher(mesh, 0u);

        MutableMesh<1,3> output_mesh_one;
        remesher.Remesh(output_mesh_one, calculator.GetBranches()[0]->GetPoiseuilleResistance()*1e7); //Key the tolerance relative to the trachea

        TS_ASSERT_EQUALS( output_mesh_one.GetNumNodes(), 168045u);
        TS_ASSERT_EQUALS( output_mesh_one.GetNumElements(), 168044u);

//        //To visualise
//
//        VtkMeshWriter<1,3> writer("TestAirwayRemesher", "Novartis002_remeshed");
//        std::vector<double> radii(output_mesh_one.GetNumElements());
//
//        for (TetrahedralMesh<1,3>::ElementIterator iter = output_mesh_one.GetElementIteratorBegin();
//            iter != output_mesh_one.GetElementIteratorEnd();
//            ++iter)
//        {
//            radii[iter->GetIndex()] = iter->GetAttribute();
//        }
//
//        writer.AddCellData("radii", radii);
//        writer.WriteFilesUsingMesh(output_mesh_one);
//
//        TrianglesMeshWriter<1,3> writer2("TestAirwayRemesher", "Novartis002_remeshed", false);
//        writer2.WriteFilesUsingMesh(output_mesh_one);
    }
};

#endif /* TESTAIRWAYREMESHER_HPP_ */
