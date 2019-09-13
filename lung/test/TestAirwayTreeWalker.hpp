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

#ifndef _TESTAIRWAYTREEWALKER_HPP_
#define _TESTAIRWAYTREEWALKER_HPP_

#include <cxxtest/TestSuite.h>
#include <queue>

#include "AirwayTreeWalker.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

class TestAirwayTreeWalker : public CxxTest::TestSuite
{
public:
    void TestGetMethods()
    {
        TetrahedralMesh<1, 3> mesh;
        TrianglesMeshReader<1, 3> mesh_reader("mesh/test/data/three_generation_branch_mesh_refined");
        mesh.ConstructFromMeshReader(mesh_reader);

        AirwayTreeWalker walker(mesh, 0u);

        TS_ASSERT_EQUALS(walker.GetOutletElementIndex(), 0u);
        TS_ASSERT_EQUALS(walker.GetNodesAreGraphOrdered(), true);

        TS_ASSERT_EQUALS(walker.GetParentElement(mesh.GetElement(2u)), mesh.GetElement(1u));
        TS_ASSERT_EQUALS(walker.GetParentElementIndex(mesh.GetElement(2u)), 1u);
        TS_ASSERT_EQUALS(walker.GetParentElement(2u), mesh.GetElement(1u));
        TS_ASSERT_EQUALS(walker.GetParentElementIndex(2u), 1u);

        TS_ASSERT(walker.GetParentElement(mesh.GetElement(0u)) == NULL);
        TS_ASSERT(walker.GetParentElement(0u) == NULL);

        std::vector<Element<1, 3>*> child_eles = walker.GetChildElements(mesh.GetElement(1u));
        TS_ASSERT_EQUALS(child_eles.size(), 2u);
        TS_ASSERT_EQUALS(child_eles[0], mesh.GetElement(2u));
        TS_ASSERT_EQUALS(child_eles[1], mesh.GetElement(6u));

        std::vector<unsigned> child_ele_indices = walker.GetChildElementIndices(mesh.GetElement(1u));
        TS_ASSERT_EQUALS(child_ele_indices.size(), 2u);
        TS_ASSERT_EQUALS(child_ele_indices[0], 2u);
        TS_ASSERT_EQUALS(child_ele_indices[1], 6u);

        TS_ASSERT_EQUALS(walker.GetNumberOfChildElements(1u), 2u);

        TS_ASSERT_EQUALS(walker.GetDistalNode(mesh.GetElement(0u)), mesh.GetNode(1u));
        TS_ASSERT_EQUALS(walker.GetDistalNodeIndex(mesh.GetElement(0u)), 1u);
        TS_ASSERT_EQUALS(walker.GetDistalNode(mesh.GetElement(9u)), mesh.GetNode(9u));
        TS_ASSERT_EQUALS(walker.GetDistalNodeIndex(mesh.GetElement(9u)), 9u);
    }

    void TestGenerations()
    {
        TetrahedralMesh<1, 3> mesh;
        TrianglesMeshReader<1, 3> mesh_reader("mesh/test/data/three_generation_branch_mesh_refined");
        mesh.ConstructFromMeshReader(mesh_reader);

        AirwayTreeWalker walker(mesh, 0u);

        TS_ASSERT_EQUALS(walker.GetElementGeneration(0u), 0u);
        TS_ASSERT_EQUALS(walker.GetElementGeneration(1u), 0u);
        TS_ASSERT_EQUALS(walker.GetElementGeneration(6u), 1u);
        TS_ASSERT_EQUALS(walker.GetElementGeneration(9u), 2u);

        TS_ASSERT_EQUALS(walker.GetElementGeneration(mesh.GetElement(0u)), 0u);
        TS_ASSERT_EQUALS(walker.GetElementGeneration(mesh.GetElement(1u)), 0u);
        TS_ASSERT_EQUALS(walker.GetElementGeneration(mesh.GetElement(6u)), 1u);
        TS_ASSERT_EQUALS(walker.GetElementGeneration(mesh.GetElement(9u)), 2u);

        TS_ASSERT_EQUALS(walker.GetMaxElementGeneration(), 2u);
    }

    void TestOrders()
    {
        TetrahedralMesh<1, 3> mesh;
        TrianglesMeshReader<1, 3> mesh_reader("mesh/test/data/three_generation_branch_mesh_refined");
        mesh.ConstructFromMeshReader(mesh_reader);

        AirwayTreeWalker walker(mesh, 0u);

        // mesh is symmetric, so Horsfield and Strahler orders are equivalent
        TS_ASSERT_EQUALS(walker.GetElementHorsfieldOrder(0u), 3u);
        TS_ASSERT_EQUALS(walker.GetElementHorsfieldOrder(1u), 3u);
        TS_ASSERT_EQUALS(walker.GetElementHorsfieldOrder(6u), 2u);
        TS_ASSERT_EQUALS(walker.GetElementHorsfieldOrder(9u), 1u);

        TS_ASSERT_EQUALS(walker.GetElementHorsfieldOrder(mesh.GetElement(0u)), 3u);
        TS_ASSERT_EQUALS(walker.GetElementHorsfieldOrder(mesh.GetElement(1u)), 3u);
        TS_ASSERT_EQUALS(walker.GetElementHorsfieldOrder(mesh.GetElement(6u)), 2u);
        TS_ASSERT_EQUALS(walker.GetElementHorsfieldOrder(mesh.GetElement(9u)), 1u);

        TS_ASSERT_EQUALS(walker.GetMaxElementHorsfieldOrder(), 3u);

        TS_ASSERT_EQUALS(walker.GetElementStrahlerOrder(0u), 3u);
        TS_ASSERT_EQUALS(walker.GetElementStrahlerOrder(1u), 3u);
        TS_ASSERT_EQUALS(walker.GetElementStrahlerOrder(6u), 2u);
        TS_ASSERT_EQUALS(walker.GetElementStrahlerOrder(9u), 1u);

        TS_ASSERT_EQUALS(walker.GetElementStrahlerOrder(mesh.GetElement(0u)), 3u);
        TS_ASSERT_EQUALS(walker.GetElementStrahlerOrder(mesh.GetElement(1u)), 3u);
        TS_ASSERT_EQUALS(walker.GetElementStrahlerOrder(mesh.GetElement(6u)), 2u);
        TS_ASSERT_EQUALS(walker.GetElementStrahlerOrder(mesh.GetElement(9u)), 1u);

        TS_ASSERT_EQUALS(walker.GetMaxElementStrahlerOrder(), 3u);
    }
};

#endif /*_TESTAIRWAYTREEWALKER_HPP_*/
