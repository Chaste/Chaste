/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef TESTSEMMESHGENERATOR_HPP_
#define TESTSEMMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "SemMeshGenerator.hpp"

class TestSemMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestConstructorExceptions() throw (Exception)
    {
        TS_ASSERT_THROWS_THIS(SemMeshGenerator<2> mesh_generator(0, 2, 0, 0), "An SemMeshGenerator must have numCellsAcross > 0");
        TS_ASSERT_THROWS_THIS(SemMeshGenerator<2> mesh_generator(2, 2, 0, 0), "An SemMeshGenerator must have more than 1 subcellular element per cell. Set numSubCellularElementsPerCell > 0");
        TS_ASSERT_THROWS_THIS(SemMeshGenerator<1> mesh_generator(2, 2, 2), "Trying to create a 1D mesh, SemMesh only defined for 2 and 3D");
        TS_ASSERT_THROWS_THIS(SemMeshGenerator<2> mesh_generator(2, 2, 2), "Trying to create a 3D SemMesh with DIM < 3");
    }

    void TestConstruct2dLinearMesh() throw (Exception)
    {
        // Create a string of 5 cells.
        SemMeshGenerator<2> generator(5);
        SemMesh<2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 500u);

        TS_ASSERT_DELTA(generator.GetEquilibriumDistance(), 0.1904, 1e-4);

        // Check some nodes are in the right place
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[1], 0.0, 1e-6);

        TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(99)->rGetLocation()[0], 9.0*generator.GetEquilibriumDistance(), 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(99)->rGetLocation()[1], 9.0*generator.GetEquilibriumDistance(), 1e-4);
    }

    void TestConstruct2dSquareMesh() throw (Exception)
    {
        // Create a string of 5 cells.
        SemMeshGenerator<2> generator(5, 5);
        SemMesh<2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 25u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 2500u);

        // Check some nodes are in the right place
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[1], 0.0, 1e-6);

        TS_ASSERT_DELTA(p_mesh->GetElement(24)->GetNode(99)->rGetLocation()[0], 49.0*generator.GetEquilibriumDistance(), 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetElement(24)->GetNode(99)->rGetLocation()[1], 49.0*generator.GetEquilibriumDistance(), 1e-4);
    }

    void TestConstruct3dCuboidMesh() throw (Exception)
    {
        // Create a string of 5 cells.
        SemMeshGenerator<3> generator(5, 5, 5, 5, 5, 5);
        SemMesh<3>* p_mesh = generator.GetMesh();

        TS_ASSERT_DELTA(generator.GetEquilibriumDistance(), 0.3619, 1e-4);

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 125u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 15625u);

        // Check some nodes are in the right place
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[2], 0.0, 1e-6);

        // Check some nodes are in the right place
        TS_ASSERT_DELTA(p_mesh->GetElement(124)->GetNode(124)->rGetLocation()[0], 24.0*generator.GetEquilibriumDistance(), 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetElement(124)->GetNode(124)->rGetLocation()[1], 24.0*generator.GetEquilibriumDistance(), 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetElement(124)->GetNode(124)->rGetLocation()[2], 24.0*generator.GetEquilibriumDistance(), 1e-6);
    }
};

#endif /*TESTSEMMESHGENERATOR_HPP_*/
