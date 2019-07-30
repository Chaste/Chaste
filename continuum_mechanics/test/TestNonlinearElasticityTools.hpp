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

#ifndef TESTNONLINEARELASTICITYTOOLSGROUP_HPP_
#define TESTNONLINEARELASTICITYTOOLSGROUP_HPP_

#include <cxxtest/TestSuite.h>

#include "NonlinearElasticityTools.hpp"
#include "QuadraticMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestNonlinearElasticityTools : public CxxTest::TestSuite
{
public:
    void TestGetNodesByComponentValue()
    {
        QuadraticMesh<2> mesh(0.1,1.0,1.0);
        mesh.Scale(1.0, 2.0); //historical reasons

        std::vector<unsigned> indices
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        TS_ASSERT_EQUALS(indices.size(), 21u);
        for (unsigned i=0; i<indices.size(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(indices[i])->rGetLocation()[0], 0.0, 1e-12);
        }

        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> mesh3d;
        mesh3d.ConstructFromMeshReader(reader);
        mesh3d.Scale(0.3423244,1.343244325,6.23435);

        std::vector<unsigned> indices3d
          = NonlinearElasticityTools<3>::GetNodesByComponentValue(mesh3d,2,6.23435);

        TS_ASSERT_EQUALS(indices3d.size(), 13u);
        for (unsigned i=0; i<indices3d.size(); i++)
        {
            TS_ASSERT_DELTA(mesh3d.GetNode(indices3d[i])->rGetLocation()[2], 6.23435, 1e-12);
        }

        TS_ASSERT_THROWS_CONTAINS(NonlinearElasticityTools<3>::GetNodesByComponentValue(mesh3d, 2, 6.234),
                "Could not find any nodes on requested surface (note: tolerance = 1e-0");
    }
};

#endif /*TESTNONLINEARELASTICITYTOOLSGROUP_HPP_*/
