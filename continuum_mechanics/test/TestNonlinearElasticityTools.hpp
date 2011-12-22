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

#ifndef TESTNONLINEARELASTICITYTOOLSGROUP_HPP_
#define TESTNONLINEARELASTICITYTOOLSGROUP_HPP_

#include <cxxtest/TestSuite.h>

#include "NonlinearElasticityTools.hpp"
#include "QuadraticMesh.hpp"
#include "TrianglesMeshReader.hpp"

class TestNonlinearElasticityTools : public CxxTest::TestSuite
{
public:
    void TestGetNodesByComponentValue() throw(Exception)
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

        TS_ASSERT_THROWS_THIS(NonlinearElasticityTools<3>::GetNodesByComponentValue(mesh3d,2,6.234),
                "Could not find any nodes on requested surface (note: tolerance = 1e-08)");
    }
};

#endif /*TESTNONLINEARELASTICITYTOOLSGROUP_HPP_*/
