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



#ifndef TESTDEFORMEDBOUNDARYELEMENT_HPP_
#define TESTDEFORMEDBOUNDARYELEMENT_HPP_

#include "DeformedBoundaryElement.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "UblasCustomFunctions.hpp"

class TestDeformedBoundaryElement : public CxxTest::TestSuite
{
public:
    void TestDeformedBoundaryElement2d() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        DeformedBoundaryElement<1,2> deformed_bdy_element;

        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        std::vector<c_vector<double,2> > displacement(2); // 2=num nodes in the element
        displacement[0] = zero_vector<double>(2);
        displacement[1] = zero_vector<double>(2);

        //// This boundary element it seems has nodes (0,1) and (0,0)
        //std::cout << (*iter)->GetNode(0)->rGetLocation() << "\n";
        //std::cout << (*iter)->GetNode(1)->rGetLocation() << "\n";

        // specify the undeformed element and the (zero) displacement
        deformed_bdy_element.ApplyUndeformedElementAndDisplacement(*iter, displacement);

        TS_ASSERT_DELTA( deformed_bdy_element.GetNode(0)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA( deformed_bdy_element.GetNode(0)->rGetLocation()[1], 1.0, 1e-12);
        TS_ASSERT_DELTA( deformed_bdy_element.GetNode(1)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA( deformed_bdy_element.GetNode(1)->rGetLocation()[1], 0.0, 1e-12);


        // normal = (-1, 0)
        TS_ASSERT_DELTA(deformed_bdy_element.CalculateNormal()(0),-1.0, 1e-12);
        TS_ASSERT_DELTA(deformed_bdy_element.CalculateNormal()(1), 0.0, 1e-12);

        // move nodes to (1.1,1.1) and (0.1,0.1)
        displacement[0](0) = 1.01;
        displacement[0](1) = 0.01;

        displacement[1](0) = 0.01;
        displacement[1](1) = 0.01;

        deformed_bdy_element.ApplyUndeformedElementAndDisplacement(*iter, displacement);

        TS_ASSERT_DELTA( deformed_bdy_element.GetNode(0)->rGetLocation()[0], 1.01, 1e-12);
        TS_ASSERT_DELTA( deformed_bdy_element.GetNode(0)->rGetLocation()[1], 1.01, 1e-12);
        TS_ASSERT_DELTA( deformed_bdy_element.GetNode(1)->rGetLocation()[0], 0.01, 1e-12);
        TS_ASSERT_DELTA( deformed_bdy_element.GetNode(1)->rGetLocation()[1], 0.01, 1e-12);

        // normal is now (-1/sqrt(2), 1/sqrt(2))
        TS_ASSERT_DELTA(deformed_bdy_element.CalculateNormal()(0),-1.0/sqrt(2), 1e-12);
        TS_ASSERT_DELTA(deformed_bdy_element.CalculateNormal()(1), 1.0/sqrt(2), 1e-12);
    }


    void TestDeformedBoundaryElement3d() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        DeformedBoundaryElement<2,3> deformed_bdy_element;

        TetrahedralMesh<3,3>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();

        std::vector<c_vector<double,3> > displacement(3); // 2=num nodes in the element
        displacement[0] = zero_vector<double>(3);
        displacement[1] = zero_vector<double>(3);
        displacement[2] = zero_vector<double>(3);

        //// This boundary element it seems has nodes (0.5,0.5,1), (0,0.5,1) and (0.25,0.75,1)
        //std::cout << (*iter)->GetNode(0)->rGetLocation() << "\n";
        //std::cout << (*iter)->GetNode(1)->rGetLocation() << "\n";
        //std::cout << (*iter)->GetNode(2)->rGetLocation() << "\n";

        // specify the undeformed element and the (zero) displacement
        deformed_bdy_element.ApplyUndeformedElementAndDisplacement(*iter, displacement);

        for(unsigned i=0; i<3; i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA( deformed_bdy_element.GetNode(i)->rGetLocation()[j], (*iter)->GetNode(i)->rGetLocation()[j], 1e-12);
            }
        }

        // normal = (0,0,1)
        TS_ASSERT_DELTA(deformed_bdy_element.CalculateNormal()(0), 0.0, 1e-12);
        TS_ASSERT_DELTA(deformed_bdy_element.CalculateNormal()(1), 0.0, 1e-12);
        TS_ASSERT_DELTA(deformed_bdy_element.CalculateNormal()(2), 1.0, 1e-12);
    }

};


#endif /* TESTDEFORMEDBOUNDARYELEMENT_HPP_ */
