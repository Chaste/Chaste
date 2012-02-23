/*

Copyright (c) 2005-2012, University of Oxford.
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
