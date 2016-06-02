/*

Copyright (c) 2005-2016, University of Oxford.
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

// Needed for test framework
#include <cxxtest/cxxtest/TestSuite.h>

#include "ImmersedBoundaryMesh.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundaryMesh : public CxxTest::TestSuite
{
public:

    void TestSolveNodeAndElementMapping() throw(Exception)
    {
    }

    void TestClear() throw(Exception)
    {

    }

    void TestSetupFluidVelocityGrids() throw(Exception)
    {

    }

    void TestArchiving() throw(Exception)
    {

    }

    void TestElementIterator() throw(Exception)
    {

    }

    void TestSetAndGetMethods() throw(Exception)
    {

    }

    void TestGetVectorFromAtoB() throw(Exception)
    {

    }

    void TestGetSkewnessOfElementMassDistributionAboutAxis() throw(Exception)
    {
        // A square should have no skewness about any axis
        {
            std::vector <Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector < ImmersedBoundaryElement < 2, 2 > * > elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2, 2> *p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

            for (unsigned i = 0 ; i < 16 ; i++)
            {
                double theta = 2.0 * M_PI * (double)i / 16.0;

                c_vector<double, 2> axis;
                axis[0] = cos(theta);
                axis[1] = sin(theta);

                TS_ASSERT_DELTA(p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(0, axis), 0.0, 1e-12);
            }

            delete(p_mesh);
        }

        // A triangle should have skewness
        {
            std::vector <Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));

            std::vector < ImmersedBoundaryElement < 2, 2 > * > elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2, 2> *p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

            c_vector<double, 2> axis;
            axis[0] = 0.0;
            axis[1] = 1.0;

            double hand_calculated_skewness = -0.5656854249;

            // Test that the skewness is equal to the hand calculated value
            TS_ASSERT_DELTA(p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(0, axis) - hand_calculated_skewness,
                            0.0, 1e-9);

            // If we flip the axis, the skewness should be minus what it was before
            axis[1] = -1.0;
            TS_ASSERT_DELTA(p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(0, axis) + hand_calculated_skewness,
                            0.0, 1e-9);
        }
    }
};
