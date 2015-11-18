/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef _TESTAIRWAYWALLMODELS_HPP_
#define _TESTAIRWAYWALLMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include "LambertAirwayWall.hpp"
#include "LambertAirwayWallFactory.hpp"

//#include "PetscSetupAndFinalize.hpp"

class TestAirwayWallModels: public CxxTest::TestSuite
{
public:

    void TestLambertAirwayWallAndFactory() throw (Exception)
    {
        //Get a simple mesh here
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> reader("lung/test/data/three_bifurcations");
        mesh.ConstructFromMeshReader(reader);

        {
            LambertAirwayWallFactory factory;
            factory.SetMesh(&mesh);

            //The three bifurcation mesh should map onto equivalent generation 0, generation 8 and generation 16 Lambert airways
            { //generation 0
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(0u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.882, 1e-3);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 0.5, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 10.0, 1e-6);
            }

            { //generation 2
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(2u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.213, 1e-3);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 1.0, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 10.0, 1e-6);
            }

            { //generation 4
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(4u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.039, 1e-6);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 1.0, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 7.0, 1e-6);
            }

            for(unsigned element_index = 0; element_index < mesh.GetNumElements(); ++element_index)
            {
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(element_index));

                p_wall->SetAirwayPressure(0.0);
                p_wall->SetPleuralPressure(0.0);

                p_wall->SolveAndUpdateState(0.0, 0.0);
                TS_ASSERT_DELTA(p_wall->GetLumenRadius(), p_wall->mRi, 1e-6);

                p_wall->SetPleuralPressure(-500000000);
                p_wall->SolveAndUpdateState(0.0, 0.0);
                TS_ASSERT_DELTA(p_wall->GetLumenRadius(), p_wall->mRiMax, 1e-6);

                p_wall->SetPleuralPressure(50000000000);
                p_wall->SolveAndUpdateState(0.0, 0.0);
                TS_ASSERT_DELTA(p_wall->GetLumenRadius(), 0.0, 1e-3);
            }
        }

        //Repeat using Strahler order (equivalent for this mesh)
        {
            LambertAirwayWallFactory factory(true);
            factory.SetMesh(&mesh);

            //The three bifurcation mesh should map onto equivalent generation 0, generation 8 and generation 16 Lambert airways
            { //generation 0
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(0u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.882, 1e-3);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 0.5, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 10.0, 1e-6);
            }

            { //generation 2
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(2u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.213, 1e-3);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 1.0, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 10.0, 1e-6);
            }

            { //generation 4
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(4u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.039, 1e-6);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 1.0, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 7.0, 1e-6);
            }

            for(unsigned element_index = 0; element_index < mesh.GetNumElements(); ++element_index)
            {
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(element_index));

                p_wall->SetAirwayPressure(0.0);
                p_wall->SetPleuralPressure(0.0);

                p_wall->SolveAndUpdateState(0.0, 0.0);
                TS_ASSERT_DELTA(p_wall->GetLumenRadius(), p_wall->mRi, 1e-6);

                p_wall->SetPleuralPressure(-500000000);
                p_wall->SolveAndUpdateState(0.0, 0.0);
                TS_ASSERT_DELTA(p_wall->GetLumenRadius(), p_wall->mRiMax, 1e-6);

                p_wall->SetPleuralPressure(50000000000);
                p_wall->SolveAndUpdateState(0.0, 0.0);
                TS_ASSERT_DELTA(p_wall->GetLumenRadius(), 0.0, 1e-3);
            }
        }
    }
};
#endif /*_TESTAIRWAYWALLMODELS_HPP_*/

