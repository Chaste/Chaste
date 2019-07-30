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

#ifndef _TESTAIRWAYWALLMODELS_HPP_
#define _TESTAIRWAYWALLMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include "LambertAirwayWall.hpp"
#include "LambertAirwayWallFactory.hpp"
#include "LaPradAirwayWall.hpp"
#include "LaPradAirwayWallFactory.hpp"
#include "HiornsAirwayWall.hpp"
#include "HiornsAirwayWallFactory.hpp"

class TestAirwayWallModels: public CxxTest::TestSuite
{
public:

    void TestLambertAirwayWallAndFactory()
    {
        // Get a simple mesh here
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> reader("lung/test/data/three_bifurcations");
        mesh.ConstructFromMeshReader(reader);

        {
            LambertAirwayWallFactory factory;
            factory.SetMesh(&mesh);

            // Coverage
            TS_ASSERT_EQUALS(factory.GetNumberOfAirways(), 7u);
            TS_ASSERT_DELTA(factory.GetAlpha0ForGeneration(0), 0.882, 1e-6);
            TS_ASSERT_DELTA(factory.GetAlpha0PrimeForGeneration(0), 0.011/98.0665, 1e-6);
            TS_ASSERT_DELTA(factory.GetN1ForGeneration(0), 0.5, 1e-6);
            TS_ASSERT_DELTA(factory.GetN2ForGeneration(0), 10.0, 1e-6);
            TS_ASSERT_DELTA(factory.GetPleuralPressureForAirway(0, NULL), 0.0, 1e-6);
            TS_ASSERT_EQUALS(factory.GetMesh()->GetNumElements(), 7u);

            // The three bifurcation mesh should map onto equivalent generation 0, generation 8 and generation 16 Lambert airways
            { // generation 0
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(0u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.882, 1e-3);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 0.5, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 10.0, 1e-6);
                delete p_wall;
            }

            { // generation 8
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(2u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.213, 1e-3);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 1.0, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 10.0, 1e-6);
                delete p_wall;
            }

            { //generation 16
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(4u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.039, 1e-6);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 1.0, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 7.0, 1e-6);
                delete p_wall;
            }

            for (unsigned element_index = 0; element_index < mesh.GetNumElements(); ++element_index)
            {
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(element_index));

                p_wall->SetTimestep(0.0); //For coverage, Lambert is a quasi-static model
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

                delete p_wall;
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
                delete p_wall;
            }

            { //generation 2
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(2u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.213, 1e-3);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 1.0, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 10.0, 1e-6);
                delete p_wall;
            }

            { //generation 4
                LambertAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(4u));
                double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(alpha0, 0.039, 1e-6);
                TS_ASSERT_DELTA(p_wall->mRiMax, 0.05, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN1, 1.0, 1e-6);
                TS_ASSERT_DELTA(p_wall->mN2, 7.0, 1e-6);
                delete p_wall;
            }

            for (unsigned element_index = 0; element_index < mesh.GetNumElements(); ++element_index)
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

                delete p_wall;
            }
        }
    }

    void TestLaPradAirwayWall()
    {
        LaPradAirwayWall airway_wall;
        airway_wall.SetTimestep(1.0); //Not used, for coverage only

        double targetPressure = -3.;
        double RIn = 2.;
        double ROut = 2.5;
        double k1 = 1.8;
        double k2 = 0.5;
        double k3 = 10.;

        airway_wall.SetAirwayPressure(0.0);
        airway_wall.SetPleuralPressure(targetPressure);
        airway_wall.SetRIn(RIn);
        airway_wall.SetROut(ROut);
        airway_wall.Setk1(k1);
        airway_wall.Setk2(k2);
        airway_wall.Setk3(k3);

        //Validated against Matlab implementation
        airway_wall.SolveAndUpdateState(0.0, 0.0);
        TS_ASSERT_DELTA(airway_wall.GetLumenRadius(), 3.02842, 1e-3)

        //Validated against Matlab implementation, consistent with plots in paper
        airway_wall.SolveAndUpdateState(0.0, 0.0);
        TS_ASSERT_DELTA(airway_wall.GetLumenRadius(), 3.02842, 1e-3)

        airway_wall.SetPleuralPressure(-1000000);
        airway_wall.SolveAndUpdateState(0.0, 0.0);

        TS_ASSERT_DELTA(airway_wall.GetLumenRadius(), 3.4917, 1e-3)

        airway_wall.SetPleuralPressure(0.0);
        airway_wall.SolveAndUpdateState(0.0, 0.0);
        TS_ASSERT_DELTA(airway_wall.GetLumenRadius(), RIn, 1e-3)

    }


    void TestLaPradAirwayWallAndFactory()
    {

         //Get a simple mesh here
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> reader("lung/test/data/three_bifurcations");
        mesh.ConstructFromMeshReader(reader);

        {
            LaPradAirwayWallFactory factory;
            factory.SetMesh(&mesh);

            // Coverage
            TS_ASSERT_DELTA(factory.Getk1ForGeneration(0), 4000.0, 1e-6);
            TS_ASSERT_DELTA(factory.Getk2ForGeneration(0), 1000.0, 1e-6);
            TS_ASSERT_DELTA(factory.Getk3ForGeneration(0), 20.0, 1e-6);
            TS_ASSERT_DELTA(factory.GetAlpha0ForGeneration(0), 0.882, 1e-6);
            TS_ASSERT_DELTA(factory.GetPleuralPressureForAirway(0, NULL), 0.0, 1e-6);

            //The three bifurcation mesh should map onto equivalent generation 0, generation 8 and generation 16 Lambert airways
            { //generation 0
                LaPradAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(0u));
                TS_ASSERT_DELTA(p_wall->mk1, 4000., 1e-6);
                delete p_wall;
            }

            { //generation 8
                LaPradAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(2u));
                TS_ASSERT_DELTA(p_wall->mk3, 15.2, 1e-6);
                delete p_wall;
            }

            { //generation 16
                LaPradAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(4u));
                //double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(p_wall->mk2, 360., 1e-6);
                delete p_wall;
            }

            for (unsigned element_index = 0; element_index < mesh.GetNumElements(); ++element_index)
            {
                LaPradAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(element_index));

                p_wall->SetTimestep(0.0);
                p_wall->SetAirwayPressure(0.0);
                p_wall->SetPleuralPressure(0.0);

                p_wall->SolveAndUpdateState(0.0, 0.0);
                TS_ASSERT_DELTA(p_wall->GetLumenRadius(), p_wall->mRIn, 1e-3);

                delete p_wall;
            }

            }


            //Repeat using Strahler order (equivalent for this mesh)
        {
            LaPradAirwayWallFactory factory(true);
            factory.SetMesh(&mesh);

            //The three bifurcation mesh should map onto equivalent generation 0, generation 8 and generation 16 Lambert airways
            { //generation 0
                LaPradAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(0u));
                //double alpha0 = p_wall->mRi*p_wall->mRi/p_wall->mRiMax/p_wall->mRiMax;
                TS_ASSERT_DELTA(p_wall->mk1, 4000., 1e-6);
                delete p_wall;
            }

            { //generation 2
                LaPradAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(2u));
                TS_ASSERT_DELTA(p_wall->mk2, 680., 1e-6);
                delete p_wall;
            }

            { //generation 4
                LaPradAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(4u));
                TS_ASSERT_DELTA(p_wall->mk3, 10.4, 1e-6);
                delete p_wall;
            }

            for (unsigned element_index = 0; element_index < mesh.GetNumElements(); ++element_index)
            {
                LaPradAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(element_index));

                p_wall->SetTimestep(0.0);
                p_wall->SetAirwayPressure(0.0);
                p_wall->SetPleuralPressure(0.0);

                p_wall->SolveAndUpdateState(0.0, 0.0);

                TS_ASSERT_DELTA(p_wall->GetLumenRadius(), p_wall->mRIn, 1e-3);

                delete p_wall;
            }
        }

    }


    void TestHiornsAirwayWall()
    {

        HiornsAirwayWall airway_wall;
        airway_wall.SetTimestep(1.0); //Not used, for coverage only

        double targetPressure = -2.;
        double RIn = 1.;
        double ROut = 1.5;
        double mu = 5.;
        double phi1 = 0.;
        double phi2 = 0.;
        double C1 = 1.;
        double C2 = 1.8;
        double A = 10.;

        airway_wall.SetAirwayPressure(0.0);
        airway_wall.SetPleuralPressure(targetPressure);
        airway_wall.SetRIn(RIn);
        airway_wall.SetROut(ROut);
        airway_wall.Setmu(mu);
        airway_wall.Setphi1(phi1);
        airway_wall.Setphi2(phi2);
        airway_wall.SetC1(C1);
        airway_wall.SetC2(C2);
        airway_wall.SetA(A);

        //Validated against Matlab implementation
        airway_wall.SolveAndUpdateState(0.0, 0.0);
        TS_ASSERT_DELTA(airway_wall.GetLumenRadius(), 0.224008, 1e-3)

        airway_wall.SetA(1.);

        airway_wall.SetPleuralPressure(0.0);
        airway_wall.SolveAndUpdateState(0.0, 0.0);
        TS_ASSERT_DELTA(airway_wall.GetLumenRadius(), 0.769301, 1e-3)

    }


    void TestHiornsAirwayWallAndFactory()
    {
        // Get a simple mesh here
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> reader("lung/test/data/three_bifurcations");
        mesh.ConstructFromMeshReader(reader);

        {
            HiornsAirwayWallFactory factory;
            factory.SetMesh(&mesh);

            // Coverage
            TS_ASSERT_DELTA(factory.GetmuForGeneration(0), 64002.0, 1e-6);
            TS_ASSERT_DELTA(factory.Getphi1ForGeneration(0), 0.0, 1e-6);
            TS_ASSERT_DELTA(factory.Getphi2ForGeneration(0), 0.7854, 1e-6);
            TS_ASSERT_DELTA(factory.GetC1ForGeneration(0), 179380.0, 1e-6);
            TS_ASSERT_DELTA(factory.GetC2ForGeneration(0), 101.9786, 1e-6);
            TS_ASSERT_DELTA(factory.GetAlpha0ForGeneration(0), 0.882, 1e-6);
            TS_ASSERT_DELTA(factory.GetPleuralPressureForAirway(0, NULL), 0.0, 1e-6);

            //The three bifurcation mesh should map onto equivalent generation 0, generation 8 and generation 16 Lambert airways
            { //generation 0
                HiornsAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(0u));
                TS_ASSERT_DELTA(p_wall->mC1, 179380., 1e-6);
                delete p_wall;
            }

            { //generation 8
                HiornsAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(2u));
                TS_ASSERT_DELTA(p_wall->mphi1, 0., 1e-6);
                delete p_wall;
            }

            { //generation 16
                HiornsAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(4u));
                TS_ASSERT_DELTA(p_wall->mA, 0., 1e-6);
                delete p_wall;
            }

            for (unsigned element_index = 0; element_index < mesh.GetNumElements(); ++element_index)
            {
                HiornsAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(element_index));

                p_wall->SetTimestep(0.0);
                p_wall->SetAirwayPressure(0.0);
                p_wall->SetPleuralPressure(0.0);

                //double testVal = 0.00907487;
                double testVal = 0.0249;
                p_wall->SolveAndUpdateState(0.0, 0.0);
                TS_ASSERT_DELTA(p_wall->GetLumenRadius(), testVal, 1e-3);

                delete p_wall;
            }

            }


            //Repeat using Strahler order (equivalent for this mesh)
        {
            HiornsAirwayWallFactory factory(true);
            factory.SetMesh(&mesh);

            //The three bifurcation mesh should map onto equivalent generation 0, generation 8 and generation 16 Lambert airways
            { //generation 0
                HiornsAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(0u));
                TS_ASSERT_DELTA(p_wall->mmu, 64002., 1e-6);
                delete p_wall;
            }

            { //generation 2
                HiornsAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(2u));
                TS_ASSERT_DELTA(p_wall->mC1, 32.9000, 1e-6);
                delete p_wall;
            }

            { //generation 4
                HiornsAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(4u));
                TS_ASSERT_DELTA(p_wall->mC2, 0.0016, 1e-3);
                delete p_wall;
            }

            for (unsigned element_index = 0; element_index < mesh.GetNumElements(); ++element_index)
            {
                HiornsAirwayWall* p_wall = factory.CreateAirwayWallForElement(mesh.GetElement(element_index));

                p_wall->SetTimestep(0.0);
                p_wall->SetAirwayPressure(0.0);
                p_wall->SetPleuralPressure(0.0);

                p_wall->SolveAndUpdateState(0.0, 0.0);

                double testVal = 0.0249;
                //double testVal = 0.00907487;
                TS_ASSERT_DELTA(p_wall->GetLumenRadius(), testVal, 1e-3);

                delete p_wall;
            }
        }

    }

};
#endif /*_TESTAIRWAYWALLMODELS_HPP_*/

