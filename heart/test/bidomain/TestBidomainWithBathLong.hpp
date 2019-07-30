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


#ifndef TESTBIDOMAINWITHBATHLONG_HPP_
#define TESTBIDOMAINWITHBATHLONG_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "LuoRudy1991.hpp"
#include "BidomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "HeartRegionCodes.hpp"
#include "TrianglesMeshReader.hpp"
#include "SimpleStimulus.hpp"

template<unsigned DIM>
class BathCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    c_vector<double,DIM> mStimulatedPoint;

public:
    BathCellFactory(double stimulusMagnitude, c_vector<double,DIM> stimulatedPoint)
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new SimpleStimulus(stimulusMagnitude, 0.5)),
          mStimulatedPoint(stimulatedPoint)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<DIM>* pNode)
    {
        // paranoia - check this is really a tissue node
        assert(HeartRegionCode::IsRegionTissue( pNode->GetRegion() ));

        // stimulate centre node normally..
        bool is_centre;

        if (DIM==1)
        {
            is_centre = (fabs(pNode->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6);
        }
        else if (DIM==2)
        {
            is_centre = (    (fabs(pNode->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6)
                          && (fabs(pNode->GetPoint()[1]-mStimulatedPoint(1)) < 1e-6) );
        }
        else
        {
            is_centre = (    (fabs(pNode->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6)
                          && (fabs(pNode->GetPoint()[1]-mStimulatedPoint(1)) < 1e-6)
                          && (fabs(pNode->GetPoint()[2]-mStimulatedPoint(2)) < 1e-6) );
        }

        if (is_centre)
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
    }
};

class TestBidomainWithBathLong : public CxxTest::TestSuite
{
public:
    void Test3dBathIntracellularStimulation()
    {
        HeartConfig::Instance()->SetSimulationDuration(1);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_3d");

        c_vector<double,3> centre;
        centre(0) = 0.05;
        centre(1) = 0.05;
        centre(2) = 0.05;
        BathCellFactory<3> cell_factory(-2.5e7, centre); // stimulates x=0.05 node

        BidomainProblem<3> bidomain_problem( &cell_factory, true );

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructRegularSlabMesh(0.01, 0.1, 0.1, 0.1);

        // Set everything outside a central sphere (radius 0.4) to be bath
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            double y = mesh.GetElement(i)->CalculateCentroid()[1];
            double z = mesh.GetElement(i)->CalculateCentroid()[2];
            if (sqrt((x-0.05)*(x-0.05) + (y-0.05)*(y-0.05) + (z-0.05)*(z-0.05)) > 0.04)
            {
                mesh.GetElement(i)->SetAttribute(HeartRegionCode::GetValidBathId());
            }
        }

        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        // test V = 0 for all bath nodes
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (HeartRegionCode::IsRegionBath( mesh.GetNode(i)->GetRegion() )) // bath
            {
                TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
            }
        }

        // a hardcoded value
        TS_ASSERT_DELTA(sol_repl[2*404], 39.6833, 1e-3);
    }

// see #1061
    void Test2dBathExtracellularStimulusOneEdgeGroundedOnOppositeEdge()
    {
        HeartConfig::Instance()->SetSimulationDuration(40);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath2dExtraStimGrounded");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_2d");

        HeartConfig::Instance()->SetOdeTimeStep(0.005);  //ms

        // need to create a cell factory but don't want any intra stim, so magnitude
        // of stim is zero.
        c_vector<double,2> centre;
        centre(0) = 0.05;
        centre(1) = 0.05;
        BathCellFactory<2> cell_factory( 0.0, centre);

        BidomainProblem<2> bidomain_problem( &cell_factory, true );

        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Set everything outside a central circle (radius 0.4) to be bath
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            double y = mesh.GetElement(i)->CalculateCentroid()[1];
            if (sqrt((x-0.05)*(x-0.05) + (y-0.05)*(y-0.05)) > 0.02)
            {
                mesh.GetElement(i)->SetAttribute(HeartRegionCode::GetValidBathId());
            }
        }

        bidomain_problem.SetMesh(&mesh);

        //boundary flux for Phi_e
        //-4e3 is enough to trigger an action potential, -3e3 is below threshold, -5e3 crashes the cell model.
        double boundary_flux = -9e3;
        double duration = 2.5; //ms

        HeartConfig::Instance()->SetElectrodeParameters(true,0,boundary_flux, 0.0, duration);
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        bool ap_triggered = false;
        /*
         * We are checking the last time step. This test will only make sure that an upstroke is triggered.
         * We ran longer simulation for 350 ms and a nice AP was observed.
         */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            // test V = 0 for all bath nodes
            if (mesh.GetNode(i)->GetRegion()==1) // bath
            {
                TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
            }
            else if (sol_repl[2*i] > 0.0)
            {
                ap_triggered = true;
            }
        }
        TS_ASSERT(ap_triggered);
    }

// see #1061
    void Test3dBathExtracellularStimulusOneEdgeGroundedOnOppositeEdge()
    {
        HeartConfig::Instance()->SetSimulationDuration(6);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath3dExtraStimGrounded");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_3d");

        HeartConfig::Instance()->SetOdeTimeStep(0.005);  //ms

        // need to create a cell factory but don't want any intra stim, so magnitude
        // of stim is zero.
        c_vector<double,3> centre;
        centre(0) = 0.1;
        centre(1) = 0.1;
        centre(2) = 0.1;
        BathCellFactory<3> cell_factory( 0.0, centre);

        BidomainProblem<3> bidomain_problem( &cell_factory, true );

        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_1016_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Set everything outside a central sphere (radius 0.4) to be bath
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            double y = mesh.GetElement(i)->CalculateCentroid()[1];
            double z = mesh.GetElement(i)->CalculateCentroid()[2];
            if (sqrt((x-0.1)*(x-0.1) + (y-0.1)*(y-0.1) + (z-0.1)*(z-0.1)) > 0.03)
            {
                mesh.GetElement(i)->SetAttribute(1);
            }
        }

        bidomain_problem.SetMesh(&mesh);

        //boundary flux for Phi_e
        double boundary_flux = -4e3;
        double duration = 2.5; //ms

        HeartConfig::Instance()->SetElectrodeParameters(true,0,boundary_flux, 0.0, duration);
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        bool ap_triggered = false;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            // test V = 0 for all bath nodes
            if (HeartRegionCode::IsRegionBath( mesh.GetNode(i)->GetRegion() )) // bath
            {
                TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
            }
            else if (sol_repl[2*i] > 0.0)
            {
                ap_triggered = true;
            }
        }
        TS_ASSERT(ap_triggered);
    }
};


#endif /*TESTBIDOMAINWITHBATHLONG_HPP_*/
