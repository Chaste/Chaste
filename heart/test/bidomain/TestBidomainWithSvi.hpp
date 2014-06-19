/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef TESTBIDOMAINWITHSVI_HPP_
#define TESTBIDOMAINWITHSVI_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "BidomainProblem.hpp"
//#include "BidomainWithBathProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PropagationPropertiesCalculator.hpp"

// stimulate a block of cells (an interval in 1d, a block in a corner in 2d)
template<unsigned DIM>
class BlockCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BlockCellFactory()
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new SimpleStimulus(-1000000.0, 0.5))
    {
        assert(DIM<3);
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<DIM>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y;
        if(DIM==2)
        {
            y = pNode->rGetLocation()[1];
        }

        if (    (DIM==1 && fabs(x)<0.02+1e-6)
             || (DIM==2 && fabs(x)<0.1+1e-6 && fabs(y)<0.1+1e-6) ) // 2D problem needs larger stimulus
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
    }
};

class BathCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStim;

public:
    BathCellFactory() : AbstractCardiacCellFactory<2>(),
    mpStim(new SimpleStimulus(-70000.0, 1.0, 0.0))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        AbstractCardiacCell* p_cell;
        ChastePoint<2> node_location = pNode->GetPoint();

        if ( node_location[0]<0.1 && node_location[1]<0.1)
        {
            p_cell = new CellLuoRudy1991FromCellML(mpSolver, mpStim);
        }
        else
        {
            p_cell = new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
        return p_cell;
    }
};

/* HOW_TO_TAG Cardiac/Solver
 * Use [ChasteGuides/StateVariableInterpolation state-variable interpolation] to improve accuracy
 */
class TestBidomainWithSvi : public CxxTest::TestSuite
{
public:
    void TestConductionVelocityConvergesFasterWithSvi1d() throw(Exception)
    {
        double h[3] = {0.001,0.01,0.02};
        std::vector<double> conduction_vel_nci(3);
        std::vector<double> conduction_vel_svi(3);

        ReplicatableVector final_voltage_ici;
        ReplicatableVector final_solution_svi;
        ReplicatableVector final_solution_svit;

        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        for(unsigned i=0; i<3; i++)
        {
            // ICI - ionic current interpolation - the default
            {
                TetrahedralMesh<1,1> mesh;
                mesh.ConstructRegularSlabMesh(h[i], 1.0);

                std::stringstream output_dir;
                output_dir << "BidomainIci_" << h[i];
                HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
                HeartConfig::Instance()->SetOutputFilenamePrefix("results");

                // need to have this for i=1,2 cases!!
                HeartConfig::Instance()->SetUseStateVariableInterpolation(false);

                BlockCellFactory<1> cell_factory;
                BidomainProblem<1> bidomain_problem( &cell_factory );
                bidomain_problem.SetMesh(&mesh);
                bidomain_problem.Initialise();

                bidomain_problem.Solve();

                final_voltage_ici.ReplicatePetscVector(bidomain_problem.GetSolution());
            }

            // SVI - state variable interpolation
            {
                DistributedTetrahedralMesh<1,1> mesh;
                mesh.ConstructRegularSlabMesh(h[i], 1.0);

                std::stringstream output_dir;
                output_dir << "BidomainSvi_" << h[i];
                HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
                HeartConfig::Instance()->SetOutputFilenamePrefix("results");

                HeartConfig::Instance()->SetUseStateVariableInterpolation();

                BlockCellFactory<1> cell_factory;
                BidomainProblem<1> bidomain_problem( &cell_factory );
                bidomain_problem.SetMesh(&mesh);
                bidomain_problem.Initialise();

                bidomain_problem.Solve();

                final_solution_svi.ReplicatePetscVector(bidomain_problem.GetSolution());
            }

            // SVIT - state variable interpolation on straight (not distributed) tetrahedral mesh
            {
                TetrahedralMesh<1,1> mesh;
                mesh.ConstructRegularSlabMesh(h[i], 1.0);

                std::stringstream output_dir;
                output_dir << "BidomainSviTet_" << h[i];
                HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
                HeartConfig::Instance()->SetOutputFilenamePrefix("results");

                HeartConfig::Instance()->SetUseStateVariableInterpolation();

                BlockCellFactory<1> cell_factory;
                BidomainProblem<1> bidomain_problem( &cell_factory );
                bidomain_problem.SetMesh(&mesh);
                bidomain_problem.Initialise();

                bidomain_problem.Solve();

                final_solution_svit.ReplicatePetscVector(bidomain_problem.GetSolution());
            }
            double voltage_at_0_03_finest_mesh;
            if(i==0) // finest mesh
            {
                for(unsigned j=0; j<final_voltage_ici.GetSize(); j++)
                {
                    // visually checked they agree at this mesh resolution, and chosen tolerance from results
                    TS_ASSERT_DELTA(final_voltage_ici[j], final_solution_svi[j], 0.35);
                    TS_ASSERT_DELTA(final_solution_svit[j], final_solution_svi[j], 1e-8);

                    if(j%2==0 /* just look at voltage */ && final_voltage_ici[j]>-80)
                    {
                        // shouldn't be exactly equal, as long as away from resting potential
                        TS_ASSERT_DIFFERS(final_voltage_ici[j], final_solution_svi[j]);
                    }
                }

                voltage_at_0_03_finest_mesh = final_voltage_ici[600];
                TS_ASSERT_DELTA(voltage_at_0_03_finest_mesh, -65.2218, 1e-2); //hardcoded value
            }
            else if(i==1)
            {
                double nci_voltage_at_0_03_middle_mesh = final_voltage_ici[60];
                double svi_voltage_at_0_03_middle_mesh = final_solution_svi[60];
                double svit_voltage_at_0_03_middle_mesh = final_solution_svit[60];
                // NCI conduction velocity > SVI conduction velocity
                // and both should be greater than CV on finesh mesh
                TS_ASSERT_DELTA(nci_voltage_at_0_03_middle_mesh, -44.3111, 1e-3);
                TS_ASSERT_DELTA(svi_voltage_at_0_03_middle_mesh, -60.7765, 1e-3);
                TS_ASSERT_DELTA(svit_voltage_at_0_03_middle_mesh, -60.7765, 1e-3);
            }
            else
            {
                double nci_voltage_at_0_03_coarse_mesh = final_voltage_ici[30];
                double svi_voltage_at_0_03_coarse_mesh = final_solution_svi[30];
                double svit_voltage_at_0_03_coarse_mesh = final_solution_svit[30];
                // NCI conduction velocity even greater than SVI conduction
                // velocity
                TS_ASSERT_DELTA(nci_voltage_at_0_03_coarse_mesh,  -6.5622, 1e-3);
                TS_ASSERT_DELTA(svi_voltage_at_0_03_coarse_mesh, -51.8848, 1e-3);
                TS_ASSERT_DELTA(svit_voltage_at_0_03_coarse_mesh, -51.8848, 1e-3);
            }
        }
    }

    void TestConductionVelocityInCrossFibreDirection2d() throw(Exception)
    {
        ReplicatableVector final_voltage_ici;
        ReplicatableVector final_solution_svi;

        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.01);

        // much lower conductivity in cross-fibre direction - NCI will struggle
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.17));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 0.7));

        // NCI - nodal current interpolation - the default
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.5, 0.3);

            HeartConfig::Instance()->SetOutputDirectory("BidomainNci2d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");
            
            HeartConfig::Instance()->SetUseStateVariableInterpolation(false);

            BlockCellFactory<2> cell_factory;
            BidomainProblem<2> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();
            bidomain_problem.Solve();

            final_voltage_ici.ReplicatePetscVector(bidomain_problem.GetSolution());
        }

        // SVI - state variable interpolation
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.5, 0.3);

            HeartConfig::Instance()->SetOutputDirectory("BidomainSvi2d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");
            
            HeartConfig::Instance()->SetUseStateVariableInterpolation(true);

            BlockCellFactory<2> cell_factory;
            BidomainProblem<2> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();

            bidomain_problem.Solve();

            final_solution_svi.ReplicatePetscVector(bidomain_problem.GetSolution());
        }

        // See comments in equivalent part of test/monodomain/TestMonodomainWithSvi.hpp
        // For bidomain with h=0.02, SVI CV is slower in both fibre and cross-fibre
        // directions

        // node 20 (for h=0.02) is on the x-axis (fibre direction)
        TS_ASSERT_DELTA(final_voltage_ici[20*2], -64.1105, 1e-3);
        TS_ASSERT_DELTA(final_solution_svi[20*2], -78.0936, 1e-3);
        // node 234 (for h=0.02) is on the y-axis (cross-fibre direction)
        TS_ASSERT_DELTA(final_voltage_ici[234*2], -57.7239, 1e-3);
        TS_ASSERT_DELTA(final_solution_svi[234*2], 38.9004, 1e-3);
    }

    void DontTestBidomainWithBathWithSvi() throw(Exception)
    {
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.04, 0.12, 0.12);

/*
        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter=mesh.GetElementIteratorBegin();
                        iter != mesh.GetElementIteratorEnd();
                        ++iter)
        {
            unsigned element_index = iter->GetIndex();
            std::cout << "Element index: " << element_index << std::endl;
            for ( unsigned i=0; i<3; ++i )
            {
                c_vector<double, 2> node_location = iter->GetNodeLocation(i);
                std::cout << "x: " << node_location[0]
                          << " y: " << node_location[1] << std::endl;
            }

            if ( element_index==8 || element_index==9 )
            {
                std::cout << "Setting as bath." << std::endl;
                iter->SetAttribute(HeartRegionCode::GetValidBathId());
            }
            else
            {
                std::cout << "Setting as tissue." << std::endl;
                iter->SetAttribute(HeartRegionCode::GetValidTissueId());
            }
        }
*/

        ChastePoint<2> bath_centre(0.12,0.12,0.12);
        ChastePoint<2> bath_radius(0.06,0.06,0.06);
        ChasteEllipsoid<2> bath_region( bath_centre, bath_radius );

        // Set some elements to bath, the rest tissue
        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter=mesh.GetElementIteratorBegin();
                iter != mesh.GetElementIteratorEnd();
                ++iter)
        {
            // Is a node within this region?
            if ( bath_region.DoesContain(iter->GetNode(0)->GetPoint())
                 //|| bath_region.DoesContain(iter->GetNode(1)->GetPoint())
                 //|| bath_region.DoesContain(iter->GetNode(2)->GetPoint())
               )
            {
                iter->SetAttribute(HeartRegionCode::GetValidBathId());

                std::cout << "Element index: " << iter->GetIndex() << std::endl;
                for ( unsigned i=0; i<3; ++i )
                {
                    c_vector<double, 2> node_location = iter->GetNodeLocation(i);
                    std::cout << "x: " << node_location[0]
                              << " y: " << node_location[1] << std::endl;
                }
            }
            else
            {
                iter->SetAttribute(HeartRegionCode::GetValidTissueId());
            }
        }

        HeartConfig::Instance()->SetSimulationDuration(10.0); // ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainWithBathSvi2d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.001, 0.025, 0.25);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.17));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 0.7));

        HeartConfig::Instance()->SetVisualizeWithMeshalyzer();
        HeartConfig::Instance()->SetVisualizeWithVtk();
        HeartConfig::Instance()->SetUseStateVariableInterpolation(true);

        BathCellFactory cell_factory;
        BidomainWithBathProblem<2> bidomain_problem( &cell_factory );
        bidomain_problem.SetMesh( &mesh );
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }
};

#endif /*TESTBIDOMAINWITHSVI_HPP_*/
