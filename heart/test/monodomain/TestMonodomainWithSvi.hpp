/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTMONODOMAINWITHSVI_HPP_
#define TESTMONODOMAINWITHSVI_HPP_


#include <cxxtest/TestSuite.h>
#include <vector>

#include "CardiacSimulationArchiver.hpp"
#include "MonodomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "TenTusscher2006Epi.hpp"
#include "Mahajan2008.hpp"
#include "PetscSetupAndFinalize.hpp"

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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNodeOrHaloNode(nodeIndex)->rGetLocation()[0];
        double y;
        if(DIM==2)
        {
            y = this->GetMesh()->GetNodeOrHaloNode(nodeIndex)->rGetLocation()[1];
        }

        if (    (DIM==1 && fabs(x)<0.02+1e-6)
             || (DIM==2 && fabs(x)<0.02+1e-6 && fabs(y)<0.02+1e-6) )
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
    }
};


// non-identical cell models
class HeterogeneousCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    HeterogeneousCellFactory()
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-1000000.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNodeOrHaloNode(nodeIndex)->rGetLocation()[0];
        if ( x<0.15 )
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        }
        else if (x < 0.65 )
        {
            return new CellMahajan2008FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
        else
        {
            return new CellTenTusscher2006EpiFromCellML(this->mpSolver, this->mpZeroStimulus);
        }
    }
};


class TestMonodomainWithSvi : public CxxTest::TestSuite
{
    void setUp()
    {
        HeartConfig::Instance()->Reset();
    }

public:
    void TestConductionVelocityConvergesFasterWithSvi1d() throw(Exception)
    {
        double h[3] = {0.001,0.01,0.02};
        unsigned probe_node_index[3] = {300, 30, 15};
        unsigned number_of_nodes[3] = {1001, 101, 51};
        std::vector<double> conduction_vel_nci(3);
        std::vector<double> conduction_vel_svi(3);

        ReplicatableVector final_voltage_ici;
        ReplicatableVector final_voltage_svi;

        //HeartConfig::Instance()->SetUseRelativeTolerance(1e-8);
        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        for(unsigned i=0; i<3; i++)
        {
            // ICI - ionic current interpolation - the default
            {
                DistributedTetrahedralMesh<1,1> mesh;
                mesh.ConstructRegularSlabMesh(h[i], 1.0);
                TS_ASSERT_EQUALS(mesh.GetNumNodes(), number_of_nodes[i]);

                //Double check (for later) that the indexing is as expected
                if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal( probe_node_index[i] ))
                {
                    TS_ASSERT_DELTA(mesh.GetNode( probe_node_index[i] )->rGetLocation()[0], 0.3, 1e-8);
                }
                std::stringstream output_dir;
                output_dir << "MonodomainIci_" << h[i];
                HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
                HeartConfig::Instance()->SetOutputFilenamePrefix("results");

                // need to have this for i=1,2 cases!!
                HeartConfig::Instance()->SetUseStateVariableInterpolation(false);

                BlockCellFactory<1> cell_factory;
                MonodomainProblem<1> monodomain_problem( &cell_factory );
                monodomain_problem.SetMesh(&mesh);
                monodomain_problem.Initialise();

                monodomain_problem.Solve();

                final_voltage_ici.ReplicatePetscVector(monodomain_problem.GetSolution());

//// see #1633
//// end time needs to be increased for these (say, to 7ms)
//                Hdf5DataReader simulation_data(OutputFileHandler::GetChasteTestOutputDirectory() + output_dir.str(),
//                                               "results", false);
//                PropagationPropertiesCalculator ppc(&simulation_data);
//                unsigned node_at_0_04 = (unsigned)round(0.04/h[i]);
//                unsigned node_at_0_40 = (unsigned)round(0.40/h[i]);
//                assert(fabs(mesh.GetNode(node_at_0_04)->rGetLocation()[0]-0.04)<1e-6);
//                assert(fabs(mesh.GetNode(node_at_0_40)->rGetLocation()[0]-0.40)<1e-6);
//                conduction_vel_nci[i] = ppc.CalculateConductionVelocity(node_at_0_04,node_at_0_40,0.36);
//                std::cout << "conduction_vel_nci = " << conduction_vel_nci[i] << "\n";
            }

            // SVI - state variable interpolation
            {
                DistributedTetrahedralMesh<1,1> mesh;
                mesh.ConstructRegularSlabMesh(h[i], 1.0);

                //Double check (for later) that the indexing is as expected
                if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal( probe_node_index[i] ))
                {
                    TS_ASSERT_DELTA(mesh.GetNode( probe_node_index[i] )->rGetLocation()[0], 0.3, 1e-8);
                }
                std::stringstream output_dir;
                output_dir << "MonodomainSvi_" << h[i];
                HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
                HeartConfig::Instance()->SetOutputFilenamePrefix("results");

                HeartConfig::Instance()->SetUseStateVariableInterpolation();

                BlockCellFactory<1> cell_factory;
                MonodomainProblem<1> monodomain_problem( &cell_factory );
                monodomain_problem.SetMesh(&mesh);
                monodomain_problem.Initialise();

                monodomain_problem.Solve();

                final_voltage_svi.ReplicatePetscVector(monodomain_problem.GetSolution());

//                Hdf5DataReader simulation_data(OutputFileHandler::GetChasteTestOutputDirectory() + output_dir.str(),
//                                               "results", false);
//                PropagationPropertiesCalculator ppc(&simulation_data);
//                unsigned node_at_0_04 = (unsigned)round(0.04/h[i]);
//                unsigned node_at_0_40 = (unsigned)round(0.40/h[i]);
//                assert(fabs(mesh.GetNode(node_at_0_04)->rGetLocation()[0]-0.04)<1e-6);
//                assert(fabs(mesh.GetNode(node_at_0_40)->rGetLocation()[0]-0.40)<1e-6);
//                conduction_vel_svi[i] = ppc.CalculateConductionVelocity(node_at_0_04,node_at_0_40,0.36);
//                std::cout << "conduction_vel_svi = " << conduction_vel_svi[i] << "\n";
            }

            if(i==0) // finest mesh
            {
                for(unsigned j=0; j<final_voltage_ici.GetSize(); j++)
                {
                    // visually checked they agree at this mesh resolution, and chosen tolerance from results
                    TS_ASSERT_DELTA(final_voltage_ici[j], final_voltage_svi[j], 0.3);

                    if(final_voltage_ici[j]>-80)
                    {
                        // shouldn't be exactly equal, as long as away from resting potential
                        TS_ASSERT_DIFFERS(final_voltage_ici[j], final_voltage_svi[j]);
                    }
                }

                double ici_voltage_at_0_03_finest_mesh = final_voltage_ici[ probe_node_index[i] ];
                double svi_voltage_at_0_03_finest_mesh = final_voltage_svi[ probe_node_index[i] ];
                TS_ASSERT_DELTA(svi_voltage_at_0_03_finest_mesh, 11.0067, 2e-4); //hardcoded value from fine svi
                TS_ASSERT_DELTA(ici_voltage_at_0_03_finest_mesh, 11.0067, 1.2e-1); //hardcoded value from fine svi
            }
            else if(i==1)
            {
                double ici_voltage_at_0_03_middle_mesh = final_voltage_ici[ probe_node_index[i] ];
                double svi_voltage_at_0_03_middle_mesh = final_voltage_svi[ probe_node_index[i] ];
                // ICI conduction velocity > SVI conduction velocity
                // and both should be greater than CV on finesh mesh
                TS_ASSERT_DELTA(ici_voltage_at_0_03_middle_mesh, 19.8924, 1e-3);
                TS_ASSERT_DELTA(svi_voltage_at_0_03_middle_mesh, 14.9579, 1e-3);
            }
            else
            {
                double ici_voltage_at_0_03_coarse_mesh = final_voltage_ici[ probe_node_index[i] ];
                double svi_voltage_at_0_03_coarse_mesh = final_voltage_svi[ probe_node_index[i] ];
                // ICI conduction velocity even greater than SVI conduction
                // velocity
                TS_ASSERT_DELTA(ici_voltage_at_0_03_coarse_mesh, 24.4938, 1e-3);
                TS_ASSERT_DELTA(svi_voltage_at_0_03_coarse_mesh, 17.3131, 1e-3);
            }
        }
    }

    void TestConductionVelocityInCrossFibreDirection2d() throw(Exception)
    {

        ReplicatableVector final_voltage_ici;
        ReplicatableVector final_voltage_svi;
        ReplicatableVector final_voltage_svit;

        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.01);

        // much lower conductivity in cross-fibre direction - ICI will struggle
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.17));

        // ICI - nodal current interpolation - the default
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.5, 0.3);

            HeartConfig::Instance()->SetOutputDirectory("MonodomainNci2d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            HeartConfig::Instance()->SetUseStateVariableInterpolation(false);

            BlockCellFactory<2> cell_factory;
            MonodomainProblem<2> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();
            monodomain_problem.Solve();

            final_voltage_ici.ReplicatePetscVector(monodomain_problem.GetSolution());
        }

        // SVI - state variable interpolation
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.5, 0.3);

            HeartConfig::Instance()->SetOutputDirectory("MonodomainSvi2d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            HeartConfig::Instance()->SetUseStateVariableInterpolation();

            BlockCellFactory<2> cell_factory;
            MonodomainProblem<2> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();

            monodomain_problem.Solve();

            final_voltage_svi.ReplicatePetscVector(monodomain_problem.GetSolution());
        }

        // SVIT - state variable interpolation on non-distributed tetrahedral mesh
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.5, 0.3);

            HeartConfig::Instance()->SetOutputDirectory("MonodomainSviTet2d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            HeartConfig::Instance()->SetUseStateVariableInterpolation();

            BlockCellFactory<2> cell_factory;
            MonodomainProblem<2> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();

            monodomain_problem.Solve();

            final_voltage_svit.ReplicatePetscVector(monodomain_problem.GetSolution());
        }

        // Visualised results with h=0.02 and h=0.01 - results looks sensible according to
        // paper:
        // 1. SVI h=0.01 and h=0.02 match more closely than ICI results - ie SVI converges faster
        // 2. CV in fibre direction faster for ICI (both values of h)
        // 3. CV in cross fibre direction: (i) h=0.01, faster for ICI; h=0.02 slower for ICI.
        // (Matches results in paper)

        // node 20 (for h=0.02) is on the x-axis (fibre direction), SVI CV is slower
        ///\todo #1462 Check that these nodes are where expected
        TS_ASSERT_DELTA(final_voltage_ici[20], -9.2270,  8e-3);  //These tolerances show difference in parallel - note that SVI is more stable in the presence of multicore...
        TS_ASSERT_DELTA(final_voltage_svi[20], -60.8510, 4e-3);
        TS_ASSERT_DELTA(final_voltage_svit[20], -60.8510, 4e-3);
        // node 130 (for h=0.02) is on the y-axis (cross-fibre direction), ICI CV is slower
        TS_ASSERT_DELTA(final_voltage_ici[130], 14.7918, 1e-3);
        TS_ASSERT_DELTA(final_voltage_svi[130], 30.5281, 1e-3);
        TS_ASSERT_DELTA(final_voltage_svit[130], 30.5281, 1e-3);
    }

    void TestCoverage3d() throw(Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(0.1); //ms
        HeartConfig::Instance()->SetUseStateVariableInterpolation(true);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.1);

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructRegularSlabMesh(0.02, 0.02, 0.02, 0.02);

        ZeroStimulusCellFactory<CellLuoRudy1991FromCellML,3> cell_factory;
        MonodomainProblem<3> monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh(&mesh);
        monodomain_problem.Initialise();
        monodomain_problem.Solve();
    }

    void TestWithHeterogeneousCellModels() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
        HeartConfig::Instance()->SetUseStateVariableInterpolation(true);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 1.0);

        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(0.01, 1.0);

        HeterogeneousCellFactory cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh(&mesh);
        monodomain_problem.Initialise();

        //// It's really very difficult to get hold of the correction assembler from here.
        //// The following, which requires 2 classes to be friends of this test compiles
        //// but fails asserts in the first line as bcc is not set up.
        //AbstractDynamicLinearPdeSolver<1,1,1>* p_solver = monodomain_problem.CreateSolver();
        //MonodomainSolver<1,1>* p_mono_solver = dynamic_cast<MonodomainSolver<1,1>*>(p_solver);
        //MonodomainCorrectionTermAssembler<1,1>* p_assembler = p_mono_solver->mpMonodomainCorrectionTermAssembler;
        //TS_ASSERT_EQUALS(p_assembler->mElementsHasIdenticalCellModels.size(), 10u);

        // Therefore, we just test that calling Solve() runs (without the checking that
        // cell models are identical, this fails).
        monodomain_problem.Solve();
    }

    /*
     * This is the same as TestConductionVelocityConvergesFasterWithSvi1d with i=2, but solves in two parts.
     * If that test changes, check the hardcoded values here!
     */
    void TestArchiving() throw (Exception)
    {
        std::string archive_dir = "monodomain_svi_archive";
        std::string archive_file = "monodomain_svi.arch";
        std::string output_dir = "monodomain_svi_output";

        { // Save
            HeartConfig::Instance()->SetSimulationDuration(0.1); //ms
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);
            HeartConfig::Instance()->SetOutputDirectory(output_dir);
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");
            HeartConfig::Instance()->SetUseStateVariableInterpolation();

            DistributedTetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(0.02, 1.0);

            BlockCellFactory<1> cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();
            monodomain_problem.Solve();

            CardiacSimulationArchiver<MonodomainProblem<1> >::Save(monodomain_problem, archive_dir);
        }

        HeartConfig::Instance()->Reset();

        { // Load
            HeartConfig::Instance()->SetUseStateVariableInterpolation(false); // Just in case...
            MonodomainProblem<1> *p_monodomain_problem = CardiacSimulationArchiver<MonodomainProblem<1> >::Load(archive_dir);

            HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
            p_monodomain_problem->Solve();

            ReplicatableVector final_voltage;
            final_voltage.ReplicatePetscVector(p_monodomain_problem->GetSolution());
            double probe_voltage = final_voltage[15];
            TS_ASSERT_DELTA(probe_voltage, 17.3131, 1e-3);

            delete p_monodomain_problem;
        }
    }
};

#endif /*TESTMONODOMAINWITHSVI_HPP_*/
