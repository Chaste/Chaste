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
#include "AbstractCvodeCell.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "LuoRudy1991Cvode.hpp"
#include "TenTusscher2006Epi.hpp"
#include "Mahajan2008.hpp"
#include "PlaneStimulusCellFactory.hpp"
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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<DIM>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y;
        if (DIM==2)
        {
            y = pNode->rGetLocation()[1];
        }

        if ((DIM==1 && fabs(x)<0.02+1e-6) || (DIM==2 && fabs(x)<0.02+1e-6 && fabs(y)<0.02+1e-6))
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
    }
};

// stimulate a block of cells (an interval in 1d, a block in a corner in 2d)
#ifdef CHASTE_CVODE
template<unsigned DIM>
class BlockCellFactoryCvode : public AbstractCardiacCellFactory<DIM>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BlockCellFactoryCvode()
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new SimpleStimulus(-1000000.0, 0.5))
    {
        assert(DIM<3);
    }

    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<DIM>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y;
        if (DIM==2)
        {
            y = pNode->rGetLocation()[1];
        }

        AbstractCvodeCell* p_cell;
        if ((DIM==1 && fabs(x)<0.02+1e-6) || (DIM==2 && fabs(x)<0.02+1e-6 && fabs(y)<0.02+1e-6))
        {
            p_cell =  new CellLuoRudy1991FromCellMLCvode(this->mpSolver, this->mpStimulus);
        }
        else
        {
            p_cell =  new CellLuoRudy1991FromCellMLCvode(this->mpSolver, this->mpZeroStimulus);
        }
        //p_cell->SetMinimalReset(true); // This is now done by AbstractCardiacCellFactory for any AbstractCvodeCell.
        p_cell->SetTolerances(1e-5,1e-7);
        return p_cell;
    }
};
#endif // CVODE

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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        if (x < 0.15)
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        }
        else if (x < 0.65)
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
    void TestConductionVelocityConvergesFasterWithSvi1d()
    {
        double h[3] = {0.001,0.01,0.02};
        unsigned probe_node_index[3] = {300, 30, 15};
        unsigned number_of_nodes[3] = {1001, 101, 51};
        std::vector<double> conduction_vel_ici(3);
        std::vector<double> conduction_vel_svi(3);

        ReplicatableVector final_voltage_ici;
        ReplicatableVector final_voltage_svi;

        //HeartConfig::Instance()->SetUseRelativeTolerance(1e-8);
        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        for (unsigned i=0; i<3; i++)
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
//                conduction_vel_ici[i] = ppc.CalculateConductionVelocity(node_at_0_04,node_at_0_40,0.36);
//                std::cout << "conduction_vel_ici = " << conduction_vel_ici[i] << "\n";
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

            if (i==0) // finest mesh
            {
                for (unsigned j=0; j<final_voltage_ici.GetSize(); j++)
                {
                    // visually checked they agree at this mesh resolution, and chosen tolerance from results
                    TS_ASSERT_DELTA(final_voltage_ici[j], final_voltage_svi[j], 0.3);

                    if (final_voltage_ici[j]>-80)
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
            else if (i==1)
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

    void TestConductionVelocityInCrossFibreDirection2d()
    {

        ReplicatableVector final_voltage_ici;
        ReplicatableVector final_voltage_svi;
        ReplicatableVector final_voltage_svit;

        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        //HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0005, 0.01, 0.01); //See comment below
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.01);

        // much lower conductivity in cross-fibre direction - ICI will struggle
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.17));

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.5, 0.3);

        // ICI - nodal current interpolation - the default
        {
            HeartConfig::Instance()->SetOutputDirectory("MonodomainIci2d");
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

#ifdef CHASTE_CVODE
        ReplicatableVector final_voltage_svi_cvode;
        // SVI - state variable interpolation with CVODE cells
        {
            HeartConfig::Instance()->SetOutputDirectory("MonodomainSvi2dCvode");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            HeartConfig::Instance()->SetUseStateVariableInterpolation();

            BlockCellFactoryCvode<2> cell_factory;
            MonodomainProblem<2> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();

            monodomain_problem.Solve();

            final_voltage_svi_cvode.ReplicatePetscVector(monodomain_problem.GetSolution());
        }
#endif //CHASTE_CVODE

        // SVIT - state variable interpolation on non-distributed tetrahedral mesh
        {

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

        double ici_20 = -17.1939; // These numbers are from a solve with ODE timestep of 0.0005,
        double svi_20 = -62.6336; // i.e. ten times smaller than that used in the test at present
        double ici_130 = 15.4282; // (for speed), hence large tolerances below.
        double svi_130 = 30.7389; // The tolerances are still nowhere near overlapping - i.e. ICI different to SVI

        // node 20 (for h=0.02) is on the x-axis (fibre direction), SVI CV is slower
        TS_ASSERT_DELTA(mesh.GetNode(20)->rGetLocation()[0],  0.4, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(20)->rGetLocation()[1],  0.0, 1e-9);

        TS_ASSERT_DELTA(final_voltage_ici[20], ici_20, 8.0);  // These tolerances show difference in parallel,
        TS_ASSERT_DELTA(final_voltage_svi[20], svi_20, 3.0);  // note that SVI is more stable in the presence of multicore...
#ifdef CHASTE_CVODE
        TS_ASSERT_DELTA(final_voltage_svi_cvode[20], svi_20, 3.0);
#endif //Cvode
        TS_ASSERT_DELTA(final_voltage_svit[20], svi_20, 3.0);

        // node 130 (for h=0.02) is on the y-axis (cross-fibre direction), ICI CV is slower
        TS_ASSERT_DELTA(mesh.GetNode(130)->rGetLocation()[0], 0.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(130)->rGetLocation()[1], 0.1, 1e-9);

        TS_ASSERT_DELTA(final_voltage_ici[130], ici_130, 1.0);
        TS_ASSERT_DELTA(final_voltage_svi[130], svi_130, 0.2);
#ifdef CHASTE_CVODE
        TS_ASSERT_DELTA(final_voltage_svi_cvode[130], svi_130, 0.3); // different CVODE versions = slightly different answer!
#endif //cvode
        TS_ASSERT_DELTA(final_voltage_svit[130], svi_130, 0.2);
    }

    void TestCoverage3d()
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

    void TestWithHeterogeneousCellModels()
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
        //TS_ASSERT_EQUALS(p_assembler->mElementsCanDoSvi.size(), 10u);

        // Therefore, we just test that calling Solve() runs (without the checking that
        // cell models are identical, this fails).
        monodomain_problem.Solve();
    }

    /*
     * This is the same as TestConductionVelocityConvergesFasterWithSvi1d with i=2, but solves in two parts.
     * If that test changes, check the hardcoded values here!
     */
    void TestArchiving()
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
    // Solve with a space step of 0.5mm.
    //
    // Note that this space step ought to be too big!
    //
    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainFailing()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1); //ms
//        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_20_elements");
//        HeartConfig::Instance()->SetSpaceDimension(1);
//        HeartConfig::Instance()->SetFibreLength(1.0, 1.0/20.0);
        HeartConfig::Instance()->SetSpaceDimension(3);
        HeartConfig::Instance()->SetSlabDimensions(0.2, 0.1, 0.1, 0.05);
        HeartConfig::Instance()->SetOutputDirectory("MonoFailing");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_SVI");
        HeartConfig::Instance()->SetUseStateVariableInterpolation(true);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory;
        MonodomainProblem<3> monodomain_problem(&cell_factory);

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
//        HeartConfig::Instance()->SetCapacitance(1.0);

        // the mesh is too coarse, and this simulation will result in cell gating
        // variables going out of range. An exception should be thrown in the
        // EvaluateYDerivatives() method of the cell model


        TS_ASSERT_THROWS_CONTAINS(monodomain_problem.Solve(),
#ifndef NDEBUG
            // VerifyStateVariables is only called when debug is on
            "State variable fast_sodium_current_m_gate__m has gone out of range. "
            "Check numerical parameters, for example time and space stepsizes");
#else
            // This test hits a later assert(!isnan) if we do a ndebug build.
            "Assertion tripped: !std::isnan(i_ionic)");
#endif // NDEBUG
    }
};

#endif /*TESTMONODOMAINWITHSVI_HPP_*/
