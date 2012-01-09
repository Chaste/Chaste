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


#ifndef _TESTEXTENDEDBIDOMAINPROBLEM_HPP_
#define _TESTEXTENDEDBIDOMAINPROBLEM_HPP_

#include "UblasIncludes.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "PetscSetupAndFinalize.hpp"

#include "HeartConfig.hpp"
#include "ArchiveOpener.hpp"
#include "SimpleStimulus.hpp"
#include "CorriasBuistSMCModified.hpp"
#include "CorriasBuistICCModified.hpp"
#include "CompareHdf5ResultsFiles.hpp"
#include "ExtendedBidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "FaberRudy2000.hpp"
#include "OutputFileHandler.hpp"
#include "Hdf5DataReader.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>



class ICC_Cell_factory : public AbstractCardiacCellFactory<1>
{
public:
    ICC_Cell_factory() : AbstractCardiacCellFactory<1>()
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        CorriasBuistICCModified *cell;
        cell = new CorriasBuistICCModified(mpSolver, mpZeroStimulus);

        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        double IP3_initial = 0.00067;
        double IP3_final = 0.00065;
        double cable_length = 10.0;

        //calculate the concentration gradient...
        double IP3_conc = IP3_initial + x*(IP3_final - IP3_initial)/cable_length;
        //..and set it
        cell->SetIP3Concentration(IP3_conc);
        cell->SetFractionOfVDDRInPU(0.04);

        return cell;
    }
};

class SMC_Cell_factory : public AbstractCardiacCellFactory<1>
{
public:
    SMC_Cell_factory() : AbstractCardiacCellFactory<1>()
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        CorriasBuistSMCModified *cell;
        cell = new CorriasBuistSMCModified(mpSolver, mpZeroStimulus);

        cell->SetFakeIccStimulusPresent(false);//it will get it from the real ICC, via gap junction
        return cell;
    }
};

class TestExtendedBidomainProblem: public CxxTest::TestSuite
{

public:

    void SetupParameters() throw (Exception)
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(5.0));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(1.0));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,10.0);

        HeartConfig::Instance()->SetKSPSolver("gmres");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-4);
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
    }

    /**
     * This test is aimed at comparing the extended bidomain implementation in Chaste with
     * the original Finite Difference code developed by Martin Buist.
     *
     * All the parameters are chosen to replicate the same conditions as in his code.
     */
    void TestExtendedProblemVsMartincCode() throw (Exception)
    {
        SetupParameters();

        TetrahedralMesh<1,1> mesh;
        unsigned number_of_elements = 100;//this is nGrid in Martin's code
        double length = 10.0;//100mm as in Martin's code
        mesh.ConstructRegularSlabMesh(length/number_of_elements, length);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), number_of_elements + 1);

        double Am_icc = 1000.0;
        double Am_smc = 1000.0;
        double Am_gap = 1.0;
        double Cm_icc = 1.0;
        double Cm_smc = 1.0;
        double G_gap = 20.0;//mS/cm^2

        HeartConfig::Instance()->SetSimulationDuration(1000.0);  //ms.

        ICC_Cell_factory icc_factory;
        SMC_Cell_factory smc_factory;

        std::string dir = "ICCandSMC";
        std::string filename = "extended1d";
        HeartConfig::Instance()->SetOutputDirectory(dir);
        HeartConfig::Instance()->SetOutputFilenamePrefix(filename);

        ExtendedBidomainProblem<1> extended_problem( &icc_factory , &smc_factory);
        extended_problem.SetMesh(&mesh);

        extended_problem.SetExtendedBidomainParameters(Am_icc,Am_smc, Am_gap, Cm_icc, Cm_smc, G_gap);
        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(1.0));

        std::vector<unsigned> outputnodes;
        outputnodes.push_back(50u);
        HeartConfig::Instance()->SetRequestedNodalTimeTraces(outputnodes);

        extended_problem.Initialise();
        extended_problem.Solve();
        HeartEventHandler::Headings();
        HeartEventHandler::Report();

        /**
         * Compare with valid data.
         * As Martin's code is an FD code, results wwill never match exactly.
         * The comparison below is done against a 'valid' h5 file.
         *
         * The h5 file is valid because, when extrapolating results from it, they look very similar
         * (excpet for a few points at the end of the upstroke) to the results taken
         * directly from Martin's code. A plot of Chaste results Vs Martin's result (at node 50) is stored
         * in the file 1DChasteVsMartin.eps for reference.
         *
         */
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/extendedbidomain", "1DValid", false,
                                                 dir, filename, true,
                                                 1.5e-3));
    }

    // Test the functionality for outputing the values of requested cell state variables
    void TestExtendedBidomainProblemPrintsMultipleVariables() throw (Exception)
    {
        // Get the singleton in a clean state
        HeartConfig::Instance()->Reset();

        // Set configuration file
        std::string dir = "ExtBidoMultiVars";
        std::string filename = "extended";
        HeartConfig::Instance()->SetOutputDirectory(dir);
        HeartConfig::Instance()->SetOutputFilenamePrefix(filename);

        HeartConfig::Instance()->SetSimulationDuration(0.1);
        HeartConfig::Instance()->SetKSPSolver("gmres");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-4);
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");

        /** Check that also the converters handle multiple variables**/
        HeartConfig::Instance()->SetVisualizeWithCmgui(true);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);

        TetrahedralMesh<1,1> mesh;
        unsigned number_of_elements = 100;
        double length = 10.0;
        mesh.ConstructRegularSlabMesh(length/number_of_elements, length);

        // Override the variables we are interested in writing.
        std::vector<std::string> output_variables;
        output_variables.push_back("calcium_dynamics__Ca_NSR");
        output_variables.push_back("ionic_concentrations__Nai");
        output_variables.push_back("fast_sodium_current_j_gate__j");
        output_variables.push_back("ionic_concentrations__Ki");

        HeartConfig::Instance()->SetOutputVariables( output_variables );

        // Set up problem
        PlaneStimulusCellFactory<CellFaberRudy2000FromCellML, 1> cell_factory_1(-60, 0.5);
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory_2(0.0,0.5);
        ExtendedBidomainProblem<1> ext_problem( &cell_factory_1, &cell_factory_2 );
        ext_problem.SetMesh(&mesh);

        // Solve
        ext_problem.Initialise();
        ext_problem.Solve();

        // Get a reference to a reader object for the simulation results
        Hdf5DataReader data_reader1 = ext_problem.GetDataReader();
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();

        // Check there is information about 11 timesteps (0, 0.01, 0.02, ...)
        unsigned num_steps = 11u;
        TS_ASSERT_EQUALS( times.size(), num_steps);
        TS_ASSERT_DELTA( times[0], 0.0, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.01, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.02, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.03, 1e-12);

        // There should be 11 values per variable and node.
        std::vector<double> node_5_v = data_reader1.GetVariableOverTime("V", 5);
        TS_ASSERT_EQUALS( node_5_v.size(), num_steps);

        std::vector<double> node_5_v_2 = data_reader1.GetVariableOverTime("V_2", 5);
        TS_ASSERT_EQUALS( node_5_v_2.size(), num_steps);

        std::vector<double> node_5_phi = data_reader1.GetVariableOverTime("Phi_e", 5);
        TS_ASSERT_EQUALS( node_5_phi.size(), num_steps);

        for (unsigned i=0; i<output_variables.size(); i++)
        {
            unsigned global_index = 2+i*2;
            std::vector<double> values = data_reader1.GetVariableOverTime(output_variables[i], global_index);
            TS_ASSERT_EQUALS( values.size(), num_steps);

            // Check the last values match the cells' state
            if (ext_problem.rGetMesh().GetDistributedVectorFactory()->IsGlobalIndexLocal(global_index))
            {
                AbstractCardiacCell* p_cell = ext_problem.GetTissue()->GetCardiacCell(global_index);
                TS_ASSERT_DELTA(values.back(), p_cell->GetAnyVariable(p_cell->GetAnyVariableIndex(output_variables[i])), 1e-12);
            }

            //check the extra files for extra variables are there (the content is tested in the converter's tests)
            std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
            std::string command = "ls " + test_output_directory + dir
                       + "/output/"+ filename +"_"+ output_variables[i] + ".dat ";
             TS_ASSERT_EQUALS(system(command.c_str()), 0);
        }
    }
};

#endif /*_TESTEXTENDEDBIDOMAINPROBLEM_HPP_*/
