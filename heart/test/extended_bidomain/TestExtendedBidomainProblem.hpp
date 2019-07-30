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


#ifndef _TESTEXTENDEDBIDOMAINPROBLEM_HPP_
#define _TESTEXTENDEDBIDOMAINPROBLEM_HPP_

#include "UblasIncludes.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "PetscSetupAndFinalize.hpp"

#include "HeartConfig.hpp"
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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        CorriasBuistICCModified *cell;
        cell = new CorriasBuistICCModified(mpSolver, mpZeroStimulus);

        double x = pNode->rGetLocation()[0];
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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
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

    void SetupParameters()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(5.0));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(1.0));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,10.0);

        //HeartConfig::Instance()->SetKSPSolver("gmres");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
    }

    /**
     * This test is aimed at comparing the extended bidomain implementation in Chaste with
     * the original Finite Difference code developed by Martin Buist.
     *
     * All the parameters are chosen to replicate the same conditions as in his code.
     */
    void TestExtendedProblemVsMartincCode()
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
         * As Martin's code is an FD code, results will never match exactly.
         * The comparison below is done against a 'valid' h5 file.
         *
         * The h5 file (1DValid.h5) is a Chaste (old phi_i formulation) file with is valid because, when extrapolating results from it, they look very similar
         * (except for a few points at the end of the upstroke) to the results taken
         * directly from Martin's code.
         * A plot of Chaste results versus Martin's result (at node 50) is stored
         * in the file 1DChasteVsMartin.eps for reference.
         *
         * A second plot comparing the old formulation (with phi_i) to the new formulation with V_m is contained in
         *.1DChasteNewFormulation.png
         *
         */
         TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/extendedbidomain", "1DValid", false,
                                                dir, filename, true,
                                                0.2));
         /*
          * Here we compare the new formulation (V_m1, V_m2, phi_e)
          *  with the previous formulation (phi_i1, phi_i2, phi_e) running with GMRES and an absolute KSP tolerance of 1e-8.
          */
         TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/extendedbidomain", "extended1d_previous_chaste_formulation_abs_tol_1e-8", false,
                                                dir, filename, true,
                                                1e-2));
    }

    // Test the functionality for outputting the values of requested cell state variables
    void TestExtendedBidomainProblemPrintsMultipleVariables()
    {
        // Get the singleton in a clean state
        HeartConfig::Instance()->Reset();

        // Set configuration file
        std::string dir = "ExtBidoMultiVars";
        std::string filename = "extended";
        HeartConfig::Instance()->SetOutputDirectory(dir);
        HeartConfig::Instance()->SetOutputFilenamePrefix(filename);

        HeartConfig::Instance()->SetSimulationDuration(0.1);
        // HeartConfig::Instance()->SetKSPSolver("gmres");
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
                AbstractCardiacCellInterface* p_cell = ext_problem.GetTissue()->GetCardiacCell(global_index);
                TS_ASSERT_DELTA(values.back(), p_cell->GetAnyVariable(output_variables[i],0), 1e-12);
            }

            //check the extra files for extra variables are there (the content is tested in the converter's tests)
            FileFinder file(dir + "/output/"+ filename +"_"+ output_variables[i] + ".dat", RelativeTo::ChasteTestOutput);
            TS_ASSERT(file.Exists());
        }
    }
};

#endif /*_TESTEXTENDEDBIDOMAINPROBLEM_HPP_*/
