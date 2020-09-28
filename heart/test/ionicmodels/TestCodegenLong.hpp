/*

Copyright (c) 2005-2020, University of Oxford.
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

#ifndef TESTCODEGENLONG_HPP_
#define TESTCODEGENLONG_HPP_

#include "CodegenLongHelperTestSuite.hpp"

#include <boost/foreach.hpp>
#include <vector>

#include "HeartConfig.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Test chaste_codegen functionality by dynamically loading (and hence converting) a wide
 * range of cell models.
 *
 * May need a test-suite setup or similar to define model-specific parameters?
 * Should we pick up the list of models by reading the folder heart/test/data/cellml?
 */
class TestCodegenLong : public CodegenLongHelperTestSuite
{
public:
    void TestNormalCells()
    {
        std::cout << "Search for 'Failure', ': ***', 'Error', or 'Failed' to find problems." << std::endl;

        std::string dirname("TestCodegenLongNormal");
        std::vector<std::string> args; 
        args.push_back("--Wu");
        std::vector<std::string> models;
        AddAllModels(models);

        std::vector<std::string> small_dt_models; // Models that need a very small dt
        small_dt_models.push_back("li_mouse_2010");
        BOOST_FOREACH (std::string small_dt_model, small_dt_models)
        {
            models.erase(std::find(models.begin(), models.end(), small_dt_model));
        }

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.1, 1.0);
        RunTests(dirname, models, args);

        // See Cooper Spiteri Mirams paper table 2
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001953125, 0.1, 1.0);
        RunTests(dirname + "-small-dt", small_dt_models, args);
    }

    void TestOptimisedCells()
    {
        std::string dirname("TestCodegenLongOpt");
        std::vector<std::string> models;
        std::vector<std::string> args;
        args.push_back("--opt");
        AddAllModels(models);

        std::vector<std::string> small_dt_models; // Models that need a very small dt
        small_dt_models.push_back("li_mouse_2010");
        BOOST_FOREACH (std::string small_dt_model, small_dt_models)
        {
            models.erase(std::find(models.begin(), models.end(), small_dt_model));
        }

        std::vector<std::string> different_lookup_table_models; // Models that need a different lookup table
        different_lookup_table_models.push_back("fink_noble_giles_model_2008");
        BOOST_FOREACH (std::string model, different_lookup_table_models)
        {
            models.erase(std::find(models.begin(), models.end(), model));
        }
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.1, 1.0);
        RunTests(dirname, models, args, true);

        // See Cooper Spiteri Mirams paper table 2
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001953125, 0.1, 1.0);
        RunTests(dirname + "-small-dt", small_dt_models, args, true);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.1, 1.0);
        RunTests(dirname + "-different_lookup_table", different_lookup_table_models,
                 {"--opt", "--lookup-table", "membrane_voltage", "-250.0005", "549.9999", "0.001"},
                 true);
    }

    void TestCvodeCells()
    {
#ifdef CHASTE_CVODE
        std::string dirname("TestCodegenLongCvodeNumericalJ");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--cvode");
        std::vector<std::string> models;
        AddAllModels(models);

        SetUseCvodeJacobian(false);
        RunTests(dirname, models, args);
        SetUseCvodeJacobian(true);
#endif
    }

    void TestAnalyticCvodeCells()
    {
#ifdef CHASTE_CVODE
        std::string dirname("TestCodegenLongCvodeAnalyticJ");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--cvode");
        args.push_back("--use-analytic-jacobian");
        std::vector<std::string> models;
        AddAllModels(models);

        // These have NaN in the jacobian due to massive exponentials
        std::vector<std::string> bad_models = boost::assign::list_of("aslanidi_model_2009")("hund_rudy_2004_a")("livshitz_rudy_2007");
        BOOST_FOREACH (std::string bad_model, bad_models)
        {
            models.erase(std::find(models.begin(), models.end(), bad_model));
        }

        RunTests(dirname, models, args);
#endif
    }

    void TestBackwardEulerCells()
    {
        std::string dirname("TestCodegenLongBE");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--backward-euler");

        std::vector<std::string> models;
        AddAllModels(models);

        std::vector<std::string> diff_models; // Models that need a smaller dt
        diff_models.push_back("iyer_model_2004");
        diff_models.push_back("iyer_model_2007");
        diff_models.push_back("jafri_rice_winslow_model_1998");
        diff_models.push_back("pandit_model_2001_epi");
        diff_models.push_back("priebe_beuckelmann_model_1998");
        diff_models.push_back("viswanathan_model_1999_epi");
        diff_models.push_back("winslow_model_1999");
        BOOST_FOREACH (std::string diff_model, diff_models)
        {
            models.erase(std::find(models.begin(), models.end(), diff_model));
        }

        // These have NaN in the jacobian due to massive exponentials
        std::vector<std::string> bad_models = boost::assign::list_of("hund_rudy_2004_a");
        BOOST_FOREACH (std::string bad_model, bad_models)
        {
            models.erase(std::find(models.begin(), models.end(), bad_model));
        }

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, models, args, true);

        dirname = dirname + "-difficult";
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.001, 0.1, 1.0);
        RunTests(dirname, diff_models, args, true);
    }

    void TestRushLarsenCells()
    {
        std::string dirname("TestCodegenLongRushLarsen");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--rush-larsen");
        std::vector<std::string> models;
        AddAllModels(models);

        std::vector<std::string> small_dt_models; // Models that need a very small dt
        small_dt_models.push_back("li_mouse_2010");
        BOOST_FOREACH (std::string small_dt_model, small_dt_models)
        {
            models.erase(std::find(models.begin(), models.end(), small_dt_model));
        }

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.001, 0.1, 1.0);
        RunTests(dirname, models, args, false, 0, false);

        // See Cooper Spiteri Mirams paper table 2
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001953125, 0.1, 1.0);
        RunTests(dirname + "-small-dt", small_dt_models, args, false, 0, false);
    }

    void TestRushLarsenOptCells()
    {
        std::string dirname("TestCodegenLongRushLarsenOpt");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--rush-larsen");
        args.push_back("--opt");
        std::vector<std::string> models;
        AddAllModels(models);

        std::vector<std::string> small_dt_models; // Models that need a very small dt
        small_dt_models.push_back("li_mouse_2010");
        small_dt_models.push_back("courtemanche_ramirez_nattel_model_1998");
        small_dt_models.push_back("demir_model_1994");
        small_dt_models.push_back("grandi2010ss");
        BOOST_FOREACH (std::string small_dt_model, small_dt_models)
        {
            models.erase(std::find(models.begin(), models.end(), small_dt_model));
        }

        std::vector<std::string> different_lookup_table_models; // Models for which we allow warnings
        different_lookup_table_models.push_back("DiFrancescoNoble1985");
        BOOST_FOREACH (std::string model, different_lookup_table_models)
        {
            models.erase(std::find(models.begin(), models.end(), model));
        }

        RunTests(dirname, models, args, true, -1000, false);

        RunTests(dirname + "-different_lookup_table", different_lookup_table_models,
                   {"--rush-larsen", "--opt", "--lookup-table", "membrane_voltage", "-250.0005", "549.9999", "0.001"},
                   true, -1000, true);

        // See Cooper Spiteri Mirams paper table 2
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001953125, 0.1, 1.0);
        RunTests(dirname + "-small-dt", small_dt_models, args, true, -1000, true);
    }

    void TestCvodeCellsOpt()
    {
#ifdef CHASTE_CVODE
        std::string dirname("TestCodegenLongCvodeNumericalJ-opt");
        std::vector<std::string> args;
        args.push_back("--cvode");
        args.push_back("--opt");

        std::vector<std::string> models;
        AddAllModels(models);

        std::vector<std::string> different_lookup_table_models; // Models we need to run it wirth finer tolerances
        different_lookup_table_models.push_back("ten_tusscher_model_2004_endo");
        BOOST_FOREACH (std::string model, different_lookup_table_models)
        {
            models.erase(std::find(models.begin(), models.end(), model));
        }

        SetUseCvodeJacobian(false);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, models, args);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001953125, 0.1, 1.0);
        RunTests(dirname, different_lookup_table_models,
                 {"--cvode", "--use-analytic-jacobian", "--opt", "--lookup-table", "membrane_voltage", "-250.0005", "549.9999", "0.001"},
                 true);

        SetUseCvodeJacobian(true);
#endif
    }

    void TestAnalyticCvodeCellsOpt()
    {
#ifdef CHASTE_CVODE
        std::string dirname("TestCodegenLongCvodeAnalyticJ-opt");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--cvode");
        args.push_back("--use-analytic-jacobian");
        args.push_back("--opt");
        std::vector<std::string> models;
        AddAllModels(models);

        // These have NaN in the jacobian due to massive exponentials
        std::vector<std::string> bad_models = boost::assign::list_of("aslanidi_model_2009")("hund_rudy_2004_a")("livshitz_rudy_2007");
        BOOST_FOREACH (std::string bad_model, bad_models)
        {
            models.erase(std::find(models.begin(), models.end(), bad_model));
        }

        std::vector<std::string> different_lookup_table_models; // Models we need to run it wirth finer tolerances
        different_lookup_table_models.push_back("ten_tusscher_model_2004_endo");
	different_lookup_table_models.push_back("noble_model_1991");
	different_lookup_table_models.push_back("luo_rudy_1994");
        BOOST_FOREACH (std::string model, different_lookup_table_models)
        {
            models.erase(std::find(models.begin(), models.end(), model));
        }

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, models, args);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001953125, 0.1, 1.0);
        RunTests(dirname, different_lookup_table_models,
                 {"--cvode", "--use-analytic-jacobian", "--opt", "--lookup-table", "membrane_voltage", "-250.0005", "549.9999", "0.001"},
                 true);

#endif
    }

    void TestBackwardEulerCellsOpt()
    {
        std::string dirname("TestCodegenLongBE-opt");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--backward-euler");
        args.push_back("--opt");

        std::vector<std::string> models;
        AddAllModels(models);

        std::vector<std::string> diff_models; // Models that need a smaller dt
        diff_models.push_back("iyer_model_2004");
        diff_models.push_back("iyer_model_2007");
        diff_models.push_back("jafri_rice_winslow_model_1998");
        diff_models.push_back("pandit_model_2001_epi");
        diff_models.push_back("priebe_beuckelmann_model_1998");
        diff_models.push_back("viswanathan_model_1999_epi");
        diff_models.push_back("winslow_model_1999");
        BOOST_FOREACH (std::string diff_model, diff_models)
        {
            models.erase(std::find(models.begin(), models.end(), diff_model));
        }

        // These have NaN in the jacobian due to massive exponentials
        std::vector<std::string> bad_models = boost::assign::list_of("hund_rudy_2004_a");
        BOOST_FOREACH (std::string bad_model, bad_models)
        {
            models.erase(std::find(models.begin(), models.end(), bad_model));
        }

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, models, args, true);

        dirname = dirname + "-difficult";
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.001, 0.1, 1.0);
        RunTests(dirname, diff_models, args, true);

    }
};

#endif // TESTCODEGENLONG_HPP_
