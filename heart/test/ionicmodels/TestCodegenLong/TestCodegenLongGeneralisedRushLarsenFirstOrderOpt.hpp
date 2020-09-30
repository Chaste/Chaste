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

#ifndef TESTCODEGENLONGGENERALISEDRUSHLARSENFIRSTORDEROPT_HPP_
#define TESTCODEGENLONGGENERALISEDRUSHLARSENFIRSTORDEROPT_HPP_

#include <cxxtest/TestSuite.h>
#include "CodegenLongHelperTestSuite.hpp"

#include <boost/foreach.hpp>
#include <vector>

#include "HeartConfig.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Test chaste_codegen functionality to generate Generalised RushLarsen First Order Optimised cells,
 * by dynamically loading (and hence converting) a wide range of cell models.
 */
class TestCodegenLongGeneralisedRushLarsenFirstOrderOpt : public CodegenLongHelperTestSuite
{
public:
    void TestGeneralizedRushLarsenFirstOrderOpt()
    {
        std::string dirname("TestCodegenLongGeneralizedRushLarsen1Opt");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--grl1");
        args.push_back("--opt");
        std::vector<std::string> models;
        AddAllModels(models);
        models.erase(std::find(models.begin(), models.end(), "iyer_model_2004"));

        std::vector<std::string> smaller_timestep_models;
        smaller_timestep_models.push_back("Shannon2004");
        smaller_timestep_models.push_back("bondarenko_model_2004_apex");
        smaller_timestep_models.push_back("courtemanche_ramirez_nattel_model_1998");
        smaller_timestep_models.push_back("demir_model_1994");
        smaller_timestep_models.push_back("dokos_model_1996");
        smaller_timestep_models.push_back("iyer_model_2007");
        smaller_timestep_models.push_back("jafri_rice_winslow_model_1998");
        smaller_timestep_models.push_back("livshitz_rudy_2007");
        smaller_timestep_models.push_back("stewart_zhang_model_2008_ss");
        smaller_timestep_models.push_back("winslow_model_1999");
        smaller_timestep_models.push_back("FaberRudy2000");
        smaller_timestep_models.push_back("li_mouse_2010");
        smaller_timestep_models.push_back("noble_model_1998");
        BOOST_FOREACH (std::string model, smaller_timestep_models)
        {
            models.erase(std::find(models.begin(), models.end(), model));
        }



        std::vector<std::string> different_lookup_table_models;
        different_lookup_table_models.push_back("fink_noble_giles_model_2008");
        different_lookup_table_models.push_back("ten_tusscher_model_2006_epi");
        BOOST_FOREACH (std::string diff_model, different_lookup_table_models)
        {
            models.erase(std::find(models.begin(), models.end(), diff_model));
        }

        std::vector<std::string> different_lookup_table_models2;
        different_lookup_table_models2.push_back("decker_2009");
        different_lookup_table_models2.push_back("ten_tusscher_model_2004_epi");
        BOOST_FOREACH (std::string diff_model, different_lookup_table_models2)
        {
            models.erase(std::find(models.begin(), models.end(), diff_model));
        }

        std::vector<std::string> different_lookup_table_models3;
        different_lookup_table_models3.push_back("grandi2010ss");
        BOOST_FOREACH (std::string diff_model, different_lookup_table_models3)
        {
            models.erase(std::find(models.begin(), models.end(), diff_model));
        }


        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, models, args, true, -1000, false);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001, 0.1, 1.0);
        RunTests(dirname, smaller_timestep_models, args, true, -1000, false);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, different_lookup_table_models,
                 {"--opt", "--grl1", "--lookup-table", "membrane_voltage", "-250.0005", "549.9999", "0.001"},
                 true, -1000, true);

        RunTests(dirname, different_lookup_table_models2,
                 {"--opt", "--grl1", "--lookup-table", "membrane_voltage", "-250.0003", "549.9999", "0.001"},
                 true, -1000, true);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001, 0.1, 1.0);
        RunTests(dirname, different_lookup_table_models3,
                 {"--opt", "--grl1", "--lookup-table", "membrane_voltage", "-250.0005", "549.9999", "0.001"},
                 true, -1000, true);
    }
};

#endif // TESTCODEGENLONGGENERALISEDRUSHLARSENFIRSTORDEROPT_HPP_
