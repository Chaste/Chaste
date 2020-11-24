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

#ifndef TESTCODEGENLONGGENERALISEDRUSHLARSENSECONDORDEROPT_HPP_
#define TESTCODEGENLONGGENERALISEDRUSHLARSENSECONDORDEROPT_HPP_

#include "CodegenLongHelperTestSuite.hpp"

#include <boost/foreach.hpp>
#include <vector>

#include "HeartConfig.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Test chaste_codegen functionality to generate Generalised RushLarsen Second Order Optimised cells,
 * by dynamically loading (and hence converting) a wide range of cell models.
 */
class TestCodegenLongGeneralisedRushLarsenSecondOrderOpt : public CodegenLongHelperTestSuite
{

public:
    void TestGeneralizedRushLarsenSecondOrderOpt()
    {
        std::string dirname("TestCodegenLongGeneralizedRushLarsen2Opt");
        std::vector<std::string> args;
        args.push_back("--grl2");
        args.push_back("--opt");

        models.erase(std::find(models.begin(), models.end(), "nygren_atrial_model_1998"));

        std::vector<std::string> smaller_timestep_models = spectail_streatment_models(models, {"bondarenko_model_2004_apex",
                                                                                               "courtemanche_ramirez_nattel_model_1998",
                                                                                               "demir_model_1994",
                                                                                               "dokos_model_1996",
                                                                                               "iyer_model_2004",
                                                                                               "iyer_model_2007",
                                                                                               "jafri_rice_winslow_model_1998",
                                                                                               "livshitz_rudy_2007",
                                                                                               "Shannon2004",
                                                                                               "stewart_zhang_model_2008_ss",
                                                                                               "winslow_model_1999",
                                                                                               "FaberRudy2000",
                                                                                               "li_mouse_2010"});

        std::vector<std::string> very_small_timestep_models = spectail_streatment_models(models, {"nygren_atrial_model_1998"});

        std::vector<std::string> different_lookup_table_models = spectail_streatment_models(models, {"fink_noble_giles_model_2008",
                                                                                                     "ten_tusscher_model_2006_epi"});

        std::vector<std::string> different_lookup_table_models2 = spectail_streatment_models(models, {"decker_2009",
                                                                                                      "ten_tusscher_model_2004_epi"});

        std::vector<std::string> different_lookup_table_models3 = spectail_streatment_models(models, {"grandi2010ss"});

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, models, args, true, -1000, false);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001, 0.1, 1.0);
        RunTests(dirname, smaller_timestep_models, args, true, -1000, false);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001953125, 0.1, 1.0);
        RunTests(dirname, very_small_timestep_models, args, true, -1000, false);


        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, different_lookup_table_models,
                 {"--opt", "--grl2", "--lookup-table", "membrane_voltage", "-250.0005", "549.9999", "0.001"},
                 true, -1000, true);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, different_lookup_table_models2,
                 {"--opt", "--grl2", "--lookup-table", "membrane_voltage", "-250.0003", "549.9999", "0.001"},
                 true, -1000, true);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001, 0.1, 1.0);
        RunTests(dirname, different_lookup_table_models3,
                 {"--opt", "--grl2", "--lookup-table", "membrane_voltage", "-250.0005", "549.9999", "0.001"},
                 true, -1000, true);

    }
};

#endif // TESTCODEGENLONGGENERALISEDRUSHLARSENSECONDORDEROPT_HPP_
