/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TESTCODEGENLONGANALYTICCVODEOPT_HPP_
#define TESTCODEGENLONGANALYTICCVODEOPT_HPP_

#include "CodegenLongHelperTestSuite.hpp"

#include <boost/foreach.hpp>
#include <vector>

#include "HeartConfig.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Test chaste_codegen functionality to generate Analytic Cvode Optimised cells,
 * by dynamically loading (and hence converting) a wide range of cell models.
 */
class TestCodegenLongAnalyticCvodeOpt : public CodegenLongHelperTestSuite
{
public:
    void TestAnalyticCvodeCellsOpt()
    {
#ifdef CHASTE_CVODE
        std::string dirname("TestCodegenLongCvodeAnalyticJ-opt");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--cvode");
        args.push_back("--use-analytic-jacobian");
        args.push_back("--opt");

        // These have NaN in the jacobian due to massive exponentials
        std::vector<std::string> bad_models = special_treatment_models(models, {"aslanidi_model_2009",
                                                                                  "hund_rudy_2004_a",
                                                                                  "livshitz_rudy_2007"});

        // Models that need a different lookup table
        std::vector<std::string> different_lookup_table_models = special_treatment_models(models, {"ten_tusscher_model_2004_endo",
                                                                                                     "noble_model_1991",
                                                                                                     "luo_rudy_1994"});


        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        // check test for negative concentrations works
        TS_ASSERT_THROWS_ANYTHING(RunTests(dirname, {"negative_concentration_paci_hyttinen_aaltosetala_severi_ventricularVersion"}, args));

        // initial value for membrane_L_type_calcium_current_f_gate in Shannon2004 is slightly above 1 (it is a probability)
        std::vector<std::string> bigger_tolerance_models = special_treatment_models(models, {"Shannon2004"});
        RunTests(dirname, models, args);

        SetUseCvOdeTolerances(1e-4);
        RunTests(dirname, bigger_tolerance_models, args);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001953125, 0.1, 1.0);
        RunTests(dirname, different_lookup_table_models, args);

#endif
    }
};

#endif // TESTCODEGENLONGANALYTICCVODEOPT_HPP_
