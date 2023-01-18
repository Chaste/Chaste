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

#ifndef TESTCODEGENLONGEASYMODELS_HPP_
#define TESTCODEGENLONGEASYMODELS_HPP_

#include "CodegenLongHelperTestSuite.hpp"

#include <boost/foreach.hpp>
#include <vector>

#include "HeartConfig.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Test chaste_codegen functionality to generate Normal cells,
 * by dynamically loading (and hence converting) a wide range of cell models.
 */
class TestCodegenLongEasyModels : public CodegenLongHelperTestSuite
{
public:
    void TestAnalyticCvodeCellsOpt()
    {
#ifdef CHASTE_CVODE
        std::string dirname("TestCodegenLongCvodeAnalyticJ-easy_models-opt");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--cvode");
        args.push_back("--use-analytic-jacobian");
        args.push_back("--opt");

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        TS_ASSERT_THROWS_ANYTHING(RunTests(dirname, {"negative_concentration_paci_hyttinen_aaltosetala_severi_ventricularVersion"}, args));

        RunTests(dirname, easy_models, args);
#endif
    }

    void TestBackwardEulerCellsOpt()
    {
        std::string dirname("TestCodegenLongBE-easy_models-opt");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--backward-euler");
        args.push_back("--opt");

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, easy_models, args, true);
    }
};

#endif // TESTCODEGENLONGEASYMODELS_HPP_
