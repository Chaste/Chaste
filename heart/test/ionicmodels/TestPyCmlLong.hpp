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

#ifndef TESTPYCMLLONG_HPP_
#define TESTPYCMLLONG_HPP_

#include "PyCmlLongHelperTestSuite.hpp"

#include <boost/foreach.hpp>
#include <vector>

#include "HeartConfig.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Test PyCml functionality by dynamically loading (and hence converting) a wide
 * range of cell models.
 *
 * May need a test-suite setup or similar to define model-specific parameters?
 * Should we pick up the list of models by reading the folder heart/test/data/cellml?
 */
class TestPyCmlLong : public PyCmlLongHelperTestSuite
{
public:
    void TestNormalCells()
    {
        std::cout << "Search for 'Failure', ': ***', 'Error', or 'Failed' to find problems." << std::endl;

        std::string dirname("TestPyCmlLongNormal");
        std::vector<std::string> args;
        args.push_back("--Wu");
        std::vector<std::string> models;
        AddAllModels(models);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.1, 1.0);
        RunTests(dirname, models, args);
    }

    void TestOptimisedCells()
    {
        std::string dirname("TestPyCmlLongOpt");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--opt");
        std::vector<std::string> models;
        AddAllModels(models);
        RunTests(dirname, models, args, true);
    }

    void TestCvodeCells()
    {
#ifdef CHASTE_CVODE
        std::string dirname("TestPyCmlLongCvodeNumericalJ");
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
        std::string dirname("TestPyCmlLongCvodeAnalyticJ");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--cvode");
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
        std::string dirname("TestPyCmlLongBE");
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
        std::string dirname("TestPyCmlLongRushLarsen");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--rush-larsen");
        std::vector<std::string> models;
        AddAllModels(models);
        // Three models require a slightly smaller timestep with RL than normal forward Euler:
        // Courtemanche 1998, Demir 1994, and Grandi 2010.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.001, 0.1, 1.0);
        RunTests(dirname, models, args, false, 0, false);
    }

    void TestRushLarsenOptCells()
    {
        std::string dirname("TestPyCmlLongRushLarsenOpt");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--rush-larsen");
        args.push_back("--opt");
        std::vector<std::string> models;
        AddAllModels(models);
        RunTests(dirname, models, args, true, -1000, false);
    }
};

#endif // TESTPYCMLLONG_HPP_
