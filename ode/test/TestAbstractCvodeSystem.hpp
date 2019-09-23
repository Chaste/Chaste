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

#ifndef _TESTABSTRACTCVODESYSTEM_HPP_
#define _TESTABSTRACTCVODESYSTEM_HPP_

#include <cmath>
#include <iostream>
#include "CheckpointArchiveTypes.hpp"

#include "Cvode1.hpp"
#include "CvodeFirstOrder.hpp"
#include "ParameterisedCvode.hpp"
#include "TwoDimCvodeSystem.hpp"

#include "OdeSolution.hpp"
#include "OutputFileHandler.hpp"
#include "VectorHelperFunctions.hpp"

#include "FakePetscSetup.hpp"

double tol = 1e-6;

class TestAbstractCvodeSystem : public CxxTest::TestSuite
{
public:
    void TestOdeSystemOne()
    {
#ifdef CHASTE_CVODE
        // Test Ode1 class
        Cvode1 ode1;
        // dy
        N_Vector y = ode1.GetInitialConditions();
        N_Vector dy = NULL;
        CreateVectorIfEmpty(dy, 1);

        ode1.EvaluateYDerivatives(1.0, y, dy);
        TS_ASSERT_DELTA(GetVectorComponent(dy, 0), 1.0, tol);

        DeleteVector(y);
        DeleteVector(dy);
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestExceptions()
    {
#ifdef CHASTE_CVODE
        Cvode1 ode;
        TS_ASSERT_EQUALS(ode.GetNumberOfStateVariables(), 1u);

        N_Vector v = NULL;
        CreateVectorIfEmpty(v, 2);

        TS_ASSERT_THROWS_THIS(ode.SetDefaultInitialConditions(v),
                              "The number of initial conditions must be that of the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.SetDefaultInitialCondition(2, -3.0),
                              "Index is greater than the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.SetStateVariables(v),
                              "The size of the passed in vector must be that of the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.ForceUseOfNumericalJacobian(false),
                              "Analytic Jacobian requested, but this ODE system doesn't have one. You can check this with HasAnalyticJacobian().");

        DeleteVector(v);
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestParameters()
    {
#ifdef CHASTE_CVODE
        ParameterisedCvode ode;
        boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();

        TS_ASSERT_EQUALS(ode.GetSystemName(), "ParameterisedCvode");
        TS_ASSERT_EQUALS(p_info->GetSystemName(), "ParameterisedCvode");

        TS_ASSERT_EQUALS(ode.GetParameter(0), 0.0);
        TS_ASSERT_EQUALS(ode.GetNumberOfParameters(), 1u);
        TS_ASSERT_EQUALS(ode.GetNumberOfStateVariables(), 1u);

        double a = 1.0;
        ode.SetParameter(0, a);
        TS_ASSERT_EQUALS(ode.GetParameter(0), a);
        ode.SetParameter("a", a);
        TS_ASSERT_EQUALS(ode.GetParameter(0), a);
        TS_ASSERT_EQUALS(ode.GetParameter("a"), a);

        TS_ASSERT_EQUALS(ode.HasParameter("a"), true);
        TS_ASSERT_EQUALS(ode.rGetParameterNames()[0], "a");
        TS_ASSERT_EQUALS(ode.rGetParameterUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(ode.GetParameterIndex("a"), 0u);
        TS_ASSERT_EQUALS(ode.GetParameterUnits(0u), "dimensionless");

        TS_ASSERT_EQUALS(p_info->rGetParameterNames()[0], "a");
        TS_ASSERT_EQUALS(p_info->rGetParameterUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetParameterIndex("a"), 0u);
        TS_ASSERT_EQUALS(p_info->GetParameterUnits(0u), "dimensionless");

        TS_ASSERT_EQUALS(ode.HasStateVariable("y"), true);
        TS_ASSERT_EQUALS(ode.rGetStateVariableNames()[0], "y");
        TS_ASSERT_EQUALS(ode.rGetStateVariableUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(ode.GetStateVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(ode.GetStateVariableUnits(0u), "dimensionless");

        TS_ASSERT_EQUALS(p_info->rGetStateVariableNames()[0], "y");
        TS_ASSERT_EQUALS(p_info->rGetStateVariableUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetStateVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(p_info->GetStateVariableUnits(0u), "dimensionless");

        double y = -1.0;
        ode.SetAnyVariable(0u, y);
        TS_ASSERT_EQUALS(ode.GetStateVariable("y"), y);
        TS_ASSERT_EQUALS(ode.GetStateVariable(0u), y);

        ode.SetStateVariable("y", y - 1);
        TS_ASSERT_EQUALS(ode.GetStateVariable("y"), y - 1);
        ode.SetAnyVariable("y", y);
        TS_ASSERT_EQUALS(ode.GetStateVariable("y"), y);

        TS_ASSERT_EQUALS(ode.HasAnalyticJacobian(), false);

        TS_ASSERT_EQUALS(ode.HasAnyVariable("y"), true);
        TS_ASSERT_EQUALS(ode.HasAnyVariable("a"), true);
        TS_ASSERT_EQUALS(ode.GetAnyVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(ode.GetAnyVariableIndex("a"), 1u);
        TS_ASSERT_EQUALS(ode.GetAnyVariable(0u), y);
        TS_ASSERT_EQUALS(ode.GetAnyVariable("y"), y);
        TS_ASSERT_EQUALS(ode.GetAnyVariable(1u), a);
        TS_ASSERT_EQUALS(ode.GetAnyVariable("a"), a);
        double new_a = 12345.6789;
        ode.SetAnyVariable(1u, new_a);
        TS_ASSERT_EQUALS(ode.GetAnyVariable(1u), new_a);

        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits(0u), "dimensionless");
        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits(1u), "dimensionless");
        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits("y"), "dimensionless");
        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits("a"), "dimensionless");

        TS_ASSERT_EQUALS(p_info->GetAnyVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(p_info->GetAnyVariableIndex("a"), 1u);
        TS_ASSERT_EQUALS(p_info->GetAnyVariableUnits(0u), "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetAnyVariableUnits(1u), "dimensionless");

        // Exceptions
        TS_ASSERT_THROWS_THIS(ode.SetParameter(1u, -a), "The index passed in must be less than the number of parameters.");

        TS_ASSERT_EQUALS(ode.HasParameter("b"), false);
        TS_ASSERT_THROWS_THIS(ode.GetParameterIndex("b"), "No parameter named 'b'.");
        TS_ASSERT_THROWS_THIS(ode.GetParameter("b"), "No parameter named 'b'.");
        TS_ASSERT_EQUALS(ode.HasStateVariable("b"), false);
        TS_ASSERT_THROWS_THIS(ode.GetStateVariableIndex("b"), "No state variable named 'b'.");
        TS_ASSERT_EQUALS(ode.HasAnyVariable("b"), false);
        TS_ASSERT_THROWS_THIS(ode.GetAnyVariableIndex("b"), "No state variable, parameter, or derived quantity named 'b'.");
        TS_ASSERT_THROWS_THIS(ode.SetAnyVariable(2u, 0.0), "Cannot set the value of a derived quantity, or invalid index.");

        TS_ASSERT_THROWS_THIS(ode.GetAnyVariable(3u), "Invalid index passed to GetAnyVariable.");
        TS_ASSERT_THROWS_THIS(ode.GetStateVariable(1u), "The index passed in must be less than the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.GetParameter(1u), "The index passed in must be less than the number of parameters.");

        TS_ASSERT_THROWS_THIS(ode.GetParameterUnits(1u), "The index passed in must be less than the number of parameters.");
        TS_ASSERT_THROWS_THIS(ode.GetStateVariableUnits(1u), "The index passed in must be less than the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.GetAnyVariableUnits(3u), "Invalid index passed to GetAnyVariableUnits.");

        TS_ASSERT_EQUALS(p_info->HasParameter("b"), false);
        TS_ASSERT_THROWS_THIS(p_info->GetParameterIndex("b"), "No parameter named 'b'.");
        TS_ASSERT_EQUALS(p_info->HasStateVariable("b"), false);
        TS_ASSERT_THROWS_THIS(p_info->GetStateVariableIndex("b"), "No state variable named 'b'.");
        TS_ASSERT_EQUALS(p_info->HasAnyVariable("b"), false);
        TS_ASSERT_THROWS_THIS(p_info->GetAnyVariableIndex("b"), "No state variable, parameter, or derived quantity named 'b'.");
        TS_ASSERT_THROWS_THIS(p_info->GetParameterUnits(1u), "The index passed in must be less than the number of parameters.");
        TS_ASSERT_THROWS_THIS(p_info->GetStateVariableUnits(1u), "The index passed in must be less than the number of state variables.");
        TS_ASSERT_THROWS_THIS(p_info->GetAnyVariableUnits(3u), "Invalid index passed to GetAnyVariableUnits.");
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestAttributes()
    {
#ifdef CHASTE_CVODE
        ParameterisedCvode ode;
        TS_ASSERT_EQUALS(ode.GetNumberOfAttributes(), 1u);
        TS_ASSERT(ode.HasAttribute("attr"));
        TS_ASSERT_DELTA(ode.GetAttribute("attr"), 1.1, 1e-12);
        TS_ASSERT(!ode.HasAttribute("missing"));
        TS_ASSERT_THROWS_THIS(ode.GetAttribute("missing"), "No attribute 'missing' found.");

        boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();
        TS_ASSERT_EQUALS(p_info->GetNumberOfAttributes(), 1u);
        TS_ASSERT(p_info->HasAttribute("attr"));
        TS_ASSERT_DELTA(p_info->GetAttribute("attr"), 1.1, 1e-12);
        TS_ASSERT(!p_info->HasAttribute("missing"));
        TS_ASSERT_THROWS_THIS(p_info->GetAttribute("missing"), "No attribute 'missing' found.");
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestDerivedQuantities()
    {
#ifdef CHASTE_CVODE
        ParameterisedCvode ode;
        boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();

        TS_ASSERT_EQUALS(ode.GetNumberOfDerivedQuantities(), 1u);
        TS_ASSERT_EQUALS(ode.HasDerivedQuantity("2a_plus_y"), true);
        TS_ASSERT_EQUALS(ode.HasDerivedQuantity("Not_there"), false);
        TS_ASSERT_EQUALS(ode.GetDerivedQuantityIndex("2a_plus_y"), 0u);
        TS_ASSERT_EQUALS(ode.GetDerivedQuantityUnits(0u), "dimensionless");
        TS_ASSERT_EQUALS(ode.rGetDerivedQuantityNames().size(), 1u);
        TS_ASSERT_EQUALS(ode.rGetDerivedQuantityUnits().size(), 1u);
        TS_ASSERT_EQUALS(ode.rGetDerivedQuantityNames()[0], "2a_plus_y");
        TS_ASSERT_EQUALS(ode.rGetDerivedQuantityUnits()[0], "dimensionless");

        TS_ASSERT_EQUALS(p_info->HasDerivedQuantity("2a_plus_y"), true);
        TS_ASSERT_EQUALS(p_info->HasDerivedQuantity("Not_there"), false);
        TS_ASSERT_EQUALS(p_info->GetDerivedQuantityIndex("2a_plus_y"), 0u);
        TS_ASSERT_EQUALS(p_info->GetDerivedQuantityUnits(0u), "dimensionless");
        TS_ASSERT_EQUALS(p_info->rGetDerivedQuantityNames().size(), 1u);
        TS_ASSERT_EQUALS(p_info->rGetDerivedQuantityUnits().size(), 1u);
        TS_ASSERT_EQUALS(p_info->rGetDerivedQuantityNames()[0], "2a_plus_y");
        TS_ASSERT_EQUALS(p_info->rGetDerivedQuantityUnits()[0], "dimensionless");

        TS_ASSERT_EQUALS(ode.HasAnyVariable("2a_plus_y"), true);
        TS_ASSERT_EQUALS(ode.HasAnyVariable("Not_there"), false);
        TS_ASSERT_EQUALS(ode.GetAnyVariableIndex("2a_plus_y"), 2u);
        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits(2u), "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetAnyVariableIndex("2a_plus_y"), 2u);
        TS_ASSERT_EQUALS(p_info->GetAnyVariableUnits(2u), "dimensionless");

        N_Vector derived = ode.ComputeDerivedQuantitiesFromCurrentState(0.0);
        double a = ode.GetParameter(0);
        TS_ASSERT_EQUALS(a, 0.0);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 0), 2 * a, 1e-4);
        TS_ASSERT_DELTA(ode.GetAnyVariable(2u, 0.0), 2 * a, 1e-4);
        TS_ASSERT_DELTA(ode.GetAnyVariable(2u, 0.0, &derived), 2 * a, 1e-4);
        a = 1.0;
        ode.SetParameter(0, a);
        DeleteVector(derived); // Have to wipe contents from memory before assigning new ones...

        N_Vector ics = ode.GetInitialConditions();
        derived = ode.ComputeDerivedQuantities(0.0, ics);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 0), 2 * a, 1e-4);
        double y = 10.0;
        ode.SetStateVariable(0, y);
        DeleteVector(derived); // Have to wipe contents from memory before assigning new ones...
        DeleteVector(ics);

        derived = ode.ComputeDerivedQuantitiesFromCurrentState(0.0);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 0), 2 * a + y, 1e-4);
        TS_ASSERT_DELTA(ode.GetAnyVariable(2u, 1.0 /* ignored for this ODE */), 2 * a + y, 1e-4);

        // Exceptions
        TS_ASSERT_THROWS_THIS(ode.GetDerivedQuantityIndex("Missing"), "No derived quantity named 'Missing'.");
        TS_ASSERT_THROWS_THIS(ode.GetDerivedQuantityUnits(1u), "The index passed in must be less than the number of derived quantities.");

        TwoDimCvodeSystem ode2;
        TS_ASSERT_THROWS_THIS(ode2.ComputeDerivedQuantitiesFromCurrentState(0.0),
                              "This ODE system does not define derived quantities.");
        N_Vector doesnt_matter;
        TS_ASSERT_THROWS_THIS(ode2.ComputeDerivedQuantities(0.0, doesnt_matter),
                              "This ODE system does not define derived quantities.");

        DeleteVector(derived);
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    // This test is mainly for coverage purposes.
    void TestDumpState()
    {
#ifdef CHASTE_CVODE
        // Create a two variable system
        TwoDimCvodeSystem ode_system;

        // Dump the state variables
        std::string state = ode_system.DumpState("This is a test.");
        TS_ASSERT_EQUALS(state, "This is a test.\nState:\n\tVariable_1:1 dimensionless\n\tVariable_2:2 dimensionless\n");

        // Dump user-supplied values
        N_Vector rY = NULL;
        CreateVectorIfEmpty(rY, 2);
        SetVectorComponent(rY, 0, 0.0);
        SetVectorComponent(rY, 1, 1.0);

        state = ode_system.DumpState("Test 2.", rY);
        TS_ASSERT_EQUALS(state, "Test 2.\nState:\n\tVariable_1:0 dimensionless\n\tVariable_2:1 dimensionless\n");

        state = ode_system.DumpState("Test 2.", rY, 12);
        TS_ASSERT_EQUALS(state, "Test 2.\nAt independent variable (usually time) = 12\nState:\n\tVariable_1:0 dimensionless\n\tVariable_2:1 dimensionless\n");

        DeleteVector(rY);

        // For coverage of a vector helper function that makes an N_Vector from a std::vector.
        std::vector<double> std_vec;
        std_vec.push_back(1.0);
        std_vec.push_back(3.0);
        std_vec.push_back(8.0);

        N_Vector n_vec = MakeNVector(std_vec);

        TS_ASSERT_EQUALS(GetVectorSize(n_vec), GetVectorSize(std_vec));
        for (unsigned i = 0; i < GetVectorSize(n_vec); i++)
        {
            TS_ASSERT_DELTA(GetVectorComponent(n_vec, i), GetVectorComponent(std_vec, i), 1e-9);
        }

        DeleteVector(n_vec);
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestSimpleSolveUsingCvode()
    {
#ifdef CHASTE_CVODE
        CvodeFirstOrder ode_system;

        double h_value = 0.01;
        ode_system.SetMaxSteps(1000);

        // Get sampled output
        OdeSolution solutions = ode_system.Solve(0.0, 2.0, h_value, 0.1);

        int last = solutions.GetNumberOfTimeSteps();
        double testvalue = solutions.rGetSolutions()[last][0];

        // The tests
        double exact_solution = exp(2.0);

        /// \todo: #890 Work out what the global error should be bounded by
        double global_error = 1e-3;

        TS_ASSERT_DELTA(testvalue, exact_solution, global_error);

        ode_system.mUseAnalyticJacobian = true;
        ode_system.FreeCvodeMemory(); // Force a re-set

        // Here Chaste throws the error "No analytic Jacobian has been defined for this system."
        // which CVODE then re-throws with its own 'headline':
        // This covers both cases in the code.
        TS_ASSERT_THROWS_CONTAINS(ode_system.Solve(0.0, 2.0, h_value, 0.1),
                                  "CVODE failed to solve system");

#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestSequentialSolveCalls()
    {
        /*
         * All of the tests in this section pass when using Sundials >= 2.4.0
         * with 'ForcedReset' defaulting to false.
         *
         * Unfortunately they don't pass when using Sundials 2.3.0, probably because it
         * isn't as smart about going back in time when it sees changes in the RHS function.
         * (done as parameter changes here).
         *
         * So this messiness means that we have to switch on 'Forced Resetting' for
         * Sundials 2.3.0 as default to make sure it gives good answers.
         */
#ifdef CHASTE_CVODE
        {
            ParameterisedCvode ode;

            // First solve with a forced reset at each Solve call
            TS_ASSERT_EQUALS(ode.GetForceReset(), false);
            ode.SetForceReset(true);
            TS_ASSERT_EQUALS(ode.GetForceReset(), true);
            for (unsigned i = 0; i < 10; i++)
            {
                ode.SetParameter("a", i); // dy/dt = a
                ode.Solve(i, i + 1.0, 1.0);
            }
            TS_ASSERT_DELTA(ode.GetStateVariable(0u), 45.0, 1e-12);
        }

        {
            ParameterisedCvode ode;
            // (back to default behaviour)

            // But we alter time which will trigger a reset to get the right answer.
            ode.SetStateVariable(0u, 0.0);
            for (unsigned i = 0; i < 10; i++)
            {
                ode.SetParameter("a", i); // dy/dt = a
                ode.Solve(i, i + 1.0, 1.0);
            }
            TS_ASSERT_DELTA(ode.GetStateVariable(0u), 45.0, 1e-12);

            // Now we change a state variable to trigger a reset (but carry on in time).
            ode.SetStateVariable(0u, 0.0);
            ode.SetParameter("a", 1.0);
            ode.Solve(10.0, 15.0, 1.0);
            TS_ASSERT_DELTA(ode.GetStateVariable(0u), 5.0, 1e-12);

            // Now we show how you can get befuddled in MinimalReset mode.
            // carrying on in time, but not resetting the solver.
            TS_ASSERT_EQUALS(ode.GetMinimalReset(), false);
            ode.SetMinimalReset(true);
            TS_ASSERT_EQUALS(ode.GetMinimalReset(), true);
            ode.SetStateVariable(0u, 0.0);
            for (unsigned i = 0; i < 10; i++)
            {
                ode.SetParameter("a", i); // dy/dt = a
                ode.Solve(15.0 + i, i + 16.0, 1.0);
            }
#if CHASTE_SUNDIALS_VERSION >= 20400
            TS_ASSERT_DELTA(ode.GetStateVariable(0u), 50.0, 1e-12); // N.B. This is wrong!
#else
            TS_ASSERT_DELTA(ode.GetStateVariable(0u), 40.8181, 1e-4); // N.B. This is also wrong!
#endif
            // (We tricked ODE system by resetting a state variable in minimal reset mode).
        }
#else
        std::cout << "Cvode is not enabled - this test was not run.\n";
#endif // CHASTE_CVODE
    }

    void TestArchiving()
    {
#ifdef CHASTE_CVODE
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "parameterised_cvode.arch";

        double param_value;
        double var_value;
        std::string param_name;
        double answer_should_be = 45.0;

        { // Save
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            AbstractCvodeSystem* const p_ode = new ParameterisedCvode;

            TS_ASSERT_EQUALS(p_ode->GetNumberOfParameters(), 1u);
            param_value = 999;
            p_ode->SetParameter(0u, param_value);
            var_value = 9999;
            p_ode->SetStateVariable(0u, var_value);
            param_name = p_ode->rGetParameterNames()[0];
            output_arch << p_ode;

            // We copy the above sequential solve test and straddle a checkpoint
            // First solve with an automatic reset at each Solve call
            AbstractCvodeSystem* const p_ode2 = new ParameterisedCvode;
            for (unsigned i = 0; i < 5; i++)
            {
                p_ode2->SetParameter("a", i); // dy/dx = i
                p_ode2->Solve(i, i + 1.0, 1.0);
            }
            output_arch << p_ode2;
            for (unsigned i = 5; i < 10; i++)
            {
                p_ode2->SetParameter("a", i); // dy/dx = i
                p_ode2->Solve(i, i + 1.0, 1.0);
            }
            TS_ASSERT_DELTA(p_ode2->GetStateVariable(0u), answer_should_be, 1e-12);

            delete p_ode;
            delete p_ode2;
        }
        { // Load
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            AbstractCvodeSystem* p_ode;
            input_arch >> p_ode;

            TS_ASSERT_EQUALS(p_ode->GetSystemName(), "ParameterisedCvode");
            TS_ASSERT_EQUALS(p_ode->GetParameter(0), param_value);
            TS_ASSERT_EQUALS(p_ode->GetStateVariable(0), var_value);
            TS_ASSERT_EQUALS(p_ode->rGetParameterNames()[0], param_name);
            TS_ASSERT_EQUALS(p_ode->GetNumberOfParameters(), 1u);

            AbstractCvodeSystem* p_ode2;
            input_arch >> p_ode2;

            // Run for a while after checkpoint, check it agrees with the above test and
            // the ode system that was saved.
            for (unsigned i = 5; i < 10; i++)
            {
                p_ode2->SetParameter("a", i); // dy/dx = i
                p_ode2->Solve(i, i + 1.0, 1.0);
            }
            TS_ASSERT_DELTA(p_ode2->GetStateVariable(0u), answer_should_be, 1e-12);

            delete p_ode;
            delete p_ode2;
        }
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }
};

#endif //_TESTABSTRACTODESYSTEM_HPP_
