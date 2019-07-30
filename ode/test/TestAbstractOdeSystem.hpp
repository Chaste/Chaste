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


#ifndef _TESTABSTRACTODESYSTEM_HPP_
#define _TESTABSTRACTODESYSTEM_HPP_

#include <cmath>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Ode1.hpp"
#include "Ode2.hpp"
#include "Ode3.hpp"
#include "TwoDimOdeSystem.hpp"
#include "VanDerPolOde.hpp"
#include "ParameterisedOde.hpp"
#include "OdeSystemForCoupledHeatEquation.hpp"

#include "OutputFileHandler.hpp"

#include "FakePetscSetup.hpp"

// Tolerance for tests
const double tol = 0.01;


class TestAbstractOdeSystem : public CxxTest::TestSuite
{
public:

    void TestOdeSystemOne()
    {
        // Test Ode1 class
        Ode1 ode1;
        // dy
        std::vector<double> dy(1);
        ode1.EvaluateYDerivatives(1.0, ode1.GetInitialConditions(),dy);
        TS_ASSERT_DELTA(dy[0], 1.0, tol);

        // System name & free var info shouldn't fall over if there isn't any
        boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode1.GetSystemInformation();

        TS_ASSERT_EQUALS(ode1.GetSystemName(), "");
        TS_ASSERT_EQUALS(p_info->GetSystemName(), "");

        TS_ASSERT_EQUALS(p_info->GetFreeVariableName(), "");
        TS_ASSERT_EQUALS(p_info->GetFreeVariableUnits(), "");

    }

    void TestOdeSystemTwo()
    {
        Ode2 ode2;
        std::vector<double> dy(1);
        ode2.EvaluateYDerivatives(2.0, ode2.GetInitialConditions(), dy);
        TS_ASSERT_DELTA(dy[0], 8.0, tol);
    }

    void TestOdeSystemThree()
    {
        Ode3 ode3;
        std::vector<double> dy(2);
        ode3.EvaluateYDerivatives(2.0, ode3.GetInitialConditions(), dy);
        TS_ASSERT_DELTA(dy[0], 8.0, tol);
        TS_ASSERT_DELTA(dy[1], 16.0, tol);
    }

    void TestExceptions()
    {
        Ode1 ode;
        TS_ASSERT_EQUALS(ode.GetNumberOfStateVariables(), 1u);
        std::vector<double> v(2);
        v[0] = -1.0;
        v[1] = -2.0;
        TS_ASSERT_THROWS_THIS(ode.SetDefaultInitialConditions(v),
                "The number of initial conditions must be that of the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.SetDefaultInitialCondition(2, -3.0),
                "Index is greater than the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.SetStateVariables(v),
                "The size of the passed in vector must be that of the number of state variables.");
    }

    void TestParameters()
    {
        ParameterisedOde ode;
        boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();

        TS_ASSERT_EQUALS(ode.GetSystemName(), "ParameterisedOde");
        TS_ASSERT_EQUALS(p_info->GetSystemName(), "ParameterisedOde");

        TS_ASSERT_EQUALS(p_info->GetFreeVariableName(), "time");
        TS_ASSERT_EQUALS(p_info->GetFreeVariableUnits(), "ms");

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

        ode.SetStateVariable("y", y-1);
        TS_ASSERT_EQUALS(ode.GetStateVariable("y"), y-1);
        ode.SetAnyVariable("y", y);
        TS_ASSERT_EQUALS(ode.GetStateVariable("y"), y);


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
    }

    void TestAttributes()
    {
        ParameterisedOde ode;
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
    }


    void TestArchivingOfParameters()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "parameterised_ode.arch";
        double param_value;
        std::string param_name;
        { // Save
            AbstractOdeSystem * const p_ode = new ParameterisedOde;

            TS_ASSERT_EQUALS(p_ode->GetNumberOfParameters(), 1u);
            param_value = p_ode->GetParameter(0);
            param_name = p_ode->rGetParameterNames()[0];

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_ode;
            delete p_ode;
        }
        { // Normal load
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            AbstractOdeSystem* p_ode;
            input_arch >> p_ode;

            TS_ASSERT_EQUALS(p_ode->GetSystemName(), "ParameterisedOde");
            TS_ASSERT_EQUALS(p_ode->GetParameter(0), param_value);
            TS_ASSERT_EQUALS(p_ode->rGetParameterNames()[0], param_name);
            TS_ASSERT_EQUALS(p_ode->GetNumberOfParameters(), 1u);

            delete p_ode;
        }
        { // Load with param name changed
            ParameterisedOde ode;
            boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();
            AbstractOdeSystemInformation* p_mod_info = const_cast<AbstractOdeSystemInformation*>(p_info.get());
            std::string new_name("new_param_name");
            TS_ASSERT_DIFFERS(p_mod_info->mParameterNames[0], new_name);
            p_mod_info->mParameterNames[0] = new_name;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            AbstractOdeSystem* p_ode = NULL;  // Shouldn't be necessary to set to NULL but seems to be!
            TS_ASSERT_THROWS_CONTAINS(input_arch >> p_ode, "Archive specifies a parameter 'a' which does not appear in this class.");
            // Mend the ODE system info for the following tests.
            p_mod_info->mParameterNames[0] = param_name;
            if (p_ode) // It seems some setups will delete if the throw happens and some won't!
            {
                delete p_ode;
            }
        }
        { // Load with a parameter added
            ParameterisedOde ode;
            boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();
            AbstractOdeSystemInformation* p_mod_info = const_cast<AbstractOdeSystemInformation*>(p_info.get());

            p_mod_info->mParameterNames.push_back("new_name");

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            AbstractOdeSystem* p_ode = NULL; // The = NULL obviously shouldn't be necessary, but seems to be on some compilers!
            TS_ASSERT_THROWS_CONTAINS(input_arch >> p_ode, "Number of ODE parameters in archive does not match number in class.");
            // Mend the ODE system info for the following tests.
            p_mod_info->mParameterNames.resize(1u);
            if (p_ode) // It seems some setups will delete if the throw happens and some won't!
            {
                delete p_ode;
            }
        }
        { // Load with a parameter added, and the constructor providing a default
            ParameterisedOde ode;
            boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();
            AbstractOdeSystemInformation* p_mod_info = const_cast<AbstractOdeSystemInformation*>(p_info.get());

            p_mod_info->mParameterNames.push_back("new_name");

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            AbstractOdeSystem* p_ode;
            ParameterisedOde::fakeSecondParameter = true;
            input_arch >> p_ode;
            delete p_ode;
            // Mend the ODE system info for the following tests.
            ParameterisedOde::fakeSecondParameter = false;
            p_mod_info->mParameterNames.resize(1u);
        }
        { // Load with no defaults provided by the constructor
            ParameterisedOde::noParameterDefaults = true;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            AbstractOdeSystem* p_ode;
            input_arch >> p_ode;
            delete p_ode;
            // Mend the ODE system info for the following tests.
            ParameterisedOde::noParameterDefaults = false;
        }
    }

    void TestDerivedQuantities()
    {
        ParameterisedOde ode;
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

        std::vector<double> derived = ode.ComputeDerivedQuantitiesFromCurrentState(0.0);
        double a = ode.GetParameter(0);
        TS_ASSERT_EQUALS(a, 0.0);
        TS_ASSERT_DELTA(derived[0], 2*a, 1e-4);
        TS_ASSERT_DELTA(ode.GetAnyVariable(2u, 0.0), 2*a, 1e-4);
        TS_ASSERT_DELTA(ode.GetAnyVariable(2u, 0.0, &derived), 2*a, 1e-4);
        a = 1.0;
        ode.SetParameter(0, a);
        derived = ode.ComputeDerivedQuantities(0.0, ode.GetInitialConditions());
        TS_ASSERT_DELTA(derived[0], 2*a, 1e-4);
        double y = 10.0;
        ode.SetStateVariable(0, y);
        derived = ode.ComputeDerivedQuantitiesFromCurrentState(0.0);
        TS_ASSERT_DELTA(derived[0], 2*a+y, 1e-4);
        TS_ASSERT_DELTA(ode.GetAnyVariable(2u, 1.0/* ignored for this ODE */), 2*a+y, 1e-4);

        // Exceptions
        TS_ASSERT_THROWS_THIS(ode.GetDerivedQuantityIndex("Missing"), "No derived quantity named 'Missing'.");
        TS_ASSERT_THROWS_THIS(ode.GetDerivedQuantityUnits(1u), "The index passed in must be less than the number of derived quantities.");

        TwoDimOdeSystem ode2;
        TS_ASSERT_THROWS_THIS(ode2.ComputeDerivedQuantitiesFromCurrentState(0.0),
                              "This ODE system does not define derived quantities.");
        std::vector<double> doesnt_matter;
        TS_ASSERT_THROWS_THIS(ode2.ComputeDerivedQuantities(0.0, doesnt_matter),
                              "This ODE system does not define derived quantities.");
    }

    void TestSetGetFunctionsInAbstractOdeSystem()
    {
        TwoDimOdeSystem ode;

        std::vector<double> initial_conditions = ode.GetInitialConditions();
        std::vector<double>& r_state_variables = ode.rGetStateVariables();

        TS_ASSERT_DELTA(initial_conditions[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 2.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[0], 3.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[1], 4.0, 1e-12);

        std::vector<double> new_initial_conditions;
        new_initial_conditions.push_back(5.0);
        new_initial_conditions.push_back(6.0);

        std::vector<double> new_state_variables;
        new_state_variables.push_back(7.0);
        new_state_variables.push_back(8.0);

        ode.SetDefaultInitialConditions(new_initial_conditions);
        ode.SetStateVariables(new_state_variables);

        initial_conditions = ode.GetInitialConditions();

        TS_ASSERT_DELTA(initial_conditions[0], 5.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 6.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[0], 7.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[1], 8.0, 1e-12);

        ode.SetDefaultInitialCondition(1, 9.0);
        initial_conditions = ode.GetInitialConditions();
        TS_ASSERT_DELTA(initial_conditions[0], 5.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 9.0, 1e-12);

        // Archive the ODE system
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "ode.arch";
        std::ofstream ofs(archive_filename.c_str());
        boost::archive::text_oarchive output_arch(ofs);

        output_arch <<  static_cast<const TwoDimOdeSystem&>(ode);
        TS_ASSERT_DELTA(r_state_variables[0], 7.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[1], 8.0, 1e-12);

        ode.SetStateVariable(0, 2.0);
        ode.SetStateVariable(1, 5.0);

        TS_ASSERT_THROWS_THIS(ode.SetStateVariable(2, 1.0),
                "The index passed in must be less than the number of state variables."); //cover exception

        TS_ASSERT_DELTA(r_state_variables[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[1], 5.0, 1e-12);
    }

    void TestLoadAbstractOdeSystem()
    {
        TwoDimOdeSystem ode;

        TS_ASSERT_EQUALS( ode.GetNumberOfStateVariables(), 2U );

        // Read archive from previous test
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "ode.arch";

        std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
        boost::archive::text_iarchive input_arch(ifs);

        input_arch >> ode;

        TS_ASSERT_EQUALS( ode.GetNumberOfStateVariables(), 2U );

        std::vector<double>& r_state_variables = ode.rGetStateVariables();

        TS_ASSERT_DELTA(r_state_variables[0], 7.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[1], 8.0, 1e-12);
    }

    void TestReadSpecificStateVariable()
    {
        // Create a VanDerPol system
        VanDerPolOde ode_system;

        // get the velocity state variable number
        unsigned var_number = ode_system.GetStateVariableIndex("v");
        TS_ASSERT_EQUALS(var_number, 1u);

        TS_ASSERT_THROWS_THIS(ode_system.GetStateVariableIndex("foo"),
                              "No state variable named 'foo'.");

        TS_ASSERT_EQUALS(ode_system.GetStateVariable(var_number), 10.0);

        TS_ASSERT_EQUALS(ode_system.GetStateVariableUnits(var_number), "m/s");
    }

    // This test is mainly for coverage purposes.
    void TestDumpState()
    {
        // Create a VanDerPol system
        VanDerPolOde ode_system;

        // Dump the state variables
        std::string state = ode_system.DumpState("This is a test.");
        TS_ASSERT_EQUALS(state, "This is a test.\nState:\n\tx:10 m\n\tv:10 m/s\n");

        // Dump user-supplied values
        std::vector<double> rY(2);
        rY[0] = 0.0;
        rY[1] = 1.0;
        state = ode_system.DumpState("Test 2.", rY, 0);
        TS_ASSERT_EQUALS(state, "Test 2.\nAt independent variable (usually time) = 0\nState:\n\tx:0 m\n\tv:1 m/s\n");

        state = ode_system.DumpState("Test 2.", rY);
        TS_ASSERT_EQUALS(state, "Test 2.\nState:\n\tx:0 m\n\tv:1 m/s\n");
    }

    void TestAbstractOdeSystemForCoupledPdeSystem()
    {
        // Create ODE system
        OdeSystemForCoupledHeatEquation ode_system(1.0);

        // Test setting/getting a PDE solution vector with the correct size
        std::vector<double> solution(1);
        solution[0] = 12.3;

        ode_system.SetPdeSolution(solution);

        TS_ASSERT_EQUALS(ode_system.GetPdeSolutionSize(), 1u);

        std::vector<double> obtained_solution = ode_system.rGetPdeSolution();
        TS_ASSERT_DELTA(obtained_solution[0], solution[0], 1e-6);

        // Test an exception is thrown if setting a PDE solution vector with incorrect size
        std::vector<double> bad_solution(2);
        bad_solution[0] = 2.1;
        bad_solution[1] = 4.9;

        TS_ASSERT_THROWS_THIS(ode_system.SetPdeSolution(bad_solution),
                              "The supplied vector is not the correct size.");
    }
};

#endif //_TESTABSTRACTODESYSTEM_HPP_
