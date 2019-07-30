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
#ifndef TESTCOMBINEDODESYSTEM_HPP_
#define TESTCOMBINEDODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "CombinedOdeSystem.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "OdeSolution.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * dy/dt = x, y(0) = 0. Here x is a parameter.
 */
class SimpleOde1 : public AbstractOdeSystem
{
public:
    SimpleOde1() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde1>::Instance();
        ResetToInitialConditions();
        mParameters.resize(1);
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = mParameters[0];
    }
};

template<>
void OdeSystemInformation<SimpleOde1>::Initialise()
{
    this->mVariableNames.push_back("Variable_y");
    this->mVariableUnits.push_back("Units_y");
    this->mInitialConditions.push_back(0.0);

    this->mParameterNames.push_back("Variable_x");
    this->mParameterUnits.push_back("Units_x");

    this->mInitialised = true;
}


/**
 * dx/dt = -y, x(0) = 1. Here y is a parameter.
 */
class SimpleOde2 : public AbstractOdeSystem
{
public:
    SimpleOde2() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde2>::Instance();
        ResetToInitialConditions();
        mParameters.resize(1);
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = -mParameters[0];
    }
};

template<>
void OdeSystemInformation<SimpleOde2>::Initialise()
{
    this->mVariableNames.push_back("Variable_x");
    this->mVariableUnits.push_back("Units_x");
    this->mInitialConditions.push_back(1.0);

    this->mParameterNames.push_back("Variable_y");
    this->mParameterUnits.push_back("Units_y");

    this->mInitialised = true;
}

/**
 * The following classes are used in the solution of
 * x'=x-y+z, y'=y-z, z'=2y-z
 * starting at (x,y,z)=(0,1,0)
 * Analytic solution is x=-sin(t), y=sin(t)+cos(t), z=2sin(t)
 */

/**
 * dx/dt = x-y+z, x(0) = 0. Here y & z are parameters.
 */
class SimpleOde3 : public AbstractOdeSystem
{
public:
    SimpleOde3() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde3>::Instance();
        ResetToInitialConditions();
        mParameters.resize(2);
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0] - mParameters[0] + mParameters[1];
    }
};

template<>
void OdeSystemInformation<SimpleOde3>::Initialise()
{
    this->mVariableNames.push_back("x");
    this->mVariableUnits.push_back("dimensionless_x");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}

/**
 * dy/dt = y-z, y(0) = 1. Here z is a parameter.
 */
class SimpleOde4 : public AbstractOdeSystem
{
public:
    SimpleOde4() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde4>::Instance();
        ResetToInitialConditions();
        mParameters.resize(1);
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0] - mParameters[0];
    }
};

template<>
void OdeSystemInformation<SimpleOde4>::Initialise()
{
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless_y");
    this->mInitialConditions.push_back(1.0);

    this->mInitialised = true;
}


/**
 * dz/dt = 2y-z, z(0) = 0. Here y is a parameter.
 */
class SimpleOde5 : public AbstractOdeSystem
{
public:
    SimpleOde5() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde5>::Instance();
        ResetToInitialConditions();
        mParameters.resize(1);
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = 2.0*mParameters[0] - rY[0];
    }
};

template<>
void OdeSystemInformation<SimpleOde5>::Initialise()
{
    this->mVariableNames.push_back("z");
    this->mVariableUnits.push_back("dimensionless_z");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}


/**
 * Double system (combination of SimpleOde4, for y' and SimpleOde5, for z') to be used with SimpleOde3
 * dy/dt = y-z, y(0) = 1.
 * dz/dt = 2y-z, z(0) = 0.
 */
class SimpleOde6 : public AbstractOdeSystem
{
public:
    SimpleOde6() : AbstractOdeSystem(2) // 2 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde6>::Instance();
        ResetToInitialConditions();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0] - rY[1];        //y'=y-z
        rDY[1] = 2.0*rY[0] - rY[1];    //z'=2y-z
    }
};

template<>
void OdeSystemInformation<SimpleOde6>::Initialise()
{
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("z");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}


class TestCombinedOdeSystem : public CxxTest::TestSuite
{
public:

    void TestSimpleCombinedOdeSystem()
    {
        // Create two ODE systems
        SimpleOde1 ode_for_y; // dy/dt = x
        SimpleOde2 ode_for_x; // dx/dt = -y

        std::vector<AbstractOdeSystem*> ode_systems;
        ode_systems.push_back(&ode_for_y);
        ode_systems.push_back(&ode_for_x);

        // Create combined ODE system
        CombinedOdeSystem combined_ode_system(ode_systems);

        // Tell the combined ODE system which state variables in the first ODE system
        // correspond to which parameters in the second ODE system...
        std::map<unsigned, unsigned> variable_parameter_map;
        variable_parameter_map[0] = 0;

        combined_ode_system.Configure(variable_parameter_map, &ode_for_y, &ode_for_x);

        // ...and vice versa (we can re-use the map in this case)
        combined_ode_system.Configure(variable_parameter_map, &ode_for_x, &ode_for_y);

        // Test number of state variables
        unsigned num_variables = combined_ode_system.GetNumberOfStateVariables();
        TS_ASSERT_EQUALS(num_variables, 2u);

        // Combined system has no parameters
        TS_ASSERT_EQUALS(combined_ode_system.GetNumberOfParameters(), 0u);
        TS_ASSERT_EQUALS(combined_ode_system.rGetParameterNames().size(), 0u);

        // Test initial conditions
        std::vector<double> initial_conditions = combined_ode_system.GetInitialConditions();
        TS_ASSERT_DELTA(initial_conditions[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 1.0, 1e-12);
        // Test variable names & units
        const std::vector<std::string>& r_names = combined_ode_system.rGetStateVariableNames();
        TS_ASSERT_EQUALS(r_names[0], ode_for_y.rGetStateVariableNames()[0]);
        TS_ASSERT_EQUALS(r_names[1], ode_for_x.rGetStateVariableNames()[0]);
        const std::vector<std::string>& r_units = combined_ode_system.rGetStateVariableUnits();
        TS_ASSERT_EQUALS(r_units[0], ode_for_y.rGetStateVariableUnits()[0]);
        TS_ASSERT_EQUALS(r_units[1], ode_for_x.rGetStateVariableUnits()[0]);

        // Test solving the combined system.
        // This is dy/dt = x, dx/dt = -y, y(0) = 0, x(0) = 1.
        // The analytic solution is y = sin(t), x = cos(t).
        EulerIvpOdeSolver solver;
        OdeSolution solutions;
        double h = 0.01;
        std::vector<double> inits = combined_ode_system.GetInitialConditions();
        solutions = solver.Solve(&combined_ode_system, inits, 0.0, 2.0, h, h);
        double global_error = 0.5*(exp(2.0)-1)*h;
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[0], sin(2.0), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[1], cos(2.0), global_error);

        // Check that if we create the same combination, we get the same information object
        boost::shared_ptr<const AbstractOdeSystemInformation> info1 = combined_ode_system.GetSystemInformation();
        CombinedOdeSystem combined_ode_system2(ode_systems);
        boost::shared_ptr<const AbstractOdeSystemInformation> info2 = combined_ode_system2.GetSystemInformation();
        TS_ASSERT_EQUALS(info1, info2);
    }

    void TestSimpleSystemWithOrderSwapped()
    {
        // The solution should be the same, but we'll have to construct a new CombinedOdeSystemInformation
        // object, because the order of subsystems has changed.

        SimpleOde1 ode_for_y; // dy/dt = x
        SimpleOde2 ode_for_x; // dx/dt = -y

        std::vector<AbstractOdeSystem*> ode_systems;
        ode_systems.push_back(&ode_for_x);
        ode_systems.push_back(&ode_for_y);

        // Create combined ODE system
        CombinedOdeSystem combined_ode_system(ode_systems);

        // Tell the combined ODE system which state variables in the first ODE system
        // correspond to which parameters in the second ODE system...
        std::map<unsigned, unsigned> variable_parameter_map;
        variable_parameter_map[0] = 0;

        combined_ode_system.Configure(variable_parameter_map, &ode_for_y, &ode_for_x);

        // ...and vice versa (we can re-use the map in this case)
        combined_ode_system.Configure(variable_parameter_map, &ode_for_x, &ode_for_y);

        // Test solving the combined system.
        // This is dy/dt = x, dx/dt = -y, y(0) = 0, x(0) = 1.
        // The analytic solution is y = sin(t), x = cos(t).
        EulerIvpOdeSolver solver;
        OdeSolution solutions;
        double h = 0.01;
        std::vector<double> inits = combined_ode_system.GetInitialConditions();
        solutions = solver.Solve(&combined_ode_system, inits, 0.0, 2.0, h, h);
        double global_error = 0.5*(exp(2.0)-1)*h;
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[1], sin(2.0), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[0], cos(2.0), global_error);
    }

    void TestWithThreeVariables()
    {
        SimpleOde3 ode_for_x; // dx/dt = x -y +z
        SimpleOde4 ode_for_y; // dy/dt = y-z
        SimpleOde5 ode_for_z; // dz/dt = 2y-z

        std::vector<AbstractOdeSystem*> ode_systems;
        ode_systems.push_back(&ode_for_x);
        ode_systems.push_back(&ode_for_y);
        ode_systems.push_back(&ode_for_z);

        // Create combined ODE system
        CombinedOdeSystem combined_ode_system(ode_systems);

        // Tell the combined ODE system which state variables in the first ODE system
        // correspond to which parameters in the second ODE system...
        std::map<unsigned, unsigned> variable_parameter_map_yx;
        variable_parameter_map_yx[0] = 0; //y in the y-ODE appears in the x-ODE
        combined_ode_system.Configure(variable_parameter_map_yx, &ode_for_y, &ode_for_x);

        std::map<unsigned, unsigned> variable_parameter_map_zx;
        variable_parameter_map_zx[0] = 1; //z in the z-ODE appears in the x-ODE
        combined_ode_system.Configure(variable_parameter_map_zx, &ode_for_z, &ode_for_x);

        //Reuse the variable_parameter_map_yx
        //y in the y-ODE appears in the z-ODE
        combined_ode_system.Configure(variable_parameter_map_yx, &ode_for_y, &ode_for_z);
        //z in the z-ODE appears in the y-ODE
        combined_ode_system.Configure(variable_parameter_map_yx, &ode_for_z, &ode_for_y);


        // Test number of state variables
        unsigned num_variables = combined_ode_system.GetNumberOfStateVariables();
        TS_ASSERT_EQUALS(num_variables, 3u);

        // Combined system has no parameters
        TS_ASSERT_EQUALS(combined_ode_system.GetNumberOfParameters(), 0u);
        TS_ASSERT_EQUALS(combined_ode_system.rGetParameterNames().size(), 0u);

        // Test initial conditions
        std::vector<double> initial_conditions = combined_ode_system.GetInitialConditions();
        TS_ASSERT_DELTA(initial_conditions[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 1.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[2], 0.0, 1e-12);
        // Test variable names & units
        const std::vector<std::string>& r_names = combined_ode_system.rGetStateVariableNames();
        TS_ASSERT_EQUALS(r_names[0], ode_for_x.rGetStateVariableNames()[0]);
        TS_ASSERT_EQUALS(r_names[1], ode_for_y.rGetStateVariableNames()[0]);
        TS_ASSERT_EQUALS(r_names[2], ode_for_z.rGetStateVariableNames()[0]);
        const std::vector<std::string>& r_units = combined_ode_system.rGetStateVariableUnits();
        TS_ASSERT_EQUALS(r_units[0], ode_for_x.rGetStateVariableUnits()[0]);
        TS_ASSERT_EQUALS(r_units[1], ode_for_y.rGetStateVariableUnits()[0]);
        TS_ASSERT_EQUALS(r_units[2], ode_for_z.rGetStateVariableUnits()[0]);
        TS_ASSERT_EQUALS(r_units[0], "dimensionless_x");
        TS_ASSERT_EQUALS(r_units[1], "dimensionless_y");
        TS_ASSERT_EQUALS(r_units[2], "dimensionless_z");

        // x'=x-y+z, y'=y-z, z'=2y-z
        // starting at (x,y,z)=(0,1,0)
        // Analytic solution is x=-sin(t), y=sin(t)+cos(t), z=2sin(t)
        EulerIvpOdeSolver solver;
        OdeSolution solutions;
        double h = 0.01;
        std::vector<double> inits = combined_ode_system.GetInitialConditions();
        solutions = solver.Solve(&combined_ode_system, inits, 0.0, 2.0, h, h);
        double global_error = 0.5*(exp(2.0)-1)*h;
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[0], -sin(2.0), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[1], sin(2.0)+cos(2.0), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[2], 2.0*sin(2.0), global_error);
    }

    void TestWithThreeVariablesTwoSystems()
    {
        SimpleOde3 ode_for_x; // dx/dt = x -y +z
        SimpleOde6 ode_for_yz; // dy/dt = y-z  and dz/dt = 2y-z

        std::vector<AbstractOdeSystem*> ode_systems;
        ode_systems.push_back(&ode_for_x);
        ode_systems.push_back(&ode_for_yz);

        // Create combined ODE system
        CombinedOdeSystem combined_ode_system(ode_systems);

        // Tell the combined ODE system which state variables in the first ODE system
        // correspond to which parameters in the second ODE system...
        std::map<unsigned, unsigned> variable_parameter_map;
        variable_parameter_map[0] = 0; //y in the yz-ODE appears in the x-ODE
        variable_parameter_map[1] = 1; //z in the yz-ODE appears in the x-ODE
        combined_ode_system.Configure(variable_parameter_map, &ode_for_yz, &ode_for_x);

        // Test number of state variables
        unsigned num_variables = combined_ode_system.GetNumberOfStateVariables();
        TS_ASSERT_EQUALS(num_variables, 3u);

        // Combined system has no parameters
        TS_ASSERT_EQUALS(combined_ode_system.GetNumberOfParameters(), 0u);
        TS_ASSERT_EQUALS(combined_ode_system.rGetParameterNames().size(), 0u);

        // Test initial conditions
        std::vector<double> initial_conditions = combined_ode_system.GetInitialConditions();
        TS_ASSERT_DELTA(initial_conditions[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 1.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[2], 0.0, 1e-12);
        // Test variable names & units
        const std::vector<std::string>& r_names = combined_ode_system.rGetStateVariableNames();
        TS_ASSERT_EQUALS(r_names[0], ode_for_x.rGetStateVariableNames()[0]);
        TS_ASSERT_EQUALS(r_names[1], ode_for_yz.rGetStateVariableNames()[0]);
        TS_ASSERT_EQUALS(r_names[2], ode_for_yz.rGetStateVariableNames()[1]);

        // x'=x-y+z, y'=y-z, z'=2y-z
        // starting at (x,y,z)=(0,1,0)
        // Analytic solution is x=-sin(t), y=sin(t)+cos(t), z=2sin(t)
        EulerIvpOdeSolver solver;
        OdeSolution solutions;
        double h = 0.01;
        std::vector<double> inits = combined_ode_system.GetInitialConditions();
        solutions = solver.Solve(&combined_ode_system, inits, 0.0, 2.0, h, h);
        double global_error = 0.5*(exp(2.0)-1)*h;
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[0], -sin(2.0), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[1], sin(2.0)+cos(2.0), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[2], 2.0*sin(2.0), global_error);
    }
};

#endif /*TESTCOMBINEDODESYSTEM_HPP_*/
