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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTSOLVINGODESTUTORIAL_HPP_
#define TESTSOLVINGODESTUTORIAL_HPP_
/*
 * = In this tutorial we show how Chaste can be used to solve an ODE system =
 *
 * The following header files need to be included.
 * First we include the header needed to define this class as a test suite. */
#include <cxxtest/TestSuite.h>

/*
 * In some early versions of Boost this file has to be included first so that the name of the ODE solver
 * used can be written to result files.
 */
#include "CheckpointArchiveTypes.hpp"
/*
 * We will use a simple forward Euler solver to solve the ODE, so the following
 * needs to be included.
 */
#include "EulerIvpOdeSolver.hpp"
/* All the ODE solvers take in a concrete ODE system class, which is user-defined
 * and must inherit from the following class, which defines an ODE interface.
 */
#include "AbstractOdeSystem.hpp"
/* In order to define useful information about the ODE system conveniently, such
 * as the names and units of variables, and suggested initial conditions, we
 * need the following header.
 */
#include "OdeSystemInformation.hpp"

/* This test doesn't support being run on multiple processes, so we need this header
 * to prevent race conditions when writing files.
 */
#include "FakePetscSetup.hpp"

/*
 * == Defining the ODE classes ==
 *
 * Let us solve the ODE dy/dt = y^2^+t^2^, with y(0) = 1. To do so, we have to define
 * our own ODE class, inheriting from {{{AbstractOdeSystem}}}, which implements the
 * {{{EvaluateYDerivatives()}}} method.
 */
class MyOde : public AbstractOdeSystem
{
public:
/* The constructor does very little.
 * It calls the base constructor, passing the number of state variables in the
 * ODE system (here, 1, i.e. y is a 1d vector).
 * It also sets the object to use to retrieve system information (see later).
 */
    MyOde() : AbstractOdeSystem(1)
    {
        mpSystemInfo = OdeSystemInformation<MyOde>::Instance();
    }

/* The ODE solvers will repeatedly call a method called `EvaluateYDerivatives`, which needs
 * to be implemented in this concrete class. This takes in the time, a {{{std::vector}}} of
 * y values (in this case, of size 1), and a reference to a {{{std::vector}}} in which the
 * derivative(s) should be filled in by the method...
 */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY)
    {
        /*...so we set {{{rDY[0]}}} to be y^2^ + t^2^. */
        rDY[0] = rY[0]*rY[0] + time*time;
    }
};

/* The following ''template specialisation'' defines the information for this
 * ODE system.  Note that we use the ODE system class that we have just defined
 * as a template parameter
 */
template<>
void OdeSystemInformation<MyOde>::Initialise()
{
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}

/* That would be all that is needed to solve this ODE. However, rather
 * than solving up to a fixed time, suppose we wanted to solve until some function
 * of y (and t) reached a certain value, e.g. let's say we wanted to solve the ODE until
 * y reached 2.5. To do this, we have to define a stopping event, by overriding
 * the method {{{CalculateStoppingEvent}}} from {{{AbstractOdeSystem}}}. For this, let us
 * define a new class, inheriting from the above class (i.e. representing the same ODE)
 * but with a stopping event defined.
 *
 * Note that we do not define separate ODE system information for this class - it uses
 * that defined in the base class `MyOde`, since it represents the same ODE.
 */
class MyOdeWithStoppingEvent : public MyOde
{
public:
    /* All we have to do is implement the following function. This is defined in
     * the base class ({{{AbstractOdeSystem}}}), where it always returns false, and here we override it
     * to return true if y>=2.5
     */
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY)
    {
        return (rY[0]>=2.5);
    }
};

/* The following class will make more sense when solving with state variables is discussed.
 * It is another ODE class which sets up a 'state variable'. Note that this is done in the
 * constructor, and the {{{EvaluateYDerivatives}}} method is identical to before. */
class MyOdeUsingStateVariables : public AbstractOdeSystem
{
public:
    MyOdeUsingStateVariables() : AbstractOdeSystem(1)
    {
        mpSystemInfo = OdeSystemInformation<MyOdeUsingStateVariables>::Instance();
        mStateVariables.push_back(1.0);
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY)
    {
        rDY[0] = rY[0]*rY[0] + time*time;
    }
};

/* This time we do need to define the ODE system information.
 */
template<>
void OdeSystemInformation<MyOdeUsingStateVariables>::Initialise()
{
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mInitialised = true;
}

/* This class is another simple ODE class, just as an example of how a 2d ODE is solved. Here
 * we solve the ODE dy,,1,,/dt = y,,2,,, dy,,2,,/dt = (y,,1,,)^2^ (which represents the second-order ODE d^2^y/dt^2^ = y^2^).
 */
class My2dOde : public AbstractOdeSystem
{
public:
    My2dOde() : AbstractOdeSystem(2)
    {
        mpSystemInfo = OdeSystemInformation<My2dOde>::Instance();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY)
    {
        rDY[0] = rY[1];
        rDY[1] = rY[0]*rY[0];
    }
};

/* Again we need to define the ODE system information.
 */
template<>
void OdeSystemInformation<My2dOde>::Initialise()
{
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("ydot");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}

/*
 * == The Tests ==
 *
 * === Standard ODE solving ===
 *
 * Now we can define the test, in which the ODEs are solved. */
class TestSolvingOdesTutorial: public CxxTest::TestSuite
{
public:
    void TestSolvingOdes()
    {
        /* First, create an instance of the ODE class to be solved. */
        MyOde my_ode;
        /* Next, create a solver. */
        EulerIvpOdeSolver euler_solver;
        /* We will need to provide an initial condition, which needs to
         * be a {{{std::vector}}}.*/
        std::vector<double> initial_condition;
        initial_condition.push_back(1.0);
        /* Then, just call `Solve`, passing in a pointer to the ODE, the
         * initial condition, the start time, end time, the solving timestep,
         * and sampling timestep (how often we want the solution stored in the returned `OdeSolution` object).
         * Here we solve from 0 to 1, with a timestep of 0.01 but a ''sampling
         * timestep'' of 0.1. The return value is an object of type {{{OdeSolution}}}
         * (which is basically just a list of times and solutions).
         */
        OdeSolution solutions = euler_solver.Solve(&my_ode, initial_condition, 0, 1, 0.01, 0.1);
        /* Let's look at the results, which can be obtained from the {{{OdeSolution}}}
         * object using the methods {{{rGetTimes()}}} and {{{rGetSolutions()}}}, which
         * return a {{{std::vector}}} and a {{{std::vector}}} of {{{std::vector}}}s
         * respectively. */
        for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            /* The {{{[0]}}} here is because we are getting the zeroth component of y (a 1-dimensional vector). */
            std::cout << solutions.rGetTimes()[i] << " " << solutions.rGetSolutions()[i][0] << "\n";
        }

        /* Alternatively, we can print the solution directly to a file, using the {{{WriteToFile}}}
         * method on the {{{OdeSolution}}} class. */
        solutions.WriteToFile("SolvingOdesTutorial", "my_ode_solution", "sec");
        /* Two files are written
         * * {{{my_ode_solution.dat}}} contains the results (a header line, then one column for time and one column per variable)
         * * {{{my_ode_solution.info}}} contains information for reading the data back, a line about the ODE solver ("{{{ODE SOLVER: EulerIvpOdeSolver}}}") and provenance information.
         */

        /* We can see from the printed out results that y goes above 2.5 somewhere just
         * before 0.6. To solve only up until y=2.5, we can solve the ODE that has the
         * stopping event defined, using the same solver as before. */
        MyOdeWithStoppingEvent my_ode_stopping;

        /* '''Note:''' ''when a {{{std::vector}}} is passed in as an initial condition
         * to a {{{Solve}}} call, it gets updated as the solve takes place''. Therefore, if
         * we want to use the same initial condition again, we have to reset it back to 1.0. */
        initial_condition[0] = 1.0;
        solutions = euler_solver.Solve(&my_ode_stopping, initial_condition, 0, 1, 0.01, 0.1);
        /* We can check with the solver that it stopped because of the stopping event, rather than because
         * it reached to end time. */
        TS_ASSERT(euler_solver.StoppingEventOccurred());
        /* Finally, let's print the time of the stopping event (to the nearest dt or so). */
        std::cout << "Stopping event occurred at t="<<solutions.rGetTimes().back()<<"\n";
    }

    /*
     * === ODE solving using state variables ===
     *
     * In this second test, we show how to do an alternative version of ODE solving, which
     * does not involve passing in initial conditions and returning an {{{OdeSolution}}}.
     * The {{{AbstractOdeSystem}}} class has a member variable called the ''state variable vector'', which can
     * be used to hold the solution, and will be updated if a particular version of `Solve`
     * is called. This can be useful for embedding ODE models in a bigger system, since
     * the ODE models will then always contain their current solution state.
     */
    void TestOdeSolvingUsingStateVariable()
    {
        /* Define an instance of the ODE. See the class definition above.
         * Note that this ODE has a variable called {{{mStateVariables}}}, which has
         * been set to be a vector of size one, containing the value 1.0. */
        MyOdeUsingStateVariables my_ode_using_state_vars;

        /* To solve updating the state variable, just call the appropriate method on
         * a chosen solver. Note that no initial condition is required, no
         * {{{OdeSolution}}} is returned, and no sampling timestep is given. */
        EulerIvpOdeSolver euler_solver;
        euler_solver.SolveAndUpdateStateVariable(&my_ode_using_state_vars, 0.0, 1.0, 0.01);

        /* To see what the solution was at the end, we have to use the state variable. */
        std::cout << "Solution at end time = " << my_ode_using_state_vars.rGetStateVariables()[0] << "\n";
    }

    /*
     * === Solving n-dimensional ODEs ===
     *
     * Finally, here's a simple test showing how to solve a 2d ODE using the first method.
     * All that is different is the initial condition has be a vector of length 2, and returned
     * solution is of length 2 at every timestep.
     */
    void TestWith2dOde()
    {
        My2dOde my_2d_ode;
        EulerIvpOdeSolver euler_solver;

        /* Define the initial condition for each state variable. */
        std::vector<double> initial_condition;
        initial_condition.push_back(1.0);
        initial_condition.push_back(0.0);

        /* Solve, and print the solution as [time, y1, y2]. */
        OdeSolution solutions = euler_solver.Solve(&my_2d_ode, initial_condition, 0, 1, 0.01, 0.1);
        for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " "
                      << solutions.rGetSolutions()[i][0] << " "
                      << solutions.rGetSolutions()[i][1] << "\n";
        }
    }
};
#endif /*TESTSOLVINGODESTUTORIAL_HPP_*/
