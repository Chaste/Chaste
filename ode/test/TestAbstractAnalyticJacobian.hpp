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


#ifndef _TESTABSTRACTANALYTICJACOBIAN_HPP_
#define _TESTABSTRACTANALYTICJACOBIAN_HPP_

// TestAbstractAnalyticJacobian.hpp

#include <cmath>
#include <iostream>
#include <vector>
#include "OdeWithJacobian1.hpp"
#include "OdeWithJacobian2.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

// Tolerance for tests
double tol = 0.01;


class TestAbstractAnalyticJacobian : public CxxTest::TestSuite
{
public:

    void TestJacobianOne()
    {
        // Pointer to TestOde1 class
        OdeWithJacobian1 ode_system;

        std::vector<double> solution_guess(1);
        solution_guess[0] = 2.0;

        // Set up a Jacobian matrix for function to put values in
        double** jacobian;

        jacobian = new double*[1];
        jacobian[0] = new double[1];

        // This is the function we are testing...
        ode_system.AnalyticJacobian(solution_guess, jacobian, 1.0, 0.01);

        TS_ASSERT_DELTA(jacobian[0][0], 0.96, tol);

        delete[] jacobian[0];
        delete[] jacobian;
    }

    void TestJacobianTwo()
    {
        // Pointer to TestOde1 class
        OdeWithJacobian2 ode_system;

        std::vector<double>  solution_guess(2);
        solution_guess[0] = 1.0;
        solution_guess[1] = 2.0;

        // Set up a Jacobian matrix for function to put values in
        double** jacobian;

        jacobian = new double* [2];
        jacobian[0] = new double[2];
        jacobian[1] = new double[2];

        // This is the function we are testing...
        ode_system.AnalyticJacobian(solution_guess, jacobian, 1.0, 0.01);

        TS_ASSERT_DELTA(jacobian[0][0], 0.98, tol);
        TS_ASSERT_DELTA(jacobian[0][1], -0.04, tol);
        TS_ASSERT_DELTA(jacobian[1][0], -0.02, tol);
        TS_ASSERT_DELTA(jacobian[1][1], 0.92, tol);

        delete[] jacobian[0];
        delete[] jacobian[1];
        delete[] jacobian;
    }
};

#endif //_TESTABSTRACTANALYTICJACOBIAN_HPP_
