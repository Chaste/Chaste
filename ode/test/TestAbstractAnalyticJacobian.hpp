/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
