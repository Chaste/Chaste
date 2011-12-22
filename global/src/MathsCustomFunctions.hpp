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

#ifndef MATHSCUSTOMFUNCTIONS_HPP_
#define MATHSCUSTOMFUNCTIONS_HPP_

#include <cfloat>

/**
 * @file
 * This file contains some utility functions and a small class for dealing with floating
 * point numbers.
 */

/**
 * Replacement "pow" function.
 *
 * @param x number to be raised to a small power
 * @param exponent small integer exponent
 * @return x^exponent a.k.a x**exponent.
 */
double SmallPow(double x, unsigned exponent);

/**
 * Uses fmod to determine if smallerNumber divides the largerNumber.
 * We expect smallerNumber/largerNumber <= 1 and therefore
 * fmod(largerNumber,smallerNumber) should be close to zero or close to smallerNumber.
 *
 * @param smallerNumber the smaller
 * @param largerNumber the larger
 */
bool Divides(double smallerNumber, double largerNumber);

/**
 * Utility static methods for comparing floating point numbers, based on
 * boost/test/floating_point_comparison.hpp.
 */
class CompareDoubles
{
public:
    /**
     * Test whether a number is close to zero, as defined by the given tolerance.
     * A number equal in magnitude to the tolerance is OK.
     *
     * @param number  the number to compare
     * @param tolerance  how close it must be
     */
    static bool IsNearZero(double number, double tolerance);

    /**
     * Test whether two numbers are close within the given relative tolerance.
     *
     * @param number1  the first number to compare
     * @param number2  the second number to compare
     * @param tolerance  the relative tolerance to use
     */
    static bool WithinRelativeTolerance(double number1, double number2, double tolerance);

    /**
     * Test whether two numbers are close within the given absolute tolerance.
     * A difference of exactly the tolerance is OK.
     *
     * @param number1  the first number to compare
     * @param number2  the second number to compare
     * @param tolerance  the absolute tolerance to use
     */
    static bool WithinAbsoluteTolerance(double number1, double number2, double tolerance);

    /**
     * Test whether two numbers are close within the given tolerance, and print an
     * error message to stdout if not.
     *
     * @param number1  the first number to compare
     * @param number2  the second number to compare
     * @param tolerance  the tolerance to use
     * @param toleranceIsAbsolute  whether the tolerance is absolute (true) or relative (false)
     */
    static bool WithinTolerance(double number1, double number2, double tolerance, bool toleranceIsAbsolute);

    /**
     * Test whether two numbers are close within the given tolerances.  If either the relative or
     * absolute tolerance is satisfied, then the result is true.
     *
     * @param number1  the first number to compare
     * @param number2  the second number to compare
     * @param relTol  the relative tolerance to compare under
     * @param absTol  the absolute tolerance to compare under
     * @param printError  whether to print an error message to stdout
     */
    static bool WithinAnyTolerance(double number1, double number2,
                                   double relTol=DBL_EPSILON, double absTol=DBL_EPSILON,
                                   bool printError=false);

    /**
     * Compute the relative or absolute difference between two numbers.
     *
     * @param number1  the first number to compare
     * @param number2  the second number to compare
     * @param toleranceIsAbsolute  whether the tolerance is absolute (true) or relative (false)
     */
    static double Difference(double number1, double number2, bool toleranceIsAbsolute);
};

#endif /*MATHSCUSTOMFUNCTIONS_HPP_*/
