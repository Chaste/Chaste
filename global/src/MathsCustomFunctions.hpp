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
 * Replacement "pow" function for unsigned integer values
 *
 * @param x number to be raised to a small power
 * @param exponent small integer exponent
 * @return x^exponent a.k.a x**exponent.
 */
unsigned SmallPow(unsigned x, unsigned exponent);

/**
 * @return true if the smaller number divides the larger (to machine precision)
 * Uses fmod to determine if smallerNumber divides the largerNumber.
 * We expect smallerNumber/largerNumber <= 1 and therefore
 * fmod(largerNumber,smallerNumber) should be close to zero or close to smallerNumber.
 *
 * @param smallerNumber the smaller
 * @param largerNumber the larger
 */
bool Divides(double smallerNumber, double largerNumber);

/**
 * @return the result of dividing (unsigned) numerator by denominator,
 * rounded up (away from 0), and overflow-safe for large integers up to UINT_MAX.
 *
 * @param numerator the numerator
 * @param denominator the denominator
 */
unsigned CeilDivide(unsigned numerator, unsigned denominator);

/**
 * @return the sign of the argument, i.e. -1 if value<0, 0 if value=0, or +1 if value>0.
 *
 * @param value  the argument value
 */
double Signum(double value);

/**
 * @return the results of dividing two non-negative floating point numbers, avoiding overflow and underflow.
 *
 * @param numerator  the number to be divided
 * @param divisor  the number to divide by
 */
double SafeDivide(double numerator, double divisor);

/**
 * Utility static methods for comparing floating point numbers, based on
 * boost/test/floating_point_comparison.hpp.
 */
class CompareDoubles
{
public:
    /**
     * @return true if a number is close to zero, as defined by the given tolerance.
     * A number equal in magnitude to the tolerance is OK.
     *
     * @param number  the number to compare
     * @param tolerance  how close it must be
     */
    static bool IsNearZero(double number, double tolerance);

    /**
     * @return true if  two numbers are close within the given relative tolerance.
     *
     * @param number1  the first number to compare
     * @param number2  the second number to compare
     * @param tolerance  the relative tolerance to use
     */
    static bool WithinRelativeTolerance(double number1, double number2, double tolerance);

    /**
     * @return true if two numbers are close within the given absolute tolerance.
     * A difference of exactly the tolerance is OK.
     *
     * @param number1  the first number to compare
     * @param number2  the second number to compare
     * @param tolerance  the absolute tolerance to use
     */
    static bool WithinAbsoluteTolerance(double number1, double number2, double tolerance);

    /**
     * @return true if two numbers are close within the given tolerance, and print an
     * error message to stdout if not.
     *
     * @param number1  the first number to compare
     * @param number2  the second number to compare
     * @param tolerance  the tolerance to use
     * @param toleranceIsAbsolute  whether the tolerance is absolute (true) or relative (false)
     */
    static bool WithinTolerance(double number1, double number2, double tolerance, bool toleranceIsAbsolute);

    /**
     * @return true if two numbers are close within the given tolerances.  If either the relative or
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
     * @return the relative or absolute difference between two numbers.
     *
     * @param number1  the first number to compare
     * @param number2  the second number to compare
     * @param toleranceIsAbsolute  whether the tolerance is absolute (true) or relative (false)
     */
    static double Difference(double number1, double number2, bool toleranceIsAbsolute);
};

#endif /*MATHSCUSTOMFUNCTIONS_HPP_*/
