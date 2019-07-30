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

#include "MathsCustomFunctions.hpp"

#include <cmath>
#include <iostream>

double SmallPow(double x, unsigned exponent)
{
    switch (exponent)
    {
        case 0:
        {
            return 1.0;
        }
        case 1:
        {
            return x;
        }
        case 2:
        {
            return x*x;
        }
        case 3:
        {
            return x*x*x;
        }
        default:
        {
            if (exponent % 2 == 0)
            {
                // Even power
                double partial_answer = SmallPow(x, exponent/2);
                return partial_answer*partial_answer;
            }
            else
            {
                // Odd power
                return SmallPow(x, exponent-1)*x;
            }
        }
    }
}
unsigned SmallPow(unsigned x, unsigned exponent)
{
    switch (exponent)
    {
        case 0:
        {
            return 1u;
        }
        case 1:
        {
            return x;
        }
        case 2:
        {
            return x*x;
        }
        case 3:
        {
            return x*x*x;
        }
        default:
        {
            if (exponent % 2 == 0)
            {
                // Even power
                unsigned partial_answer = SmallPow(x, exponent/2);
                return partial_answer*partial_answer;
            }
            else
            {
                // Odd power
                return SmallPow(x, exponent-1)*x;
            }
        }
    }
}

bool Divides(double smallerNumber, double largerNumber)
{
    double remainder = fmod(largerNumber, smallerNumber);
    /*
     * Is the remainder close to zero? Note that the comparison is scaled
     * with respect to the larger of the numbers.
     */
    if (remainder < DBL_EPSILON*largerNumber)
    {
        return true;
    }
    /*
     * Is the remainder close to smallerNumber? Note that the comparison
     * is scaled with respect to the larger of the numbers.
     */
    if (fabs(remainder-smallerNumber) < DBL_EPSILON*largerNumber)
    {
        return true;
    }

    return false;
}

unsigned CeilDivide(unsigned numerator, unsigned denominator)
{
    if (numerator == 0u)
    {
        return 0u;
    }
    else
    {
        // Overflow-safe for large numbers, but not valid for numerator==0.
        return ((numerator - 1u) / denominator) + 1u;
    }
}

double Signum(double value)
{
    return (0.0 < value) - (value < 0.0);
}

bool CompareDoubles::IsNearZero(double number, double tolerance)
{
    return fabs(number) <= fabs(tolerance);
}

double SafeDivide(double numerator, double divisor)
{
    // Avoid overflow
    if (divisor < 1.0 && numerator > divisor*DBL_MAX)
    {
        return DBL_MAX;
    }

    // Avoid underflow
    if (numerator == 0.0 || (divisor > 1.0 && numerator < divisor*DBL_MIN))
    {
        return 0.0;
    }

    return numerator/divisor;
}

bool CompareDoubles::WithinRelativeTolerance(double number1, double number2, double tolerance)
{
    double difference = fabs(number1 - number2);
    double d1 = SafeDivide(difference, fabs(number1));
    double d2 = SafeDivide(difference, fabs(number2));

    return d1 <= tolerance && d2 <= tolerance;
}

bool CompareDoubles::WithinAbsoluteTolerance(double number1, double number2, double tolerance)
{
    return fabs(number1 - number2) <= tolerance;
}

bool CompareDoubles::WithinAnyTolerance(double number1, double number2, double relTol, double absTol, bool printError)
{
    bool ok = WithinAbsoluteTolerance(number1, number2, absTol) || WithinRelativeTolerance(number1, number2, relTol);
    if (printError && !ok)
    {
        std::cout << "CompareDoubles::WithinAnyTolerance: " << number1 << " and " << number2
                  << " differ by more than relative tolerance of " << relTol
                  << " and absolute tolerance of " << absTol << std::endl;
    }
    return ok;
}

bool CompareDoubles::WithinTolerance(double number1, double number2, double tolerance, bool toleranceIsAbsolute)
{
    bool ok;
    if (toleranceIsAbsolute)
    {
        ok = WithinAbsoluteTolerance(number1, number2, tolerance);
    }
    else
    {
        ok = WithinRelativeTolerance(number1, number2, tolerance);
    }
    if (!ok)
    {
        std::cout << "CompareDoubles::WithinTolerance: " << number1 << " and " << number2
                  << " differ by more than " << (toleranceIsAbsolute ? "absolute" : "relative")
                  << " tolerance of " << tolerance << std::endl;
    }
    return ok;
}

double CompareDoubles::Difference(double number1, double number2, bool toleranceIsAbsolute)
{
    if (toleranceIsAbsolute)
    {
        return fabs(number1 - number2);
    }
    else
    {
        double difference = fabs(number1 - number2);
        double d1 = SafeDivide(difference, fabs(number1));
        double d2 = SafeDivide(difference, fabs(number2));
        return d1 > d2 ? d1 : d2;
    }
}
