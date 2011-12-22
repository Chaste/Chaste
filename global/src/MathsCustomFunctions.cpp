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

bool CompareDoubles::IsNearZero(double number, double tolerance)
{
    return fabs(number) <= fabs(tolerance);
}

/**
 * Divide two non-negative floating point numbers, avoiding overflow and underflow.
 * @param number  the number to be divided
 * @param divisor  the number to divide by
 */
double SafeDivide(double number, double divisor)
{
    // Avoid overflow
    if (divisor < 1.0 && number > divisor*DBL_MAX)
    {
        return DBL_MAX;
    }

    // Avoid underflow
    if (number == 0.0 || (divisor > 1.0 && number < divisor*DBL_MIN))
    {
        return 0.0;
    }

    return number/divisor;

}

bool CompareDoubles::WithinRelativeTolerance(double number1, double number2, double tolerance)
{
    double diff = fabs(number1 - number2);
    double d1 = SafeDivide(diff, fabs(number1));
    double d2 = SafeDivide(diff, fabs(number2));

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
        double diff = fabs(number1 - number2);
        double d1 = SafeDivide(diff, fabs(number1));
        double d2 = SafeDivide(diff, fabs(number2));
        return d1 > d2 ? d1 : d2;
    }
}
