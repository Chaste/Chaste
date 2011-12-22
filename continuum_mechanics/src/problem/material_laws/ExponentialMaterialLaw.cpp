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

#include "ExponentialMaterialLaw.hpp"

template<unsigned DIM>
ExponentialMaterialLaw<DIM>::ExponentialMaterialLaw(double a, double b)
    : mA(a),
      mB(b)
{
    assert(DIM==2 || DIM==3);
    if (a < 0.0)
    {
        EXCEPTION("a must be positive");
    }
}

template<unsigned DIM>
double ExponentialMaterialLaw<DIM>::GetA()
{
    return mA;
}

template<unsigned DIM>
double ExponentialMaterialLaw<DIM>::GetB()
{
    return mB;
}

template<unsigned DIM>
double ExponentialMaterialLaw<DIM>::Get_dW_dI1(double I1, double I2)
{
    return mA * mB * exp(mB*(I1-DIM));
}

template<unsigned DIM>
double ExponentialMaterialLaw<DIM>::Get_dW_dI2(double I1, double I2)
{
    /*
     * This is covered, but gcov doesn't see this as being covered
     * for some reason, maybe because of optimisations.
     */
    #define COVERAGE_IGNORE
    assert(DIM == 3);
    #undef COVERAGE_IGNORE

    return 0.0;
}

template<unsigned DIM>
double ExponentialMaterialLaw<DIM>::Get_d2W_dI1(double I1, double I2)
{
    return mA * mB * mB * exp(mB*(I1-DIM));
}

template<unsigned DIM>
double ExponentialMaterialLaw<DIM>::Get_d2W_dI2(double I1, double I2)
{
    /*
     * This is covered, but gcov doesn't see this as being covered
     * for some reason, maybe because of optimisations.
     */
    #define COVERAGE_IGNORE
    assert(DIM == 3);
    #undef COVERAGE_IGNORE

    return 0.0;
}

template<unsigned DIM>
double ExponentialMaterialLaw<DIM>::Get_d2W_dI1I2(double I1, double I2)
{
    /*
     * This is covered, but gcov doesn't see this as being covered
     * for some reason, maybe because of optimisations.
     */
    #define COVERAGE_IGNORE
    assert(DIM == 3);
    #undef COVERAGE_IGNORE

    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

template class ExponentialMaterialLaw<2>;
template class ExponentialMaterialLaw<3>;
