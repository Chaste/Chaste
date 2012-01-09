/*

Copyright (C) University of Oxford, 2005-2012

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

#include "MooneyRivlinMaterialLaw.hpp"

template<unsigned DIM>
double MooneyRivlinMaterialLaw<DIM>::Get_dW_dI1(double I1, double I2)
{
    return mC1;
}

template<unsigned DIM>
double MooneyRivlinMaterialLaw<DIM>::Get_dW_dI2(double I1, double I2)
{
    /*
     * This is covered, but gcov doesn't see this as being covered
     * for some reason, maybe because of optimisations.
     */
    #define COVERAGE_IGNORE
    assert(DIM == 3);
    #undef COVERAGE_IGNORE
    return mC2;
}

template<unsigned DIM>
double MooneyRivlinMaterialLaw<DIM>::Get_d2W_dI1(double I1, double I2)
{
    return 0.0;
}

template<unsigned DIM>
double MooneyRivlinMaterialLaw<DIM>::Get_d2W_dI2(double I1, double I2)
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
double MooneyRivlinMaterialLaw<DIM>::Get_d2W_dI1I2(double I1, double I2)
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
double MooneyRivlinMaterialLaw<DIM>::GetC1()
{
    return mC1;
}

template<unsigned DIM>
double MooneyRivlinMaterialLaw<DIM>::GetC2()
{
    assert(DIM==3);
    return mC2;
}

template<unsigned DIM>
MooneyRivlinMaterialLaw<DIM>::MooneyRivlinMaterialLaw(double c1, double c2)
    : mC1(c1),
      mC2(c2)
{
    assert(DIM==2 || DIM ==3);

    // If dim==3, check that c2 was passed in, ie c2 isn't the default value
    if ((DIM==3) && (c2<MINUS_LARGE+1))
    {
        EXCEPTION("Two parameters needed for 3d Mooney-Rivlin");
    }

    if (c1 < 0.0)
    {
        EXCEPTION("c1 must be positive in mooney-rivlin"); // is this correct?
    }
}

template<unsigned DIM>
void MooneyRivlinMaterialLaw<DIM>::ScaleMaterialParameters(double scaleFactor)
{
    assert(scaleFactor > 0.0);
    mC1 /= scaleFactor;
    mC2 /= scaleFactor;
}

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

template class MooneyRivlinMaterialLaw<2>;
template class MooneyRivlinMaterialLaw<3>;
