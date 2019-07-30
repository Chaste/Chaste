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

#include "MooneyRivlinMaterialLaw.hpp"

template<unsigned DIM>
double MooneyRivlinMaterialLaw<DIM>::Get_dW_dI1(double I1, double I2)
{
    return mC1;
}

template<unsigned DIM>
double MooneyRivlinMaterialLaw<DIM>::Get_dW_dI2(double I1, double I2)
{
    assert(DIM == 3); // LCOV_EXCL_LINE
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
    assert(DIM == 3); // LCOV_EXCL_LINE
    return 0.0;
}

template<unsigned DIM>
double MooneyRivlinMaterialLaw<DIM>::Get_d2W_dI1I2(double I1, double I2)
{
    assert(DIM == 3); // LCOV_EXCL_LINE
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
    assert(DIM == 3); // LCOV_EXCL_LINE
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

// Explicit instantiation
template class MooneyRivlinMaterialLaw<2>;
template class MooneyRivlinMaterialLaw<3>;
