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

#include "NashHunterPoleZeroLaw.hpp"

template<>
NashHunterPoleZeroLaw<3>::NashHunterPoleZeroLaw()
        : PoleZeroMaterialLaw<3>()
{
    std::vector<std::vector<double> > k(3),a(3),b(3);
    for (unsigned i=0; i<3; i++)
    {
        k[i].resize(3);
        a[i].resize(3);
        b[i].resize(3);
    }

    /////////////////////////////////////////////////////////////////
    // Everything here has been entered in kPa.
    // All contraction models should return Ta in kPa.
    /////////////////////////////////////////////////////////////////
    k[0][0] = 2; //ff
    k[1][0] = k[0][1] = 1; //fs
    k[0][2] = k[2][0] = 1; //fn
    k[1][1] = 2; //ss
    k[1][2] = k[2][1] = 1; //sn
    k[2][2] = 2; //nn

    // dimensionless
    a[0][0] = 0.475; //ff
    a[1][0] = a[0][1] = 0.8; //fs
    a[2][0] = a[0][2] = 0.8; //fn
    a[1][1] = 0.619; //ss
    a[2][1] = a[1][2] = 0.8; //sn
    a[2][2] = 0.943; //nn

    // dimensionless
    b[0][0] = 1.5;
    b[1][0] = b[0][1] = 1.2;
    b[2][0] = b[0][2] = 1.2;
    b[1][1] = 1.5;
    b[2][1] = b[1][2] = 1.2;
    b[2][2] = 0.442;

    this->SetParameters(k,a,b);
}

template<>
NashHunterPoleZeroLaw<2>::NashHunterPoleZeroLaw()
        : PoleZeroMaterialLaw<2>()
{
    std::vector<std::vector<double> > k(2),a(2),b(2);
    for (unsigned i=0; i<2; i++)
    {
        k[i].resize(2);
        a[i].resize(2);
        b[i].resize(2);
    }

    k[0][0] = 2; //ff
    k[1][0] = k[0][1] = 1; //fs
    k[1][1] = 2; //ss

    a[0][0] = 0.475; //ff
    a[1][0] = a[0][1] = 0.8; //fs
    a[1][1] = 0.619; //ss

    b[0][0] = 1.5;
    b[1][0] = b[0][1] = 1.2;
    b[1][1] = 1.5;

    this->SetParameters(k,a,b);
}

// Explicit instantiation
template class NashHunterPoleZeroLaw<2>;
template class NashHunterPoleZeroLaw<3>;
