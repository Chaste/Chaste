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

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

template class NashHunterPoleZeroLaw<2>;
template class NashHunterPoleZeroLaw<3>;
