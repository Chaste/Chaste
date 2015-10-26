
/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "ImmersedBoundary2dArrays.hpp"
#include <assert.h>

ImmersedBoundary2dArrays::ImmersedBoundary2dArrays(unsigned numGridPtsX, unsigned numGridPtsY)
{
    // The complex grids are reduced in size due to redundancy in the fourier domain
    unsigned reduced_y = 1 + (numGridPtsY/2);

    // We require an even number of grid points
    assert(numGridPtsY % 2 == 0);

    // 2 X by Y grids, one for each dimension
    mForceGrids.resize(extents[2][numGridPtsX][numGridPtsY]);

    // 2 X by Y grids, one for each dimension
    mRightHandSideGrids.resize(extents[2][numGridPtsX][numGridPtsY]);

    // 2 X by Y grids, one for each dimension, but some redundancy due to FFT of real data
    mFourierGrids.resize(extents[2][numGridPtsX][reduced_y]);

    // 1 X by Y grid
    mPressureGrid.resize(extents[numGridPtsX][reduced_y]);
}

ImmersedBoundary2dArrays::~ImmersedBoundary2dArrays()
{
}

multi_array<double, 3>& ImmersedBoundary2dArrays::rGetModifiableForceGrids()
{
    return mForceGrids;
}

multi_array<double, 3>& ImmersedBoundary2dArrays::rGetModifiableRightHandSideGrids()
{
    return mRightHandSideGrids;
}

multi_array<std::complex<double>, 3>& ImmersedBoundary2dArrays::rGetModifiableFourierGrids()
{
    return mFourierGrids;
}

multi_array<std::complex<double>, 2>& ImmersedBoundary2dArrays::rGetModifiablePressureGrid()
{
    return mPressureGrid;
}