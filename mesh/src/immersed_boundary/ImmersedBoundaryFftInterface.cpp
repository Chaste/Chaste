/*

Copyright (c) 2005-2025, University of Oxford.
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

#include "ImmersedBoundaryFftInterface.hpp"
#include <assert.h>
#include "FileFinder.hpp"
#include "Warnings.hpp"

template<unsigned DIM>
ImmersedBoundaryFftInterface<DIM>::ImmersedBoundaryFftInterface(ImmersedBoundaryMesh<DIM,DIM>* pMesh,
                                                                double* pIn,
                                                                std::complex<double>* pComplex,
                                                                double* pOut,
                                                                bool activeSources)
    : mpMesh(pMesh),
      mpInputArray(pIn),
      mpComplexArray(pComplex),
      mpOutputArray(pOut)
{
    int num_grid_pts_x = (int)mpMesh->GetNumGridPtsX();
    int num_grid_pts_y = (int)mpMesh->GetNumGridPtsY();

    // We require an even number of grid points
    assert(num_grid_pts_y % 2 == 0);

    /*
     * Resize the grids.  All complex grids are half-sized in the y-coordinate due to redundancy inherent in the
     * fast-Fourier method for solving Navier-Stokes.
     */
    int reduced_y = 1 + (num_grid_pts_y/2);

    // Plan variables
    mRealDims = {(long unsigned int)num_grid_pts_x, (long unsigned int)num_grid_pts_y};   // Dimensions of each real array
    mCompDims = {(long unsigned int)num_grid_pts_x, (long unsigned int)reduced_y};       // Dimensions of each complex array

    mHowManyForward = 2 + (unsigned)activeSources;      // Number of forward transforms (one more if sources are active)
    mHowManyInverse = 2;                           // Number of inverse transforms (always 2)

    mRealSep = num_grid_pts_x * num_grid_pts_y;       // How many doubles between start of first array and start of second
    mCompSep = num_grid_pts_x * reduced_y;           // How many fftw_complex between start of first array and start of second

    mRealStride = sizeof(double);                                // Each real array is contiguous in memory
    mCompStride = sizeof(std::complex<double>);                                // Each complex array is contiguous in memory

/*

 multi-dimensional arrays are stored row-major
 real_dims = dimensions of matrix
 real_sep = pointer offset between matrices
 how_many_forward = number of matrices to operate on
 stride = 1 => matrices are continguously stored one after the other
 */
}

template<unsigned DIM>
ImmersedBoundaryFftInterface<DIM>::~ImmersedBoundaryFftInterface()
{
}

template<unsigned DIM>
void ImmersedBoundaryFftInterface<DIM>::FftExecuteForward()
{
    // Real to complex
    pocketfft::stride_t rStride = {mRealStride*static_cast<long int>(mRealDims[1]), mRealStride};
    pocketfft::stride_t cStride = {mCompStride*static_cast<long int>(mCompDims[1]), mCompStride};
    pocketfft::shape_t axes = {0, 1};

    for (unsigned i = 0; i < mHowManyForward; i++)
    {
        pocketfft::r2c<double>(mRealDims, rStride, cStride, axes, true, mpInputArray + i*mRealSep, mpComplexArray + i*mCompSep, 1.0, 1);
    }
}

template<unsigned DIM>
void ImmersedBoundaryFftInterface<DIM>::FftExecuteInverse()
{
    // Complex to real
    pocketfft::stride_t rStride = {mRealStride*static_cast<long int>(mRealDims[1]), mRealStride};
    pocketfft::stride_t cStride = {mCompStride*static_cast<long int>(mCompDims[1]), mCompStride};
    pocketfft::shape_t axes = {0, 1};

    for (unsigned i = 0; i < mHowManyInverse; i++)
    {
        // For c2r, output array dims are supplied
        pocketfft::c2r(mRealDims, cStride, rStride, axes, false, mpComplexArray + i*mCompSep, mpOutputArray + i*mRealSep, 1.0, 1);
    }
}

// Explicit instantiation
template class ImmersedBoundaryFftInterface<1>;
template class ImmersedBoundaryFftInterface<2>;
template class ImmersedBoundaryFftInterface<3>;
