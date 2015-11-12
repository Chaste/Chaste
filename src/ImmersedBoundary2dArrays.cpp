
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

#include "Debug.hpp"

template<unsigned DIM>
ImmersedBoundary2dArrays<DIM>::ImmersedBoundary2dArrays(ImmersedBoundaryMesh<DIM,DIM>* p_mesh, double dt, double reynoldsNumber)
        : mpMesh(p_mesh)
{
    unsigned num_gridpts_x = mpMesh->GetNumGridPtsX();
    unsigned num_gridpts_y = mpMesh->GetNumGridPtsY();

    // We require an even number of grid points
    assert(num_gridpts_y % 2 == 0);

    /*
     * Resize the grids.  All complex grids are half-sized in the y-coordinate due to redundancy inherent in the
     * fast-Fourier method for solving Navier-Stokes.
     */
    unsigned reduced_y = 1 + (num_gridpts_y/2);

    // Resize real arrays to X by Y
    mForceGrids.resize(extents[2][num_gridpts_x][num_gridpts_y]);
    mRightHandSideGrids.resize(extents[2][num_gridpts_x][num_gridpts_y]);

    // Resize Fourier-domain arrays
    mOperator1.resize(extents[num_gridpts_x][reduced_y]);
    mOperator2.resize(extents[num_gridpts_x][reduced_y]);
    mFourierGrids.resize(extents[2][num_gridpts_x][reduced_y]);
    mPressureGrid.resize(extents[num_gridpts_x][reduced_y]);

    mSin2x.resize(num_gridpts_x);
    mSin2y.resize(reduced_y);

    /*
     * There are several constants used in the Fourier-domain as part of the Navier-Stokes solution which are constant
     * once grid sizes are known.  We pre-calculate these to eliminate re-calculation at every timestep.
     */
    double x_spacing = 1.0 / (double) num_gridpts_x;
    double y_spacing = 1.0 / (double) num_gridpts_y;

    for (unsigned x = 0 ; x < num_gridpts_x ; x++)
    {
        mSin2x[x] = sin(2 * M_PI * (double) x * x_spacing);
    }

    for (unsigned y = 0 ; y < reduced_y ; y++)
    {
        mSin2y[y] = sin(2 * M_PI * (double) y * y_spacing);
    }

    for (unsigned x = 0 ; x < num_gridpts_x ; x++)
    {
        for (unsigned y = 0 ; y < reduced_y ; y++)
        {
            mOperator1[x][y] = (mSin2x[x] * mSin2x[x] / (x_spacing * x_spacing)) + (mSin2y[y] * mSin2y[y] / (y_spacing * y_spacing));
            mOperator1[x][y] *= dt / reynoldsNumber;

            double sin_x = sin(M_PI * (double) x * x_spacing);
            double sin_y = sin(M_PI * (double) y * y_spacing);

            mOperator2[x][y] = (sin_x * sin_x / (x_spacing * x_spacing)) + (sin_y * sin_y / (y_spacing * y_spacing));
            mOperator2[x][y] *= 4.0 * dt / reynoldsNumber;
            mOperator2[x][y] += 1.0;
        }
    }

    /*
     * Finally, plan the discrete Fourier transforms:
     *
     *  * Forward are real-to-complex and out-of-place
     *  * Backward are complex-to-real and out-of-place
     *
     * Because the above arrays are created once and stay in the same place in memory throughout the simulation, we may
     * plan the transforms with references to locations in these arrays, and execute the plans as necessary each
     * timestep of the simulation.
     */

    // We first set the pointers to arrays where the data will be stored
    double* p_in = &(mRightHandSideGrids[0][0][0]);
    fftw_complex* p_complex = reinterpret_cast<fftw_complex*>(&(mFourierGrids[0][0][0]));
    double* p_out = &(mpMesh->rGetModifiable2dVelocityGrids()[0][0][0]);

    // Plan variables
    int rank = 2;                                       // Number of dimensions for each array
    int real_dims[] = {num_gridpts_x, num_gridpts_y};   // Dimensions of each real array
    int comp_dims[] = {num_gridpts_x, reduced_y};       // Dimensions of each complex array
    int how_many = 2;                                   // Number of transforms
    int real_sep = num_gridpts_x * num_gridpts_y;       // How many doubles between start of first array and start of second
    int comp_sep = num_gridpts_x * reduced_y;           // How many fftw_complex between start of first array and start of second
    int real_stride = 1;                                // Each real array is contiguous in memory
    int comp_stride = 1;                                // Each complex array is contiguous in memory
    int* real_nembed = real_dims;
    int* comp_nembed = comp_dims;

    mFftwForwardPlan = fftw_plan_many_dft_r2c(rank, real_dims, how_many,
                                              p_in,      real_nembed, real_stride, real_sep,
                                              p_complex, comp_nembed, comp_stride, comp_sep,
                                              FFTW_PATIENT);

    mFftwInversePlan = fftw_plan_many_dft_c2r(rank, real_dims, how_many,
                                              p_complex, comp_nembed, comp_stride, comp_sep,
                                              p_out,     real_nembed, real_stride, real_sep,
                                              FFTW_PATIENT);
}

template<unsigned DIM>
ImmersedBoundary2dArrays<DIM>::~ImmersedBoundary2dArrays()
{
    fftw_destroy_plan(mFftwForwardPlan);
    fftw_destroy_plan(mFftwInversePlan);
}

template<unsigned DIM>
multi_array<double, 3>& ImmersedBoundary2dArrays<DIM>::rGetModifiableForceGrids()
{
    return mForceGrids;
}

template<unsigned DIM>
multi_array<double, 3>& ImmersedBoundary2dArrays<DIM>::rGetModifiableRightHandSideGrids()
{
    return mRightHandSideGrids;
}

template<unsigned DIM>
multi_array<std::complex<double>, 3>& ImmersedBoundary2dArrays<DIM>::rGetModifiableFourierGrids()
{
    return mFourierGrids;
}

template<unsigned DIM>
multi_array<std::complex<double>, 2>& ImmersedBoundary2dArrays<DIM>::rGetModifiablePressureGrid()
{
    return mPressureGrid;
}

template<unsigned DIM>
const multi_array<double, 2>& ImmersedBoundary2dArrays<DIM>::rGetOperator1() const
{
    return mOperator1;
}

template<unsigned DIM>
const multi_array<double, 2>& ImmersedBoundary2dArrays<DIM>::rGetOperator2() const
{
    return mOperator2;
}

template<unsigned DIM>
const std::vector<double>& ImmersedBoundary2dArrays<DIM>::rGetSin2x() const
{
    return mSin2x;
}

template<unsigned DIM>
const std::vector<double>& ImmersedBoundary2dArrays<DIM>::rGetSin2y() const
{
    return mSin2y;
}

template<unsigned DIM>
void ImmersedBoundary2dArrays<DIM>::FftwExecuteForward()
{
    fftw_execute(mFftwForwardPlan);
}

template<unsigned DIM>
void ImmersedBoundary2dArrays<DIM>::FftwExecuteInverse()
{
    fftw_execute(mFftwInversePlan);
}

// Explicit instantiation
template class ImmersedBoundary2dArrays<1>;
template class ImmersedBoundary2dArrays<2>;
template class ImmersedBoundary2dArrays<3>;