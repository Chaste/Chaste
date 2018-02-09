/*

Copyright (c) 2005-2018, University of Oxford.
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
      mpComplexArray(reinterpret_cast<fftw_complex*>(pComplex)),
      mpOutputArray(pOut)
{
    /*
     * Set up fftw routines
     */

    // Forget all wisdom; the correct wisdom will be loaded from file
    void fftw_forget_wisdom();

    // Load wisdom from file
    std::string wisdom_filename = "fftw.wisdom";
    FileFinder file_finder(wisdom_filename, RelativeTo::ChasteTestOutput);

    std::string wisdom_path = file_finder.GetAbsolutePath();
    int wisdom_flag;

    if (file_finder.IsFile())
    {
        wisdom_flag = fftw_import_wisdom_from_filename(wisdom_path.c_str());

        // 1 means success, 0 indicates a failure
        if (wisdom_flag != 1)
        {
            WARNING("fftw wisdom not imported correctly from " + wisdom_path);
        }
    }
    else // file in test output folder not found
    {
        WARNING("Cannot find fftw wisdom file at " + wisdom_path +
                ". It is strongly recommended to run TestGenerateFftwWisdom.hpp.");

        wisdom_flag = fftw_import_system_wisdom();

        // 1 means success, 0 indicates a failure
        if (wisdom_flag != 1)
        {
            WARNING("fftw system wisdom not imported correctly");
        }
    }

    int num_gridpts_x = (int)mpMesh->GetNumGridPtsX();
    int num_gridpts_y = (int)mpMesh->GetNumGridPtsY();

    // We require an even number of grid points
    assert(num_gridpts_y % 2 == 0);

    /*
     * Resize the grids.  All complex grids are half-sized in the y-coordinate due to redundancy inherent in the
     * fast-Fourier method for solving Navier-Stokes.
     */
    int reduced_y = 1 + (num_gridpts_y/2);

    /*
     * Plan the discrete Fourier transforms:
     *
     *  * Forward are real-to-complex and out-of-place
     *  * Backward are complex-to-real and out-of-place
     *
     * Because the above arrays are created once and stay in the same place in memory throughout the simulation, we may
     * plan the transforms with references to locations in these arrays, and execute the plans as necessary each
     * timestep of the simulation.
     */

    // Plan variables
    int rank = 2;                                       // Number of dimensions for each array
    int real_dims[] = {num_gridpts_x, num_gridpts_y};   // Dimensions of each real array
    int comp_dims[] = {num_gridpts_x, reduced_y};       // Dimensions of each complex array
    int how_many_forward = 2 + (int)activeSources;      // Number of forward transforms (one more if sources are active)
    int how_many_inverse = 2;                           // Number of inverse transforms (always 2)
    int real_sep = num_gridpts_x * num_gridpts_y;       // How many doubles between start of first array and start of second
    int comp_sep = num_gridpts_x * reduced_y;           // How many fftw_complex between start of first array and start of second
    int real_stride = 1;                                // Each real array is contiguous in memory
    int comp_stride = 1;                                // Each complex array is contiguous in memory
    int* real_nembed = real_dims;
    int* comp_nembed = comp_dims;

    mFftwForwardPlan = fftw_plan_many_dft_r2c(rank, real_dims, how_many_forward,
                                              mpInputArray,   real_nembed, real_stride, real_sep,
                                              mpComplexArray, comp_nembed, comp_stride, comp_sep,
                                              FFTW_PATIENT);

    mFftwInversePlan = fftw_plan_many_dft_c2r(rank, real_dims, how_many_inverse,
                                              mpComplexArray, comp_nembed, comp_stride, comp_sep,
                                              mpOutputArray,  real_nembed, real_stride, real_sep,
                                              FFTW_PATIENT);
}

template<unsigned DIM>
ImmersedBoundaryFftInterface<DIM>::~ImmersedBoundaryFftInterface()
{
    fftw_destroy_plan(mFftwForwardPlan);
    fftw_destroy_plan(mFftwInversePlan);
}

template<unsigned DIM>
void ImmersedBoundaryFftInterface<DIM>::FftExecuteForward()
{
    fftw_execute(mFftwForwardPlan);
}

template<unsigned DIM>
void ImmersedBoundaryFftInterface<DIM>::FftExecuteInverse()
{
    fftw_execute(mFftwInversePlan);
}

// Explicit instantiation
template class ImmersedBoundaryFftInterface<1>;
template class ImmersedBoundaryFftInterface<2>;
template class ImmersedBoundaryFftInterface<3>;
