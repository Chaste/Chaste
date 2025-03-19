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

#ifndef IMMERSEDBOUNDARYFFTINTERFACE_HPP_
#define IMMERSEDBOUNDARYFFTINTERFACE_HPP_

#include <complex>
#include <fftw3.h>
#include "ImmersedBoundaryMesh.hpp"
#include <pocketfft_hdronly.h>

/**
 * A class to interface with discrete Fourier transform libraries and perform
 * the necessary transforms for immersed boundary simulations.
 */
template<unsigned DIM>
class ImmersedBoundaryFftInterface
{
protected:

    /** The immersed boundary mesh. */
    ImmersedBoundaryMesh<DIM,DIM>* mpMesh;

    /** Pointer to the start of the input arrays. */
    double* mpInputArray;

    /** Pointer to the start of Fourier domain. */
    std::complex<double>* mpComplexArray;

    /** Pointer to the start output array. */
    double* mpOutputArray;

    /** Dimensions of the real arrays */
    pocketfft::shape_t mRealDims;

    /** Dimensions of the complex arrays */
    pocketfft::shape_t mCompDims;

    /** How many forward transforms to execute */
    unsigned mHowManyForward;

    /** How many inverse transofrms to execute */
    unsigned mHowManyInverse;

    /** Number of elements in each real array */
    int mRealSep;

    /** Number of elements in each complex array */
    int mCompSep;

    /** Distance between each element in the real arrays in bytes */
    long int mRealStride;

    /** Distance between each element in the complex arrays in bytes */
    long int mCompStride;

public:

    /**
     * Default constructor.
     *
     * @param pMesh the immersed boundary mesh
     * @param pIn pointer to the input array
     * @param pComplex pointer to the complex number array
     * @param pOut pointer to the output array
     * @param activeSources whether the population has active fluid sources
     */
    ImmersedBoundaryFftInterface(ImmersedBoundaryMesh<DIM,DIM>* pMesh,
                                 double* pIn,
                                 std::complex<double>* pComplex,
                                 double* pOut,
                                 bool activeSources);

    /**
     * Empty constructor.
     */
    ImmersedBoundaryFftInterface() = delete;

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundaryFftInterface();

    /** Performs forward fourier transforms */
    void FftExecuteForward();

    /** Performs inverse fourier transforms */
    void FftExecuteInverse();

    friend class TestImmersedBoundaryFftInterface;
};

#endif /*IMMERSEDBOUNDARYFFTINTERFACE_HPP_*/
