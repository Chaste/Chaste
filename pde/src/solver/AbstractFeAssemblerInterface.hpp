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

#ifndef ABSTRACTFEASSEMBLERINTERFACE_HPP_
#define ABSTRACTFEASSEMBLERINTERFACE_HPP_

#include <cassert>
#include "UblasCustomFunctions.hpp"
#include "PetscTools.hpp"

/**
 *   A common bass class for AbstractFeVolumeIntegralAssembler (the main abstract assembler class), and other assembler classes
 *   (including continuum mechanics assemblers, which is why this class is separate to AbstractFeAssemblerInterface).
 *
 *   See AbstractFeVolumeIntegralAssembler documentation for info on these assembler classes.
 */

template <bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX>
class AbstractFeAssemblerInterface  : private boost::noncopyable
{
protected:
    /** The vector to be assembled (only used if CAN_ASSEMBLE_VECTOR == true). */
    Vec mVectorToAssemble;

    /** The matrix to be assembled (only used if CAN_ASSEMBLE_MATRIX == true). */
    Mat mMatrixToAssemble;

    /**
     * Whether to assemble the matrix (an assembler may be able to assemble matrices
     * (CAN_ASSEMBLE_MATRIX==true), but may not want to do so each timestep, hence
     * this second boolean.
     */
    bool mAssembleMatrix;

    /** Whether to assemble the vector. */
    bool mAssembleVector;

    /** Whether to zero the given matrix before assembly, or just add to it. */
    bool mZeroMatrixBeforeAssembly;

    /** Whether to zero the given vector before assembly, or just add to it. */
    bool mZeroVectorBeforeAssembly;

    /** Ownership range of the vector/matrix - lowest component owned. */
    PetscInt mOwnershipRangeLo;

    /** Ownership range of the vector/matrix - highest component owned +1. */
    PetscInt mOwnershipRangeHi;

    /**
     * The main assembly method. Protected, should only be called through Assemble(),
     * AssembleMatrix() or AssembleVector() which set mAssembleMatrix, mAssembleVector
     * accordingly. Pure and therefore is implemented in child classes. Will involve looping
     * over elements (which may be volume, surface or cable elements), and computing
     * integrals and adding them to the vector or matrix
     */
    virtual void DoAssemble()=0;

public:

    /**
     * Constructor.
     */
    AbstractFeAssemblerInterface();

    /**
     * Set the matrix that needs to be assembled. Requires CAN_ASSEMBLE_MATRIX==true.
     *
     * @param rMatToAssemble Reference to the matrix
     * @param zeroMatrixBeforeAssembly Whether to zero the matrix before assembling
     *  (otherwise it is just added to)
     */
    void SetMatrixToAssemble(Mat& rMatToAssemble, bool zeroMatrixBeforeAssembly=true);

    /**
     * Set the vector that needs to be assembled. Requires CAN_ASSEMBLE_VECTOR==true.
     *
     * @param rVecToAssemble Reference to the vector
     * @param zeroVectorBeforeAssembly Whether to zero the vector before assembling
     *  (otherwise it is just added to)
     */
    void SetVectorToAssemble(Vec& rVecToAssemble, bool zeroVectorBeforeAssembly);

    /**
     * Assemble everything that the class can assemble.
     */
    void Assemble()
    {
        mAssembleMatrix = CAN_ASSEMBLE_MATRIX;
        mAssembleVector = CAN_ASSEMBLE_VECTOR;
        DoAssemble();
    }

    /**
     * Assemble the matrix. Requires CAN_ASSEMBLE_MATRIX==true
     */
    void AssembleMatrix()
    {
        assert(CAN_ASSEMBLE_MATRIX);
        mAssembleMatrix = true;
        mAssembleVector = false;
        DoAssemble();
    }

    /**
     * Assemble the vector. Requires CAN_ASSEMBLE_VECTOR==true
     */
    void AssembleVector()
    {
        assert(CAN_ASSEMBLE_VECTOR);
        mAssembleMatrix = false;
        mAssembleVector = true;
        DoAssemble();
    }

    /**
     * Destructor.
     */
    virtual ~AbstractFeAssemblerInterface()
    {
    }
};

template <bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX>
AbstractFeAssemblerInterface<CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX>::AbstractFeAssemblerInterface()
    : mVectorToAssemble(nullptr),
      mMatrixToAssemble(nullptr),
      mZeroMatrixBeforeAssembly(true),
      mZeroVectorBeforeAssembly(true)
{
    assert(CAN_ASSEMBLE_VECTOR || CAN_ASSEMBLE_MATRIX);
}

template <bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX>
void AbstractFeAssemblerInterface<CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX>::SetMatrixToAssemble(Mat& rMatToAssemble, bool zeroMatrixBeforeAssembly)
{
    assert(rMatToAssemble);
    MatGetOwnershipRange(rMatToAssemble, &mOwnershipRangeLo, &mOwnershipRangeHi);

    mMatrixToAssemble = rMatToAssemble;
    mZeroMatrixBeforeAssembly = zeroMatrixBeforeAssembly;
}

template <bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX>
void AbstractFeAssemblerInterface<CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX>::SetVectorToAssemble(Vec& rVecToAssemble, bool zeroVectorBeforeAssembly)
{
    assert(rVecToAssemble);
    VecGetOwnershipRange(rVecToAssemble, &mOwnershipRangeLo, &mOwnershipRangeHi);

    mVectorToAssemble = rVecToAssemble;
    mZeroVectorBeforeAssembly = zeroVectorBeforeAssembly;
}

#endif // ABSTRACTFEASSEMBLERINTERFACE_HPP_
