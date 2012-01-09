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
class AbstractFeAssemblerInterface  : boost::noncopyable
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
    : mVectorToAssemble(NULL),
      mMatrixToAssemble(NULL),
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
