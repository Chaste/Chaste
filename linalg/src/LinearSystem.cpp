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

#include "LinearSystem.hpp"

#include <iostream>
#include <cassert>
#include <algorithm>
#include <boost/scoped_array.hpp>

#include "PetscException.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "HeartEventHandler.hpp"
#include "Timer.hpp"
#include "Warnings.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

LinearSystem::LinearSystem(PetscInt lhsVectorSize, unsigned rowPreallocation)
   :mPrecondMatrix(nullptr),
    mSize(lhsVectorSize),
    mMatNullSpace(nullptr),
    mDestroyMatAndVec(true),
    mKspIsSetup(false),
    mNonZerosUsed(0.0),
    mMatrixIsConstant(false),
    mTolerance(1e-6),
    mUseAbsoluteTolerance(false),
    mDirichletBoundaryConditionsVector(nullptr),
    mpBlockDiagonalPC(nullptr),
    mpLDUFactorisationPC(nullptr),
    mpTwoLevelsBlockDiagonalPC(nullptr),
    mpBathNodes( boost::shared_ptr<std::vector<PetscInt> >() ),
    mPrecondMatrixIsNotLhs(false),
    mRowPreallocation(rowPreallocation),
    mUseFixedNumberIterations(false),
    mEvaluateNumItsEveryNSolves(UINT_MAX),
    mpConvergenceTestContext(nullptr),
    mEigMin(DBL_MAX),
    mEigMax(DBL_MIN),
    mForceSpectrumReevaluation(false)
{
    assert(lhsVectorSize > 0);
    if (mRowPreallocation == UINT_MAX)
    {
        // Automatic preallocation if it's a small matrix
        if (lhsVectorSize<15)
        {
            mRowPreallocation=lhsVectorSize;
        }
        else
        {
            EXCEPTION("You must provide a rowPreallocation argument for a large sparse system");
        }
    }

    mRhsVector = PetscTools::CreateVec(mSize);
    PetscTools::SetupMat(mLhsMatrix, mSize, mSize, mRowPreallocation, PETSC_DECIDE, PETSC_DECIDE);

    VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);

    /// \todo: if we create a linear system object outside a cardiac solver, these are gonna
    /// be the default solver and preconditioner. Not consitent with ChasteDefaults.xml though...
    mKspType = "gmres";
    mPcType = "jacobi";

    mNumSolves = 0;
#ifdef TRACE_KSP
    mTotalNumIterations = 0;
    mMaxNumIterations = 0;
#endif
}

LinearSystem::LinearSystem(PetscInt lhsVectorSize, Mat lhsMatrix, Vec rhsVector)
   :mPrecondMatrix(nullptr),
    mSize(lhsVectorSize),
    mMatNullSpace(nullptr),
    mDestroyMatAndVec(true),
    mKspIsSetup(false),
    mNonZerosUsed(0.0),
    mMatrixIsConstant(false),
    mTolerance(1e-6),
    mUseAbsoluteTolerance(false),
    mDirichletBoundaryConditionsVector(nullptr),
    mpBlockDiagonalPC(nullptr),
    mpLDUFactorisationPC(nullptr),
    mpTwoLevelsBlockDiagonalPC(nullptr),
    mpBathNodes( boost::shared_ptr<std::vector<PetscInt> >() ),
    mPrecondMatrixIsNotLhs(false),
    mUseFixedNumberIterations(false),
    mEvaluateNumItsEveryNSolves(UINT_MAX),
    mpConvergenceTestContext(nullptr),
    mEigMin(DBL_MAX),
    mEigMax(DBL_MIN),
    mForceSpectrumReevaluation(false)
{
    assert(lhsVectorSize > 0);
    // Conveniently, PETSc Mats and Vecs are actually pointers
    mLhsMatrix = lhsMatrix;
    mRhsVector = rhsVector;

    VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);

    mNumSolves = 0;
#ifdef TRACE_KSP
    mTotalNumIterations = 0;
    mMaxNumIterations = 0;
#endif
}

LinearSystem::LinearSystem(Vec templateVector, unsigned rowPreallocation, bool newAllocationError)
   :mPrecondMatrix(nullptr),
    mMatNullSpace(nullptr),
    mDestroyMatAndVec(true),
    mKspIsSetup(false),
    mMatrixIsConstant(false),
    mTolerance(1e-6),
    mUseAbsoluteTolerance(false),
    mDirichletBoundaryConditionsVector(nullptr),
    mpBlockDiagonalPC(nullptr),
    mpLDUFactorisationPC(nullptr),
    mpTwoLevelsBlockDiagonalPC(nullptr),
    mpBathNodes( boost::shared_ptr<std::vector<PetscInt> >() ),
    mPrecondMatrixIsNotLhs(false),
    mRowPreallocation(rowPreallocation),
    mUseFixedNumberIterations(false),
    mEvaluateNumItsEveryNSolves(UINT_MAX),
    mpConvergenceTestContext(nullptr),
    mEigMin(DBL_MAX),
    mEigMax(DBL_MIN),
    mForceSpectrumReevaluation(false)
{
    VecDuplicate(templateVector, &mRhsVector);
    VecGetSize(mRhsVector, &mSize);
    VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);
    PetscInt local_size = mOwnershipRangeHi - mOwnershipRangeLo;

    PetscTools::SetupMat(mLhsMatrix, mSize, mSize, mRowPreallocation, local_size, local_size, true, newAllocationError);

    /// \todo: if we create a linear system object outside a cardiac solver, these are gonna
    /// be the default solver and preconditioner. Not consistent with ChasteDefaults.xml though...
    mKspType = "gmres";
    mPcType = "jacobi";

    mNumSolves = 0;
#ifdef TRACE_KSP
    mTotalNumIterations = 0;
    mMaxNumIterations = 0;
#endif
}

LinearSystem::LinearSystem(Vec residualVector, Mat jacobianMatrix)
   :mPrecondMatrix(nullptr),
    mMatNullSpace(nullptr),
    mDestroyMatAndVec(false),
    mKspIsSetup(false),
    mMatrixIsConstant(false),
    mTolerance(1e-6),
    mUseAbsoluteTolerance(false),
    mDirichletBoundaryConditionsVector(nullptr),
    mpBlockDiagonalPC(nullptr),
    mpLDUFactorisationPC(nullptr),
    mpTwoLevelsBlockDiagonalPC(nullptr),
    mpBathNodes( boost::shared_ptr<std::vector<PetscInt> >() ),
    mPrecondMatrixIsNotLhs(false),
    mRowPreallocation(UINT_MAX),
    mUseFixedNumberIterations(false),
    mEvaluateNumItsEveryNSolves(UINT_MAX),
    mpConvergenceTestContext(nullptr),
    mEigMin(DBL_MAX),
    mEigMax(DBL_MIN),
    mForceSpectrumReevaluation(false)
{
    assert(residualVector || jacobianMatrix);
    mRhsVector = residualVector;
    mLhsMatrix = jacobianMatrix;

    PetscInt mat_size=0, vec_size=0;
    if (mRhsVector)
    {
        VecGetSize(mRhsVector, &vec_size);
        mSize = (unsigned)vec_size;
        VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);
    }
    if (mLhsMatrix)
    {
        PetscInt mat_cols;
        MatGetSize(mLhsMatrix, &mat_size, &mat_cols);
        assert(mat_size == mat_cols);
        mSize = (unsigned)mat_size;
        MatGetOwnershipRange(mLhsMatrix, &mOwnershipRangeLo, &mOwnershipRangeHi);

        MatInfo matrix_info;
        MatGetInfo(mLhsMatrix, MAT_GLOBAL_MAX, &matrix_info);

        /*
         * Assuming that mLhsMatrix was created with PetscTools::SetupMat, the value
         * below should be equivalent to what was used as preallocation in that call.
         */
        mRowPreallocation = (unsigned) matrix_info.nz_allocated / mSize;
    }
    assert(!mRhsVector || !mLhsMatrix || vec_size == mat_size);

    /// \todo: if we create a linear system object outside a cardiac solver, these are gonna
    /// be the default solver and preconditioner. Not consitent with ChasteDefaults.xml though...
    mKspType = "gmres";
    mPcType = "jacobi";

    mNumSolves = 0;
#ifdef TRACE_KSP
    mTotalNumIterations = 0;
    mMaxNumIterations = 0;
#endif
}

LinearSystem::~LinearSystem()
{
    delete mpBlockDiagonalPC;
    delete mpLDUFactorisationPC;
    delete mpTwoLevelsBlockDiagonalPC;

    if (mDestroyMatAndVec)
    {
        PetscTools::Destroy(mRhsVector);
        PetscTools::Destroy(mLhsMatrix);
    }

    if (mPrecondMatrixIsNotLhs)
    {
        PetscTools::Destroy(mPrecondMatrix);
    }

    if (mMatNullSpace)
    {
        MatNullSpaceDestroy(PETSC_DESTROY_PARAM(mMatNullSpace));
    }

    if (mKspIsSetup)
    {
        KSPDestroy(PETSC_DESTROY_PARAM(mKspSolver));
    }

    if (mDirichletBoundaryConditionsVector)
    {
        ///\todo Never tested in linalg component
        PetscTools::Destroy(mDirichletBoundaryConditionsVector);
    }

#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
    if (mpConvergenceTestContext)
    {
#if (PETSC_VERSION_MINOR >= 5) //PETSc 3.5 or later
        KSPConvergedDefaultDestroy(mpConvergenceTestContext);
#else
        KSPDefaultConvergedDestroy(mpConvergenceTestContext);
#endif
    }
#endif

#ifdef TRACE_KSP
    if (PetscTools::AmMaster())
    {
        if (mNumSolves > 0)
        {
            double ave_num_iterations = mTotalNumIterations/(double)mNumSolves;

            std::cout << std::endl << "KSP iterations report:" << std::endl;
            std::cout << "mNumSolves" << "\t" << "mTotalNumIterations" << "\t" << "mMaxNumIterations" << "\t" << "mAveNumIterations" << std::endl;
            std::cout << mNumSolves << "\t" << mTotalNumIterations << "\t" << mMaxNumIterations << "\t" << ave_num_iterations << std::endl;
        }
    }
#endif
}

void LinearSystem::SetMatrixElement(PetscInt row, PetscInt col, double value)
{
    PetscMatTools::SetElement(mLhsMatrix, row, col, value);
}

void LinearSystem::AddToMatrixElement(PetscInt row, PetscInt col, double value)
{
    PetscMatTools::AddToElement(mLhsMatrix, row, col, value);
}

void LinearSystem::AssembleFinalLinearSystem()
{
    FinaliseLhsMatrix();
    FinaliseRhsVector();
}

void LinearSystem::AssembleIntermediateLinearSystem()
{
    SwitchWriteModeLhsMatrix();
    FinaliseRhsVector();
}

void LinearSystem::FinaliseLhsMatrix()
{
    PetscMatTools::Finalise(mLhsMatrix);
}

void LinearSystem::SwitchWriteModeLhsMatrix()
{
    PetscMatTools::SwitchWriteMode(mLhsMatrix);
}

void LinearSystem::FinalisePrecondMatrix()
{
    PetscMatTools::Finalise(mPrecondMatrix);
}

void LinearSystem::FinaliseRhsVector()
{
    PetscVecTools::Finalise(mRhsVector);
}

void LinearSystem::SetRhsVectorElement(PetscInt row, double value)
{
    PetscVecTools::SetElement(mRhsVector, row, value);
}

void LinearSystem::AddToRhsVectorElement(PetscInt row, double value)
{
    PetscVecTools::AddToElement(mRhsVector, row, value);
}

void LinearSystem::DisplayMatrix()
{
    PetscMatTools::Display(mLhsMatrix);
}

void LinearSystem::DisplayRhs()
{
    PetscVecTools::Display(mRhsVector);
}

void LinearSystem::SetMatrixRow(PetscInt row, double value)
{
    PetscMatTools::SetRow(mLhsMatrix, row, value);
}

Vec LinearSystem::GetMatrixRowDistributed(unsigned rowIndex)
{
    return PetscMatTools::GetMatrixRowDistributed(mLhsMatrix, rowIndex);
}

void LinearSystem::ZeroMatrixRowsWithValueOnDiagonal(std::vector<unsigned>& rRows, double diagonalValue)
{
    PetscMatTools::ZeroRowsWithValueOnDiagonal(mLhsMatrix, rRows, diagonalValue);
}

void LinearSystem::ZeroMatrixRowsAndColumnsWithValueOnDiagonal(std::vector<unsigned>& rRowColIndices, double diagonalValue)
{
    PetscMatTools::ZeroRowsAndColumnsWithValueOnDiagonal(mLhsMatrix, rRowColIndices, diagonalValue);
}

void LinearSystem::ZeroMatrixColumn(PetscInt col)
{
    PetscMatTools::ZeroColumn(mLhsMatrix, col);
}

void LinearSystem::ZeroRhsVector()
{
    PetscVecTools::Zero(mRhsVector);
}

void LinearSystem::ZeroLhsMatrix()
{
    PetscMatTools::Zero(mLhsMatrix);
}

void LinearSystem::ZeroLinearSystem()
{
    ZeroRhsVector();
    ZeroLhsMatrix();
}

unsigned LinearSystem::GetSize() const
{
    return (unsigned) mSize;
}

void LinearSystem::SetNullBasis(Vec nullBasis[], unsigned numberOfBases)
{
#ifndef NDEBUG
    // Check all the vectors of the base are normal
    for (unsigned vec_index=0; vec_index<numberOfBases; vec_index++)
    {
        PetscReal l2_norm;
        VecNorm(nullBasis[vec_index], NORM_2, &l2_norm);
        if (fabs(l2_norm-1.0) > 1e-08)
        {
            EXCEPTION("One of the vectors in the null space is not normalised");
        }
    }

    // Check all the vectors of the base are orthogonal
    for (unsigned vec_index=1; vec_index<numberOfBases; vec_index++)
    {
        // The strategy is to check the (i-1)-th vector against vectors from i to n with VecMDot. This should be the most efficient way of doing it.
        unsigned num_vectors_ahead = numberOfBases-vec_index;
        boost::scoped_array<PetscScalar> dot_products(new PetscScalar[num_vectors_ahead]);
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
        VecMDot(num_vectors_ahead, nullBasis[vec_index-1], &nullBasis[vec_index], dot_products.get());
#else
        VecMDot(nullBasis[vec_index-1], num_vectors_ahead, &nullBasis[vec_index], dot_products.get());
#endif
        for (unsigned index=0; index<num_vectors_ahead; index++)
        {
            if (fabs(dot_products[index]) > 1e-08 )
            {
                EXCEPTION("The null space is not orthogonal.");
            }
        }
    }

#endif

    PETSCEXCEPT( MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, numberOfBases, nullBasis, &mMatNullSpace) );
}

void LinearSystem::RemoveNullSpace()
{
    // Only remove if previously set
    if (mMatNullSpace)
    {
        PETSCEXCEPT( MatNullSpaceDestroy(PETSC_DESTROY_PARAM(mMatNullSpace)) );
        PETSCEXCEPT( MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 0, nullptr, &mMatNullSpace) );
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3) //PETSc 3.3 or later
        // Setting null space in the KSP was deprecated in PETSc 3.6, but setting the null space
        // for the matrix appeared in PETSc 3.3 so 3.3, 3.4, 3.5 can do either
        PETSCEXCEPT( MatSetNullSpace(mLhsMatrix, mMatNullSpace) );
        /* This is heavy-handed:
         * Changing the Null-space appears to only work if we reset the KSP
         * Adding null-space to the matrix has to happen *before* KSPSetOperators
         */
        ResetKspSolver();
#else
        if (mKspIsSetup)
        {
            PETSCEXCEPT( KSPSetNullSpace(mKspSolver, mMatNullSpace) );
        }
        //else: it will be set next time Solve() is called
#endif
    }
}

void LinearSystem::GetOwnershipRange(PetscInt& lo, PetscInt& hi)
{
    lo = mOwnershipRangeLo;
    hi = mOwnershipRangeHi;
}

double LinearSystem::GetMatrixElement(PetscInt row, PetscInt col)
{
    return PetscMatTools::GetElement(mLhsMatrix, row, col);
}

double LinearSystem::GetRhsVectorElement(PetscInt row)
{
    return PetscVecTools::GetElement(mRhsVector, row);
}

unsigned LinearSystem::GetNumIterations() const
{
    assert(this->mKspIsSetup);

    PetscInt num_its;
    KSPGetIterationNumber(this->mKspSolver, &num_its);

    return (unsigned) num_its;
}

Vec& LinearSystem::rGetRhsVector()
{
    return mRhsVector;
}

Vec LinearSystem::GetRhsVector() const
{
    return mRhsVector;
}

Mat& LinearSystem::rGetLhsMatrix()
{
    return mLhsMatrix;
}

Mat LinearSystem::GetLhsMatrix() const
{
    return mLhsMatrix;
}

Mat& LinearSystem::rGetPrecondMatrix()
{
    if (!mPrecondMatrixIsNotLhs)
    {
        EXCEPTION("LHS matrix used for preconditioner construction");
    }

    return mPrecondMatrix;
}

Vec& LinearSystem::rGetDirichletBoundaryConditionsVector()
{
    return mDirichletBoundaryConditionsVector;
}

void LinearSystem::SetMatrixIsSymmetric(bool isSymmetric)
{
    /// \todo: shall we allow modifying the symmetry flag anytime?

    if (isSymmetric)
    {
        PetscMatTools::SetOption(mLhsMatrix, MAT_SYMMETRIC);
        PetscMatTools::SetOption(mLhsMatrix, MAT_SYMMETRY_ETERNAL);
    }
    else
    {
// don't have a PetscMatTools method for setting options to false
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        MatSetOption(mLhsMatrix, MAT_SYMMETRIC, PETSC_FALSE);
        MatSetOption(mLhsMatrix, MAT_STRUCTURALLY_SYMMETRIC, PETSC_FALSE);
        MatSetOption(mLhsMatrix, MAT_SYMMETRY_ETERNAL, PETSC_FALSE);
#else
        MatSetOption(mLhsMatrix, MAT_NOT_SYMMETRIC);
        MatSetOption(mLhsMatrix, MAT_NOT_STRUCTURALLY_SYMMETRIC);
        MatSetOption(mLhsMatrix, MAT_NOT_SYMMETRY_ETERNAL);
#endif
    }
}

bool LinearSystem::IsMatrixSymmetric()
{
    PetscBool symmetry_flag_is_set;
    PetscBool symmetry_flag;

    MatIsSymmetricKnown(mLhsMatrix, &symmetry_flag_is_set, &symmetry_flag);

    // If the flag is not set we assume is a non-symmetric matrix
    return symmetry_flag_is_set && symmetry_flag;
}

void LinearSystem::SetMatrixIsConstant(bool matrixIsConstant)
{
    mMatrixIsConstant = matrixIsConstant;
}

void LinearSystem::SetRelativeTolerance(double relativeTolerance)
{
    mTolerance = relativeTolerance;
    mUseAbsoluteTolerance = false;
    if (mKspIsSetup)
    {
        KSPSetTolerances(mKspSolver, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    }
}

void LinearSystem::SetAbsoluteTolerance(double absoluteTolerance)
{
    mTolerance = absoluteTolerance;
    mUseAbsoluteTolerance = true;
    if (mKspIsSetup)
    {
        KSPSetTolerances(mKspSolver, DBL_EPSILON, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT);
    }
}

void LinearSystem::SetKspType(const char *kspType)
{
    mKspType = kspType;
    if (mKspIsSetup)
    {
        KSPSetType(mKspSolver, kspType);
        KSPSetFromOptions(mKspSolver);
    }
}

void LinearSystem::SetPcType(const char* pcType, boost::shared_ptr<std::vector<PetscInt> > pBathNodes)
{
    mPcType = pcType;
    mpBathNodes = pBathNodes;

    if (mKspIsSetup)
    {
        if (mPcType == "blockdiagonal")
        {
            // If the previous preconditioner was purpose-built we need to free the appropriate pointer.
            /// \todo: #1082 use a single pointer to abstract class
            delete mpBlockDiagonalPC;
            mpBlockDiagonalPC = nullptr;
            delete mpLDUFactorisationPC;
            mpLDUFactorisationPC = nullptr;
            delete mpTwoLevelsBlockDiagonalPC;
            mpTwoLevelsBlockDiagonalPC = nullptr;

            mpBlockDiagonalPC = new PCBlockDiagonal(mKspSolver);
        }
        else if (mPcType == "ldufactorisation")
        {
            // If the previous preconditioner was purpose-built we need to free the appropriate pointer.
            /// \todo: #1082 use a single pointer to abstract class
            delete mpBlockDiagonalPC;
            mpBlockDiagonalPC = nullptr;
            delete mpLDUFactorisationPC;
            mpLDUFactorisationPC = nullptr;
            delete mpTwoLevelsBlockDiagonalPC;
            mpTwoLevelsBlockDiagonalPC = nullptr;

            mpLDUFactorisationPC = new PCLDUFactorisation(mKspSolver);
        }
        else if (mPcType == "twolevelsblockdiagonal")
        {
            // If the previous preconditioner was purpose-built we need to free the appropriate pointer.
            /// \todo: #1082 use a single pointer to abstract class
            delete mpBlockDiagonalPC;
            mpBlockDiagonalPC = nullptr;
            delete mpLDUFactorisationPC;
            mpLDUFactorisationPC = nullptr;
            delete mpTwoLevelsBlockDiagonalPC;
            mpTwoLevelsBlockDiagonalPC = nullptr;

            if (!mpBathNodes)
            {
                TERMINATE("You must provide a list of bath nodes when using TwoLevelsBlockDiagonalPC"); // LCOV_EXCL_LINE
            }
            mpTwoLevelsBlockDiagonalPC = new PCTwoLevelsBlockDiagonal(mKspSolver, *mpBathNodes);
        }
        else
        {
            PC prec;
            KSPGetPC(mKspSolver, &prec);
            PCSetType(prec, pcType);
        }
        KSPSetFromOptions(mKspSolver);
    }
}

Vec LinearSystem::Solve(Vec lhsGuess)
{
    /*
     * The following lines are very useful for debugging:
     *    MatView(mLhsMatrix,    PETSC_VIEWER_STDOUT_WORLD);
     *    VecView(mRhsVector,    PETSC_VIEWER_STDOUT_WORLD);
     */

    // Double check that the non-zero pattern hasn't changed
    MatInfo mat_info;
    MatGetInfo(mLhsMatrix, MAT_GLOBAL_SUM, &mat_info);

    if (!mKspIsSetup)
    {
        // Create PETSc Vec that may be required if we use a Chebyshev solver
        Vec chebyshev_lhs_vector = nullptr;

        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        mNonZerosUsed=mat_info.nz_used;
        //MatNorm(mLhsMatrix, NORM_FROBENIUS, &mMatrixNorm);
        PC prec; //Type of pre-conditioner

        KSPCreate(PETSC_COMM_WORLD, &mKspSolver);

        if (mMatNullSpace) // Adding null-space to the matrix (new style) has to happen *before* KSPSetOperators
        {
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3) //PETSc 3.3 or later
            // Setting null space in the KSP was deprecated in PETSc 3.6, but setting the null space
            // for the matrix appeared in PETSc 3.3 so 3.3, 3.4, 3.5 can do either

            PETSCEXCEPT(MatSetNullSpace(mLhsMatrix, mMatNullSpace));
#else
            PETSCEXCEPT(KSPSetNullSpace(mKspSolver, mMatNullSpace));
#endif
        }
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5) //PETSc 3.5 or later
        // Do nothing.  Note that reusing the pre-conditioner in later PETSc versions is done after the pre-conditioner is formed
        // (This comment is retained here so that the #if logic is consistent.)
#else
        /*
         * The preconditioner flag (last argument) in the following calls says
         * how to reuse the preconditioner on subsequent iterations.
         */
        MatStructure preconditioner_over_successive_calls;

        if (mMatrixIsConstant)
        {
            preconditioner_over_successive_calls = SAME_PRECONDITIONER;
        }
        else
        {
            preconditioner_over_successive_calls = SAME_NONZERO_PATTERN;
        }
#endif

        if (mPrecondMatrixIsNotLhs)
        {
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5) //PETSc 3.5 or later
            KSPSetOperators(mKspSolver, mLhsMatrix, mPrecondMatrix);
#else
            KSPSetOperators(mKspSolver, mLhsMatrix, mPrecondMatrix, preconditioner_over_successive_calls);
#endif
        }
        else
        {
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5) //PETSc 3.5 or later
            KSPSetOperators(mKspSolver, mLhsMatrix, mLhsMatrix);
#else
            KSPSetOperators(mKspSolver, mLhsMatrix, mLhsMatrix, preconditioner_over_successive_calls);
#endif
        }
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 5) //PETSc 3.5 or later
        if (mMatrixIsConstant)
        {
            // Emulate SAME_PRECONDITIONER as above
            KSPSetReusePreconditioner(mKspSolver, PETSC_TRUE);
        }
#endif

        // Set either absolute or relative tolerance of the KSP solver.
        // The default is to use relative tolerance (1e-6)
        if (mUseAbsoluteTolerance)
        {
            KSPSetTolerances(mKspSolver, DBL_EPSILON, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT);
        }
        else
        {
            KSPSetTolerances(mKspSolver, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
        }

        // Set ksp and pc types
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3) //PETSc 3.3 or later
        if (mKspType == "chebychev")
        {
            KSPSetType(mKspSolver, "chebyshev");
        }
        else
        {
            KSPSetType(mKspSolver, mKspType.c_str());
        }
#else
        KSPSetType(mKspSolver, mKspType.c_str());
#endif
        KSPGetPC(mKspSolver, &prec);

        // Turn off pre-conditioning if the system size is very small
        const bool is_small = (mSize <= 6); ///\todo This is a magic number.  Do we want a warning here?
        if (is_small)
        {
            PCSetType(prec, PCNONE);
        }
        else
        {
#ifdef TRACE_KSP
            Timer::Reset();
#endif
            if (mPcType == "blockdiagonal")
            {
                mpBlockDiagonalPC = new PCBlockDiagonal(mKspSolver);
#ifdef TRACE_KSP
                if (PetscTools::AmMaster())
                {
                    Timer::Print("Purpose-build preconditioner creation");
                }
#endif
            }
            else if (mPcType == "ldufactorisation")
            {
                mpLDUFactorisationPC = new PCLDUFactorisation(mKspSolver);
#ifdef TRACE_KSP
                if (PetscTools::AmMaster())
                {
                    Timer::Print("Purpose-build preconditioner creation");
                }
#endif
            }
            else if (mPcType == "twolevelsblockdiagonal")
            {
                if (!mpBathNodes)
                {
                    TERMINATE("You must provide a list of bath nodes when using TwoLevelsBlockDiagonalPC"); // LCOV_EXCL_LINE
                }
                mpTwoLevelsBlockDiagonalPC = new PCTwoLevelsBlockDiagonal(mKspSolver, *mpBathNodes);
#ifdef TRACE_KSP
                if (PetscTools::AmMaster())
                {
                    Timer::Print("Purpose-build preconditioner creation");
                }
#endif

            }
            else
            {
                PCSetType(prec, mPcType.c_str());
            }
        }

        KSPSetFromOptions(mKspSolver);

        if (lhsGuess)
        {
            // Assume that the user of this method will always be kind enough to give us a reasonable guess.
            KSPSetInitialGuessNonzero(mKspSolver,PETSC_TRUE);
        }
        /*
         * Non-adaptive Chebyshev: the required spectrum approximation is computed just once
         * at the beginning of the simulation. This is done with two extra CG solves.
         */
        if (mKspType == "chebychev" && !mUseFixedNumberIterations)
        {
#ifdef TRACE_KSP
            Timer::Reset();
#endif
            // Preconditioned matrix spectrum is approximated with CG
            KSPSetType(mKspSolver,"cg");
            KSPSetComputeEigenvalues(mKspSolver, PETSC_TRUE);
            KSPSetUp(mKspSolver);

            VecDuplicate(mRhsVector, &chebyshev_lhs_vector);
            if (lhsGuess)
            {
                VecCopy(lhsGuess, chebyshev_lhs_vector);
            }

            // Smallest eigenvalue is approximated to default tolerance
            KSPSolve(mKspSolver, mRhsVector, chebyshev_lhs_vector);

            PetscReal *r_eig = new PetscReal[mSize];
            PetscReal *c_eig = new PetscReal[mSize];
            PetscInt eigs_computed;
            KSPComputeEigenvalues(mKspSolver, mSize, r_eig, c_eig, &eigs_computed);

            mEigMin = r_eig[0];

            // Largest eigenvalue is approximated to machine precision
            KSPSetTolerances(mKspSolver, DBL_EPSILON, DBL_EPSILON, PETSC_DEFAULT, PETSC_DEFAULT);
            KSPSetUp(mKspSolver);
            if (lhsGuess)
            {
                VecCopy(lhsGuess, chebyshev_lhs_vector);
            }

            KSPSolve(mKspSolver, mRhsVector, chebyshev_lhs_vector);
            KSPComputeEigenvalues(mKspSolver, mSize, r_eig, c_eig, &eigs_computed);

            mEigMax = r_eig[eigs_computed-1];

#ifdef TRACE_KSP
            /*
             * Under certain circumstances (see Golub&Overton 1988), underestimating
             * the spectrum of the preconditioned operator improves convergence rate.
             * See publication for a discussion and for definition of alpha and sigma_one.
             */
            if (PetscTools::AmMaster())
            {
                std::cout << "EIGS: ";
                for (int index=0; index<eigs_computed; index++)
                {
                    std::cout << r_eig[index] << ", ";
                }
                std::cout << std::endl;
            }

            if (PetscTools::AmMaster()) std::cout << "EIGS "<< mEigMax << " " << mEigMin <<std::endl;
            double alpha = 2/(mEigMax+mEigMin);
            double sigma_one = 1 - alpha*mEigMin;
            if (PetscTools::AmMaster()) std::cout << "sigma_1 = 1 - alpha*mEigMin = "<< sigma_one <<std::endl;
#endif

            // Set Chebyshev solver and max/min eigenvalues
            assert(mKspType == "chebychev");
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3) //PETSc 3.3 or later
            KSPSetType(mKspSolver, "chebyshev");
            KSPChebyshevSetEigenvalues(mKspSolver, mEigMax, mEigMin);
#else
            KSPSetType(mKspSolver, mKspType.c_str());
            KSPChebychevSetEigenvalues(mKspSolver, mEigMax, mEigMin);
#endif
            KSPSetComputeEigenvalues(mKspSolver, PETSC_FALSE);
            if (mUseAbsoluteTolerance)
            {
                KSPSetTolerances(mKspSolver, DBL_EPSILON, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT);
            }
            else
            {
                KSPSetTolerances(mKspSolver, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
            }

            delete[] r_eig;
            delete[] c_eig;

#ifdef TRACE_KSP
            if (PetscTools::AmMaster())
            {
                Timer::Print("Computing extremal eigenvalues");
            }
#endif
        }

#ifdef TRACE_KSP
        Timer::Reset();
#endif

        KSPSetUp(mKspSolver);

        if (chebyshev_lhs_vector)
        {
            PetscTools::Destroy(chebyshev_lhs_vector);
        }

#ifdef TRACE_KSP
        if (PetscTools::AmMaster())
        {
            Timer::Print("KSPSetUP (contains preconditioner creation for PETSc preconditioners)");
        }
#endif

        mKspIsSetup = true;

        HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
    }
    else
    {
        if (mNonZerosUsed!=mat_info.nz_used)
        {
            WARNING("LinearSystem doesn't like the non-zero pattern of a matrix to change. (I think you changed it).");
            mNonZerosUsed = mat_info.nz_used;
        }
//        PetscScalar norm;
//        MatNorm(mLhsMatrix, NORM_FROBENIUS, &norm);
//        if (fabs(norm - mMatrixNorm) > 0)
//        {
//            EXCEPTION("LinearSystem doesn't allow the matrix norm to change");
//        }
    }

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    // Create solution vector
    ///\todo Should it be compulsory for the caller to supply this and manage the memory?
    Vec lhs_vector;
    VecDuplicate(mRhsVector, &lhs_vector); // Sets the same size (doesn't copy)
    if (lhsGuess)
    {
        VecCopy(lhsGuess, lhs_vector);
        // If this wasn't done at construction time then it may be too late for this:
        // KSPSetInitialGuessNonzero(mKspSolver, PETSC_TRUE);
        // Is it possible to warn the user?
    }

    // Check if the right hand side is small (but non-zero), PETSc can diverge immediately
    // with a non-zero initial guess. Here we check for this and alter the initial guess to zero.
    PetscReal rhs_norm;
    VecNorm(mRhsVector, NORM_1, &rhs_norm);
    double rhs_zero_tol = 1e-15;

    if (rhs_norm < rhs_zero_tol && rhs_norm > 0.0)
    {
        WARNING("Using zero initial guess due to small right hand side vector");
        PetscVecTools::Zero(lhs_vector);
    }

    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
//    // Double check that the mRhsVector contains sensible values
//    double* p_rhs,* p_guess;
//    VecGetArray(mRhsVector, &p_rhs);
//    VecGetArray(lhsGuess, &p_guess);
//    for (int global_index=mOwnershipRangeLo; global_index<mOwnershipRangeHi; global_index++)
//    {
//        int local_index = global_index - mOwnershipRangeLo;
//        assert(!std::isnan(p_rhs[local_index]));
//        assert(!std::isnan(p_guess[local_index]));
//        if (p_rhs[local_index] != p_rhs[local_index])
//        {
//            std::cout << "********* PETSc nan\n";
//            assert(0);
//        }
//    }
//    std::cout << "b[0] = " << p_rhs[0] << ", guess[0] = " << p_guess[0] << "\n";
//    VecRestoreArray(mRhsVector, &p_rhs);
//    VecRestoreArray(lhsGuess, &p_guess);
//    // Test A*guess
//    Vec temp;
//    VecDuplicate(mRhsVector, &temp);
//    MatMult(mLhsMatrix, lhs_vector, temp);
//    double* p_temp;
//    VecGetArray(temp, &p_temp);
//    std::cout << "temp[0] = " << p_temp[0] << "\n";
//    VecRestoreArray(temp, &p_temp);
//    PetscTools::Destroy(temp);
//    // Dump the matrix to file
//    PetscViewer viewer;
//    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.output",&viewer);
//    MatView(mLhsMatrix, viewer);
//    PetscViewerFlush(viewer);
//    PetscViewerDestroy(PETSC_DESTROY_PARAM(viewer);
//    // Dump the rhs vector to file
//    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vec.output",&viewer);
//    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
//    VecView(mRhsVector, viewer);
//    PetscViewerFlush(viewer);
//    PetscViewerDestroy(PETSC_DESTROY_PARAM(viewer);

    try
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::SOLVE_LINEAR_SYSTEM);

#ifdef TRACE_KSP
        Timer::Reset();
#endif

        // Current solve has to be done with tolerance-based stop criteria in order to record iterations taken
        if (mUseFixedNumberIterations && (mNumSolves%mEvaluateNumItsEveryNSolves==0 || mForceSpectrumReevaluation))
        {
#if ((PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR <= 2))
            KSPSetNormType(mKspSolver, KSP_PRECONDITIONED_NORM);
#else
            KSPSetNormType(mKspSolver, KSP_NORM_PRECONDITIONED);
#endif

#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
            if (!mpConvergenceTestContext)
            {
    #if (PETSC_VERSION_MINOR >= 5) //PETSc 3.5 or later
                KSPConvergedDefaultCreate(&mpConvergenceTestContext);
    #else
                KSPDefaultConvergedCreate(&mpConvergenceTestContext);
    #endif
            }
    #if (PETSC_VERSION_MINOR >= 5) //PETSc 3.5 or later
            // From PETSc 3.5, KSPDefaultConverged became KSPConvergedDefault.
            KSPSetConvergenceTest(mKspSolver, KSPConvergedDefault, &mpConvergenceTestContext, PETSC_NULL);
    #else
            KSPSetConvergenceTest(mKspSolver, KSPDefaultConverged, &mpConvergenceTestContext, PETSC_NULL);
    #endif

#else
            KSPSetConvergenceTest(mKspSolver, KSPDefaultConverged, PETSC_NULL);
#endif

            if (mUseAbsoluteTolerance)
            {
                KSPSetTolerances(mKspSolver, DBL_EPSILON, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT);
            }
            else
            {
                KSPSetTolerances(mKspSolver, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
            }

            /// \todo #1695 Store this number in a member variable.
            std::stringstream num_it_str;
            num_it_str << 1000;
            PetscTools::SetOption("-ksp_max_it", num_it_str.str().c_str());

            // Adaptive Chebyshev: reevaluate spectrum with cg
            if (mKspType == "chebychev")
            {
                // You can estimate preconditioned matrix spectrum with CG
                KSPSetType(mKspSolver,"cg");
                KSPSetComputeEigenvalues(mKspSolver, PETSC_TRUE);
            }

            KSPSetFromOptions(mKspSolver);
            KSPSetUp(mKspSolver);
        }

        PETSCEXCEPT(KSPSolve(mKspSolver, mRhsVector, lhs_vector));
        HeartEventHandler::EndEvent(HeartEventHandler::SOLVE_LINEAR_SYSTEM);

#ifdef TRACE_KSP
        PetscInt num_it;
        KSPGetIterationNumber(mKspSolver, &num_it);
        if (PetscTools::AmMaster())
        {
            std::cout << "++ Solve: " << mNumSolves << " NumIterations: " << num_it << " "; // don't add std::endl so we get Timer::Print output in the same line (better for grep-ing)
            Timer::Print("Solve");
        }

        mTotalNumIterations += num_it;
        if ((unsigned) num_it > mMaxNumIterations)
        {
            mMaxNumIterations = num_it;
        }
#endif

        // Check that solver converged and throw if not
        KSPConvergedReason reason;
        KSPGetConvergedReason(mKspSolver, &reason);

        if (mUseFixedNumberIterations && PETSC_VERSION_MAJOR < 3)
        {
            WARNING("Not explicitly checking convergence reason when using fixed number of iterations and PETSc 2");
        }
        else
        {
            KSPEXCEPT(reason);
        }

        if (mUseFixedNumberIterations && (mNumSolves%mEvaluateNumItsEveryNSolves==0 || mForceSpectrumReevaluation))
        {
            // Adaptive Chebyshev: reevaluate spectrum with cg
            if (mKspType == "chebychev")
            {
                PetscReal *r_eig = new PetscReal[mSize];
                PetscReal *c_eig = new PetscReal[mSize];
                PetscInt eigs_computed;
                KSPComputeEigenvalues(mKspSolver, mSize, r_eig, c_eig, &eigs_computed);

                mEigMin = r_eig[0];
                /*
                 * Using max() covers a borderline case found in TestChasteBenchmarksForPreDiCT where there's a big
                 * gap in the spectrum between ~1.2 and ~2.5. Some reevaluations pick up 2.5 and others don't. If it
                 * is not picked up, Chebyshev will diverge after 10 solves or so.
                 */
                mEigMax = std::max(mEigMax,r_eig[eigs_computed-1]);

                delete[] r_eig;
                delete[] c_eig;

                assert(mKspType == "chebychev");
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3) //PETSc 3.3 or later
                KSPSetType(mKspSolver, "chebyshev");
                KSPChebyshevSetEigenvalues(mKspSolver, mEigMax, mEigMin);
#else
                KSPSetType(mKspSolver, mKspType.c_str());
                KSPChebychevSetEigenvalues(mKspSolver, mEigMax, mEigMin);
#endif
                KSPSetComputeEigenvalues(mKspSolver, PETSC_FALSE);
            }

#if ((PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR <= 2))
            if (mKspType == "chebychev")
            {
                // See #1695 for more details.
                EXCEPTION("Chebyshev with fixed number of iterations is known to be broken in PETSc <= 2.3.2");
            }

            KSPSetNormType(mKspSolver, KSP_NO_NORM);
#elif (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
            KSPSetNormType(mKspSolver, KSP_NORM_NONE);
    #if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 7) //PETSc 3.7 or later
              /*
             * Up to PETSc 3.7.2 the above call also turned off the default convergence test.
             * However, in PETSc 3.7.3 (subminor release) this behaviour was removed and so, here,
             * we explicitly add it back again.
             * See
             * https://bitbucket.org/petsc/petsc/commits/eb70c44be3430b039effa3de7e1ca2fab9f75a57
             * This following line of code is actually valid from PETSc 3.5.
             */
            KSPSetConvergenceTest(mKspSolver, KSPConvergedSkip, PETSC_NULL, PETSC_NULL);
    #endif
#else
            KSPSetNormType(mKspSolver, KSP_NORM_NO);
#endif

#if (PETSC_VERSION_MAJOR == 2)
            KSPSetConvergenceTest(mKspSolver, KSPSkipConverged, PETSC_NULL);
#endif

            PetscInt num_it;
            KSPGetIterationNumber(mKspSolver, &num_it);
            std::stringstream num_it_str;
            num_it_str << num_it;
            PetscTools::SetOption("-ksp_max_it", num_it_str.str().c_str());

            KSPSetFromOptions(mKspSolver);
            KSPSetUp(mKspSolver);

            mForceSpectrumReevaluation=false;
        }

        mNumSolves++;

    }
    catch (const Exception& e)
    {
        // Destroy solution vector on error to avoid memory leaks
        PetscTools::Destroy(lhs_vector);
        throw e;
    }

    return lhs_vector;
}

void LinearSystem::SetPrecondMatrixIsDifferentFromLhs(bool precondIsDifferent)
{
    mPrecondMatrixIsNotLhs = precondIsDifferent;

    if (mPrecondMatrixIsNotLhs)
    {
        if (mRowPreallocation == UINT_MAX)
        {
            /*
             * At the time of writing, this line will be reached if the constructor
             * with signature LinearSystem(Vec residualVector, Mat jacobianMatrix) is
             * called with jacobianMatrix=NULL and preconditioning matrix different
             * from lhs is used.
             *
             * If this combination is ever required you will need to work out
             * matrix allocation (mRowPreallocation) here.
             */
            NEVER_REACHED;
        }

        PetscInt local_size = mOwnershipRangeHi - mOwnershipRangeLo;
        PetscTools::SetupMat(mPrecondMatrix, mSize, mSize, mRowPreallocation, local_size, local_size);
    }
}

void LinearSystem::SetUseFixedNumberIterations(bool useFixedNumberIterations, unsigned evaluateNumItsEveryNSolves)
{

    mUseFixedNumberIterations = useFixedNumberIterations;
    mEvaluateNumItsEveryNSolves = evaluateNumItsEveryNSolves;
}

void LinearSystem::ResetKspSolver()
{
    if (mKspIsSetup)
    {
        KSPDestroy(PETSC_DESTROY_PARAM(mKspSolver));
    }

    mKspIsSetup = false;
    mForceSpectrumReevaluation = true;

    /*
     * Reset max number of iterations. This option is stored in the configuration database and
     * explicitely read in with KSPSetFromOptions() everytime a KSP object is created. Therefore,
     * destroying the KSP object will not ensure that it is set back to default.
     */
    /// \todo #1695 Store this number in a member variable.
    std::stringstream num_it_str;
    num_it_str << 1000;
    PetscTools::SetOption("-ksp_max_it", num_it_str.str().c_str());
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(LinearSystem)
