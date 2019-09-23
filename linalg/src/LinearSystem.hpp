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

#ifndef _LINEARSYSTEM_HPP_
#define _LINEARSYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include "UblasCustomFunctions.hpp" // needs to be 'first'

#include "PetscTools.hpp"
#include "PetscVecTools.hpp"
#include "PetscMatTools.hpp"
#include "OutputFileHandler.hpp"
#include "PCBlockDiagonal.hpp"
#include "PCLDUFactorisation.hpp"
#include "PCTwoLevelsBlockDiagonal.hpp"
#include "ArchiveLocationInfo.hpp"
#include <boost/serialization/shared_ptr.hpp>

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscviewer.h>

#include <string>
#include <cassert>

/**
 * Linear System class. Stores and solves a linear equation of the form Ax=b,
 * where A is a square matrix and x and b are column vectors.
 * The class uses PETSc.
 */
class LinearSystem
{
    friend class TestLinearSystem;
    friend class TestPCBlockDiagonal;
    friend class TestPCTwoLevelsBlockDiagonal;
    friend class TestPCLDUFactorisation;
    friend class TestChebyshevIteration;

private:

    Mat mLhsMatrix;  /**< The left-hand side matrix. */
    Mat mPrecondMatrix; /**< The matrix used for preconditioning. */
    Vec mRhsVector;  /**< The right-hand side vector. */
    PetscInt mSize;  /**< The size of the linear system. */

    /**
     * \todo Verify claim that ownership range for Vec and Mat is same.
     * This should only matter for efficiency if the claim is false.
     */
    PetscInt mOwnershipRangeLo; /**< For parallel code.  Stores lowest index of vectors and lowest row of matrix stored locally. */
    PetscInt mOwnershipRangeHi; /**< Stores <b>one more than</b> the highest index stored locally. */

    MatNullSpace mMatNullSpace; /**< PETSc null matrix. */

    /** Whether we need to destroy the PETSc matrix and vector in our destructor */
    bool mDestroyMatAndVec;

    KSP mKspSolver;   /**< The PETSc linear solver object */
    bool mKspIsSetup; /**< Used by Solve method to track whether KSP has been used. */
    double mNonZerosUsed;  /**< Yes, it really is stored as a double. */
    bool mMatrixIsConstant; /**< Whether the matrix is unchanged each time Solve() is called */
    double mTolerance; /**< absolute or relative tolerance of the KSP solver */
    /**
     * Sets either absolute or relative tolerance of the KSP solver.
     * Default is to false
     */
    bool mUseAbsoluteTolerance;
    std::string mKspType;/**< KSP solver type (see PETSc KSPSetType() ) */
    std::string mPcType;/**< Preconditioner type (see PETSc PCSetType() ) */

    Vec mDirichletBoundaryConditionsVector; /**< Storage for efficient application of Dirichlet BCs, see AbstractBoundaryConditionsContainer */

    /// \todo: #1082 Create an abstract class for the preconditioners and use a single pointer
    /** Stores a pointer to a purpose-build preconditioner*/
    PCBlockDiagonal* mpBlockDiagonalPC;
    /** Stores a pointer to a purpose-build preconditioner*/
    PCLDUFactorisation* mpLDUFactorisationPC;
    /** Stores a pointer to a purpose-build preconditioner*/
    PCTwoLevelsBlockDiagonal* mpTwoLevelsBlockDiagonalPC;

    /** Pointer to vector containing a list of bath nodes*/
    boost::shared_ptr<std::vector<PetscInt> > mpBathNodes;

    /** Whether the matrix used for preconditioning is the same as the LHS*/
    bool mPrecondMatrixIsNotLhs;

    /** The max number of nonzero entries expected on a LHS row */
    unsigned mRowPreallocation;

    /** Whether to use fixed number of iterations */
    bool mUseFixedNumberIterations;

    /**
     * When using fixed number of iterations, a solve with residual-based
     * stop criteria will be performed every mEvaluateNumItsEveryNSolves solves
     * to decide how many iterations perform for the next mEvaluateNumItsEveryNSolves-1 solves
     */
    unsigned mEvaluateNumItsEveryNSolves;

    /**  Context for KSPDefaultConverged() */
    void* mpConvergenceTestContext;

    /** Number of solves performed since the current object was created */
    unsigned mNumSolves;

    /** Preconditioned operator smallest eigenvalue */
    PetscReal mEigMin;

    /** Preconditioned operator largest eigenvalue */
    PetscReal mEigMax;

    /** Under certain circunstances you have to reevaluate the spectrum before the k*n-th, k=0,1,..., iteration*/
    bool mForceSpectrumReevaluation;

#ifdef TRACE_KSP
    unsigned mTotalNumIterations;
    unsigned mMaxNumIterations;
#endif
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        //MatNullSpace mMatNullSpace; // Gets re-created by calling code on load

        archive & mNonZerosUsed;
        archive & mMatrixIsConstant;
        archive & mTolerance;
        archive & mUseAbsoluteTolerance;
        archive & mKspType;
        archive & mPcType;

        //Vec mDirichletBoundaryConditionsVector; // Gets re-created by calling code on load
    }

public:

    /**
     * Constructor.
     *
     * @param lhsVectorSize  the size of the LHS vector
     * @param rowPreallocation the max number of nonzero entries expected on a row
     *  - a value of 0 is allowed: no preallocation is then done and the user must
     *    preallocate the memory for the matrix themselves.
     *  - the default value allows for small size systems to be set as dense matrices
     *    automatically.
     */
    LinearSystem(PetscInt lhsVectorSize, unsigned rowPreallocation=UINT_MAX);

    /**
     * Alternative constructor.
     *
     * Create a linear system, where the size is based on the size of a given
     * PETSc vec.
     *
     * The LHS & RHS vectors will be created by duplicating this vector's
     * settings.  This should avoid problems with using VecScatter on
     * bidomain simulation results.
     *
     * @param templateVector  a PETSc vec
     * @param rowPreallocation the max number of nonzero entries expected on a row
     * @param newAllocationError tells PETSc whether to set the MAT_NEW_NONZERO_ALLOCATION_ERR.
     *        ** currently only used in PETSc 3.3 and later **
     *        in PETSc 3.2 and earlier MAT_NEW_NONZERO_ALLOCATION_ERR defaults to false
     *        in PETSc 3.3 MAT_NEW_NONZERO_ALLOCATION_ERR defaults to true
     */
    LinearSystem(Vec templateVector, unsigned rowPreallocation, bool newAllocationError=true);

    /**
     * Alternative constructor.
     *
     * Create a linear system which wraps the provided PETSc objects so we can
     * access them using our API.  Either of the objects may be NULL, but at
     * least one of them must not be.
     *
     * Useful for storing residuals and jacobians when solving nonlinear PDEs.
     *
     * @param residualVector  the residual vector
     * @param jacobianMatrix  the Jacobian matrix
     */
    LinearSystem(Vec residualVector, Mat jacobianMatrix);

    /**
     * Alternative constructor for archiving.
     *
     * @param lhsVectorSize  the size of the LHS vector
     * @param lhsMatrix  the RHS matrix
     * @param rhsVector  the RHS vector
     */
    LinearSystem(PetscInt lhsVectorSize, Mat lhsMatrix, Vec rhsVector);

    /**
     * Destructor.
     */
    ~LinearSystem();

//    bool IsMatrixEqualTo(Mat testMatrix);
//    bool IsRhsVectorEqualTo(Vec testVector);

    /**
     * Change one of the entires of the matrix to the specified value.
     *
     * @param row  the row index
     * @param col  the column index
     * @param value  the value for this entry
     */
    void SetMatrixElement(PetscInt row, PetscInt col, double value);

    /**
     * Add the specified value to an entry of the matrix.
     *
     * @param row  the row index
     * @param col  the column index
     * @param value  the value for this entry
     */
    void AddToMatrixElement(PetscInt row, PetscInt col, double value);

    /**
     * Call this before Solve().
     *
     * This calls FinaliseLhsMatrix() and FinaliseRhsVector().
     */
    void AssembleFinalLinearSystem();

    /**
     * Should be called before AddToMatrixElement.
     *
     * This calls SwitchWriteModeLhsMatrix() and FinaliseRhsVector().
     */
    void AssembleIntermediateLinearSystem();

    /**
     * Sets up the PETSc matrix left-hand-side mLhsMatrix
     */
    void FinaliseLhsMatrix();

    /**
     * Sets up the PETSc matrix used for preconditioning.
     */
    void FinalisePrecondMatrix();

    /**
     * Sets up the PETSc matrix left-hand-side mLhsMatrix
     */
    void SwitchWriteModeLhsMatrix();

    /**
     * Sets up the PETSc vector right-hand-side mRhsVector
     */
    void FinaliseRhsVector();

    /**
     * Force PETSc to treat the matrix in this linear system as symmetric from now on.
     *
     * @param isSymmetric  whether the matrix is symmetric or not (defaults to true)
     */
    void SetMatrixIsSymmetric(bool isSymmetric=true);

    /**
     * Get whether PETSc considers the matrix in this linear system as symmetric or not.
     *
     * @return whether the matrix is symmetric or not.
     */
    bool IsMatrixSymmetric();

    /**
     * Set mMatrixIsConstant.
     *
     * @param matrixIsConstant  whether the matrix is constant
     */
    void SetMatrixIsConstant(bool matrixIsConstant);

    /**
     * Set the relative tolerance.
     *
     * @param relativeTolerance  the relative tolerance
     */
    void SetRelativeTolerance(double relativeTolerance);

    /**
     * Set the absolute tolerance.
     *
     * @param absoluteTolerance  the absolute tolerance
     */
    void SetAbsoluteTolerance(double absoluteTolerance);

    /**
     * Set the KSP solver type (see PETSc KSPSetType() for valid arguments).
     *
     * @param kspType  the KSP solver type
     */
    void SetKspType(const char* kspType);

    /**
     * Set the preconditioner type  (see PETSc PCSetType() for valid arguments).
     *
     * @param pcType the preconditioner type
     * @param pBathNodes the list of nodes defining the bath
     */
    /// \todo: #1082 is this the way of defining a null pointer as the default value of pBathNodes?
    void SetPcType(const char* pcType, boost::shared_ptr<std::vector<PetscInt> > pBathNodes=boost::shared_ptr<std::vector<PetscInt> >() );

    /**
     * Display the left-hand side matrix.
     */
    void DisplayMatrix();

    /**
     * Display the right-hand side vector.
     */
    void DisplayRhs();

    /**
     * Set all entries in a given row of a matrix to a certain value.
     * This must be called by the process who owns the row, (but other
     * processors will treat it as a null-op)
     *
     * @param row  the row index
     * @param value  the value to set each entry in this row
     */
    void SetMatrixRow(PetscInt row, double value);

    /**
     * Returns the i-th row of the LHS matrix as a distributed PETSc Vec
     *
     * @param rowIndex the row index
     * @return rowIndex-th row of the matrix in distributed format
     */
    Vec GetMatrixRowDistributed(unsigned rowIndex);

    /**
     * Zero several rows of the matrix, putting a given value in the diagonal entries.
     *
     * *Massively* less expensive than zeroing each matrix row individually
     *
     * @param rRows std::vector of rows to be zeroed
     * @param diagonalValue value to put in the diagonal entries (of the zeroed rows)
     */
    void ZeroMatrixRowsWithValueOnDiagonal(std::vector<unsigned>& rRows, double diagonalValue);


    /**
     * Zero several rows and columns of the matrix, putting a given value on the diagonal.
     *
     * @param rRowColIndices A list of indices. All the rows with these indices, and all the columns
     * with these indices, will be zeroed.
     * @param diagonalValue value to put in the diagonal entries (of the zeroed rows)
     */
    void ZeroMatrixRowsAndColumnsWithValueOnDiagonal(std::vector<unsigned>& rRowColIndices, double diagonalValue);


    /**
     * Zero a column of the left-hand side matrix.
     *
     * Unfortunately there is no equivalent method in Petsc, so this has to be
     * done carefully to ensure that the sparsity structure of the matrix
     * is not broken. Only owned entries which are non-zero are zeroed.
     *
     * @param col  the column index
     */
    void ZeroMatrixColumn(PetscInt col);

    /**
     * Zero all entries of the left-hand side matrix.
     */
    void ZeroLhsMatrix();

    /**
     * Zero all entries of the right-hand side vector.
     */
    void ZeroRhsVector();

    /**
     * Zero all entries of the left-hand side matrix and right-hand side vector.
     */
    void ZeroLinearSystem();

    /**
     * Solve the linear system.
     * @return the solution vector
     * @param lhsGuess  an optional initial guess for the solution (defaults to NULL)
     */
    Vec Solve(Vec lhsGuess=nullptr);

    /**
     * Set an element of the right-hand side vector to a given value.
     *
     * @param row  the row index
     * @param value  the value to set this entry
     */
    void SetRhsVectorElement(PetscInt row, double value);

    /**
     * Add a value to an element of the right-hand side vector.
     *
     * @param row  the row index
     * @param value  the value to set this entry
     */
    void AddToRhsVectorElement(PetscInt row, double value);

    /**
     * @return #mSize.
     */
    unsigned GetSize() const;

    /**
     * Set the null basis of the linear system.
     * In debug mode we test for orthonormality and throw EXCEPTIONs:
     *  - each direction is of unit length (to within a tolerance)
     *  - the directions are mutually orthogonal (to within a tolerance)
     * Note that in NDEBUG (optimised mode) these tests are skipped and the parameters
     * are passed directly into the LinearSystem.
     *
     * @param nullbasis  an array PETSc vectors containing orthogonal directions in the nullspace
     * @param numberOfBases the number of directions (size of nullbasis array)
     */
    void SetNullBasis(Vec nullbasis[], unsigned numberOfBases);

    /**
     * Remove the null space from the linear system.
     *
     * Use for example if Dirichlet BC are applied to a singular system and, therefore, there's
     * no null space anymore.
     */
    void RemoveNullSpace();

    /**
     * @return access to the RHS vector directly. Shouldn't generally need to be called.
     */
    Vec& rGetRhsVector();

    /**
     * @return access to the RHS vector for archiving
     */
    Vec GetRhsVector() const;

    /**
     * @return access to the LHS matrix directly. Shouldn't generally need to be called.
     */
    Mat& rGetLhsMatrix();

    /**
     * @return access to the LHS matrix for archiving
     */
    Mat GetLhsMatrix() const;

    /**
     * @return access to the matrix used for preconditioning.
     */
    Mat& rGetPrecondMatrix();

    /**
     * @return access to the dirichlet boundary conditions vector.
     *
     * Should only be used by the BoundaryConditionsContainer.
     */
    Vec& rGetDirichletBoundaryConditionsVector();

    // DEBUGGING CODE:
    /**
     * @return this process's ownership range of the contents of the system.
     *
     * @param lo  lowest index owned by this process
     * @param hi  highest index owned by this process
     */
    void GetOwnershipRange(PetscInt& lo, PetscInt& hi);

    /**
     * @return an element of the matrix.
     * May only be called for elements you own.
     *
     * @param row  the row index
     * @param col  the column index
     */
    double GetMatrixElement(PetscInt row, PetscInt col);

    /**
     * @return an element of the RHS vector.
     * May only be called for elements you own.
     *
     * @param row  the row index
     */
    double GetRhsVectorElement(PetscInt row);

    /**
     * @return the number of iterations taken by the last Solve()
     */
    unsigned GetNumIterations() const;

    /**
     * Add multiple values to the matrix of linear system.
     *
     * @param matrixRowAndColIndices mapping from index of the ublas matrix (see param below)
     *  to index of the PETSc matrix of this linear system
     * @param rSmallMatrix Ublas matrix containing the values to be added
     *
     * N.B. Values which are not local (ie the row is not owned) will be skipped.
     */
    template<size_t MATRIX_SIZE>
    void AddLhsMultipleValues(unsigned* matrixRowAndColIndices, c_matrix<double, MATRIX_SIZE, MATRIX_SIZE>& rSmallMatrix)
    {
        PetscMatTools::AddMultipleValues(mLhsMatrix, matrixRowAndColIndices, rSmallMatrix);
    }

    /**
     * Add multiple values to the RHS vector.
     *
     * @param vectorIndices mapping from index of the ublas vector (see param below)
     *  to index of the vector of this linear system
     * @param smallVector Ublas vector containing the values to be added
     *
     * N.B. Values which are not local (ie the row is not owned) will be skipped.
     */
    template<size_t VECTOR_SIZE>
    void AddRhsMultipleValues(unsigned* vectorIndices, c_vector<double, VECTOR_SIZE>& smallVector)
    {
        PetscVecTools::AddMultipleValues(mRhsVector, vectorIndices, smallVector);
    }

    /**
     * Set method for #mPrecondMatrixIsNotLhs.
     * @param precondIsDifferent  whether the matrix used for preconditioning is the same as the LHS.
     */
    void SetPrecondMatrixIsDifferentFromLhs(bool precondIsDifferent = true);

    /**
     * Set method for #mUseFixedNumberIterations
     * @param useFixedNumberIterations whether to use fixed number of iterations
     * @param evaluateNumItsEveryNSolves tells LinearSystem to perform a solve with convergence-based stop criteria every n solves to decide how many iterations perform for the next n-1 solves. Default is perfoming a single evaluation at the beginning of the simulation.
     */
    void SetUseFixedNumberIterations(bool useFixedNumberIterations = true, unsigned evaluateNumItsEveryNSolves = UINT_MAX);

    /**
     * Method to regenerate all KSP objects, including the solver and the preconditioner (e.g. after
     * changing the PDE time step when using time adaptivity).
     */
    void ResetKspSolver();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(LinearSystem)

namespace boost
{
namespace serialization
{
template<class Archive>
inline void save_construct_data(
    Archive & ar, const LinearSystem * t, const unsigned int file_version)
{

    std::string archive_filename_lhs = ArchiveLocationInfo::GetArchiveDirectory() + "lhs.mat";
    std::string archive_filename_rhs = ArchiveLocationInfo::GetArchiveDirectory() + "rhs.vec";
    const unsigned size = t->GetSize();
    ar << size;

    PetscTools::DumpPetscObject(t->GetRhsVector(), archive_filename_rhs);
    PetscTools::DumpPetscObject(t->GetLhsMatrix(), archive_filename_lhs);

    //Is the matrix structurally symmetric?
    PetscBool symm_set, is_symmetric;
    is_symmetric = PETSC_FALSE;
    //Note that the following call only changes is_symmetric when symm_set is true
    MatIsSymmetricKnown(t->GetLhsMatrix(), &symm_set, &is_symmetric);
    assert(symm_set == is_symmetric);
    ar << symm_set;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor)
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, LinearSystem * t, const unsigned int file_version)
{
     std::string archive_filename_lhs = ArchiveLocationInfo::GetArchiveDirectory() + "lhs.mat";
     std::string archive_filename_rhs = ArchiveLocationInfo::GetArchiveDirectory() + "rhs.vec";

     PetscInt size;
     ar >> size;

     Vec new_vec;
     PetscTools::ReadPetscObject(new_vec, archive_filename_rhs);

     Mat new_mat;
     PetscTools::ReadPetscObject(new_mat, archive_filename_lhs);

     //This has to occur after the call to MatLoad as the matrix does not exist until MatLoad is called.
     //The property will be copied & set correctly in the LinearSystem constructor.
     PetscBool symm_set;

     ar >> symm_set;
     if (symm_set == PETSC_TRUE)
     {
        PetscMatTools::SetOption(new_mat, MAT_SYMMETRIC);
        PetscMatTools::SetOption(new_mat, MAT_SYMMETRY_ETERNAL);
     }

     ::new(t)LinearSystem(size, new_mat, new_vec);
}
}
} // namespace ...

#endif //_LINEARSYSTEM_HPP_
