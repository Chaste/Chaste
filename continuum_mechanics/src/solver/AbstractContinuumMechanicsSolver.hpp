/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef ABSTRACTCONTINUUMMECHANICSSOLVER_HPP_
#define ABSTRACTCONTINUUMMECHANICSSOLVER_HPP_

#include "ContinuumMechanicsProblemDefinition.hpp"
#include "CompressibilityType.hpp"
#include "LinearSystem.hpp"
#include "OutputFileHandler.hpp"
#include "ReplicatableVector.hpp"
#include "QuadraticMesh.hpp"
#include "Warnings.hpp"
#include "PetscException.hpp"
#include "GaussianQuadratureRule.hpp"
#include "PetscTools.hpp"
#include "MechanicsEventHandler.hpp"

/**
 *  Simple enumeration for options that can be passed into
 *  AbstractContinuumMechanicsSolver::ApplyDirichletBoundaryConditions().
 *  See documentation for this method.
 */
typedef enum _ApplyDirichletBcsType
{
    LINEAR_PROBLEM,
    NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY,
    NONLINEAR_PROBLEM_APPLY_TO_EVERYTHING,
} ApplyDirichletBcsType;


/**
 *  General base class for continuum mechanics solvers
 */
template<unsigned DIM>
class AbstractContinuumMechanicsSolver
{
protected:
    /**
     * The mesh to be solved on. Requires 6 nodes per triangle (or 10 per tetrahedron)
     * as quadratic bases are used.
     */
    QuadraticMesh<DIM>& mrQuadMesh;

    /** Problem definition class - contains info on boundary conditions, etc */
    ContinuumMechanicsProblemDefinition<DIM>& mrProblemDefinition;

    /** Whether to write any output. */
    bool mWriteOutput;

    /** Where to write output, relative to CHASTE_TESTOUTPUT. */
    std::string mOutputDirectory;

    /** Output file handler. */
    OutputFileHandler* mpOutputFileHandler;

    /** Spatial solution:
     *  For solids problems, mSpatialSolution[i](j) = x_j (new position) for node i.
     *  For fluids problems, mSpatialSolution[i](j) = u_j (flow) for node i. */
    std::vector<c_vector<double,DIM> > mSpatialSolution;

    /** Pressures solution at each vertex of the mesh. Only valid if mCompressibilityType==INCOMPRESSIBLE. */
    std::vector<double> mPressureSolution;


    /**
     * The current solution, in the form (assuming 2d):
     *   Incompressible problem: [u1 v1 u2 v2 ... uN vN p1 p2 .. pM]
     *   Compressible problem:   [u1 v1 u2 v2 ... uN vN]
     * where there are N total nodes and M vertices.
     */
    std::vector<double> mCurrentSolution;

    /** Gaussian quadrature rule. */
    GaussianQuadratureRule<DIM>* mpQuadratureRule;

    /** Boundary Gaussian quadrature rule. */
    GaussianQuadratureRule<DIM-1>* mpBoundaryQuadratureRule;

    /**
     * This is equal to either COMPRESSIBLE or INCOMPRESSIBLE (see enumeration defined at top of file)
     * and is only used in computing mNumDofs and allocating matrix memory.
     */
    CompressibilityType mCompressibilityType;

    /**
     * Number of degrees of freedom (equal to, in the incompressible case:
     * DIM*N + M if quadratic-linear bases are used, where there are N total
     * nodes and M vertices; or DIM*N in the compressible case).
     */
    unsigned mNumDofs;

    /**
     * Residual vector nonlinear problems.
     *
     * Since the residual in nonlinear problems is usually also the RHS vector in the linear
     * system, it may seem unncessary to also have the member variable mLinearSystemRhsVector.
     *
     * However: Newton's method is Ju = f, where J is the Jacobian, u the (negative of the) update
     * and f the residual, but when applying Dirichlet boundary conditions in
     * the compressible case, we alter the rows of the matrix and also alter the columns in order to
     * maintain symmetry. This requires making further changes to the right-hand vector, meaning that
     * it no longer properly represents the residual. Hence, we have to use two vectors.
     *
     * Overall, this can be represents as
     *  - compute residual f
     *  - compute Jacobian J
     *  - apply BCs to f.
     *  - alter the linear system from Ju=f to (J*)u=f* which enforces the dirichlet boundary conditions but enforces them symmetrically.
     *
     * mLinearSystemRhsVector below represents f*.
     */
    Vec mResidualVector;

    /**
     * The RHS side in the linear system that is solved each Newton iteration.
     */
    Vec mLinearSystemRhsVector;

    /**
     * Jacobian matrix of the nonlinear system, LHS matrix for the linear system.
     */
    Mat mSystemLhsMatrix;

    /**
     * Helper vector (see ApplyDirichletBoundaryConditions code).
     */
    Vec mDirichletBoundaryConditionsVector;

    /**
     * Precondition matrix for the linear system.
     */
    Mat mPreconditionMatrix;


    /**
     * Allocates memory for the matrices and vectors
     */
    void AllocateMatrixMemory();


    /**
     * Apply the Dirichlet boundary conditions to the linear system.
     *
     * The first input parameter should be one of the following
     *   LINEAR_PROBLEM -- indicating the overall problem is linear, and in which
     *     case the BCs will be applied to both the matrix and vector of the
     *     linear system
     *   NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY -- indicating the overall problem is nonlinear
     *     and here only the residual vector will be altered (apply the Dirichlet boundary
     *     conditions to the residual vector involves setting appropriate components to the
     *     difference between the current value and the correct value).
     *
     *   NONLINEAR_PROBLEM_APPLY_TO_EVERYTHING -- indicating the overall problem is nonlinear,
     *     and here, the residual vector will be altered, as will the matrix and RHS vector
     *     (see documentation on mResidualVector for why there is a separate residual vector
     *     and RHS vector).
     *
     *  The second parameter should be true if the overall problem is symmetric, in which case
     *  boundary conditions will be applied keeping the matrix symmetric (if the matrix is being
     *  altered). See in-code comments for how this is done.
     *
     *  @param type see above
     *  @param symmetricProblem see above
     */
    void ApplyDirichletBoundaryConditions(ApplyDirichletBcsType type, bool symmetricProblem);


public:
    /**
     *  Constructor
     *  @param rQuadMesh the mesh
     *  @param rProblemDefinition problem definition object
     *  @param outputDirectory output directory name
     *  @param compressibilityType 'INCOMPRESSIBLE' or 'COMPRESSIBLE'
     */
    AbstractContinuumMechanicsSolver(QuadraticMesh<DIM>& rQuadMesh,
                                     ContinuumMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                     std::string outputDirectory,
                                     CompressibilityType compressibilityType);

    /** Destructor */
    virtual ~AbstractContinuumMechanicsSolver();

    /**
     * Write the spatial solution (deformed position if solids, flow if fluids) at the nodes
     *
     * @param fileName (stem)
     * @param fileExtension to append at end.
     * @param counterToAppend append a counter to the file name (defaults to nothing appended).
     *
     * For example:
     * WriteCurrentSpatialSolution("solution","nodes") --> file called "solution.nodes"
     * WriteCurrentSpatialSolution("solution","nodes",3) --> file called "solution_3.nodes"
     */
    void WriteCurrentSpatialSolution(std::string fileName, std::string fileExtension, int counterToAppend=-1);

    /**
     * Write the pressure solution. Only valid if mCompressibilityType==INCOMPRESSIBLE.
     *
     * @param counterToAppend append a counter to the file name
     *
     * WriteCurrentPressureSolution() --> file called "pressure.txt"
     * WriteCurrentPressureSolution(3) --> file called "pressure_3.txt"
     */
    void WriteCurrentPressureSolution( int counterToAppend=-1);

    /**
     * Set whether to write any output.
     *
     * @param writeOutput (defaults to true)
     */
    void SetWriteOutput(bool writeOutput=true);

    /**
     * Get the current solution vector (advanced use only - for getting the deformed position use
     * rGetDeformedPosition()).
     */
    std::vector<double>& rGetCurrentSolution()
    {
        return mCurrentSolution;
    }

    /**
     * Get the spatial solution. For solids problems this will be the deformed position,
     * for fluids problems this will be the flow.
     */
    virtual std::vector<c_vector<double,DIM> >& rGetSpatialSolution()=0;

    /**
     *  Get the pressure, for each vertex in the mesh. Only valid if mCompressibilityType==INCOMPRESSIBLE
     */
    std::vector<double>& rGetPressures();

};


template<unsigned DIM>
AbstractContinuumMechanicsSolver<DIM>::AbstractContinuumMechanicsSolver(QuadraticMesh<DIM>& rQuadMesh,
                                                                        ContinuumMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                                                        std::string outputDirectory,
                                                                        CompressibilityType compressibilityType)
    : mrQuadMesh(rQuadMesh),
      mrProblemDefinition(rProblemDefinition),
      mOutputDirectory(outputDirectory),
      mpOutputFileHandler(NULL),
      mpQuadratureRule(NULL),
      mpBoundaryQuadratureRule(NULL),
      mCompressibilityType(compressibilityType),
      mResidualVector(NULL),
      mSystemLhsMatrix(NULL),
      mPreconditionMatrix(NULL)
{
    assert(DIM==2 || DIM==3);

    mWriteOutput = (mOutputDirectory != "");
    if (mWriteOutput)
    {
        mpOutputFileHandler = new OutputFileHandler(mOutputDirectory);
    }

    if (mCompressibilityType==COMPRESSIBLE)
    {
        mNumDofs = DIM*mrQuadMesh.GetNumNodes();
    }
    else
    {
        mNumDofs = DIM*mrQuadMesh.GetNumNodes() + mrQuadMesh.GetNumVertices();
    }

    AllocateMatrixMemory();

    mpQuadratureRule = new GaussianQuadratureRule<DIM>(3);
    mpBoundaryQuadratureRule = new GaussianQuadratureRule<DIM-1>(3);

    mCurrentSolution.resize(mNumDofs, 0.0);
}


template<unsigned DIM>
AbstractContinuumMechanicsSolver<DIM>::~AbstractContinuumMechanicsSolver()
{
    if (mpOutputFileHandler)
    {
        delete mpOutputFileHandler;
    }

    if (mpQuadratureRule)
    {
        delete mpQuadratureRule;
        delete mpBoundaryQuadratureRule;
    }

    if (mResidualVector)
    {
        PetscTools::Destroy(mResidualVector);
        PetscTools::Destroy(mLinearSystemRhsVector);
        PetscTools::Destroy(mSystemLhsMatrix);
        PetscTools::Destroy(mPreconditionMatrix);
    }
    if (mDirichletBoundaryConditionsVector)
    {
        PetscTools::Destroy(mDirichletBoundaryConditionsVector);
    }
}

template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::WriteCurrentSpatialSolution(std::string fileName,
                                                                        std::string fileExtension,
                                                                        int counterToAppend)
{
    // Only write output if the flag mWriteOutput has been set
    if (!mWriteOutput)
    {
        return;
    }

    std::stringstream file_name;
    file_name << fileName;
    if (counterToAppend >= 0)
    {
        file_name << "_" << counterToAppend;
    }
    file_name << "." << fileExtension;

    out_stream p_file = mpOutputFileHandler->OpenOutputFile(file_name.str());

    std::vector<c_vector<double,DIM> >& r_spatial_solution = rGetSpatialSolution();
    for (unsigned i=0; i<r_spatial_solution.size(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            *p_file << r_spatial_solution[i](j) << " ";
        }
        *p_file << "\n";
    }
    p_file->close();
}

template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::WriteCurrentPressureSolution(int counterToAppend)
{
    // Only write output if the flag mWriteOutput has been set
    if (!mWriteOutput)
    {
        return;
    }

    std::stringstream file_name;
    file_name << "pressure";
    if (counterToAppend >= 0)
    {
        file_name << "_" << counterToAppend;
    }
    file_name << ".txt";

    out_stream p_file = mpOutputFileHandler->OpenOutputFile(file_name.str());

    std::vector<double>& r_pressure = rGetPressures();
    for (unsigned i=0; i<r_pressure.size(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            *p_file << mrQuadMesh.GetNode(i)->rGetLocation()[j] << " ";
        }

        *p_file << r_pressure[i] << "\n";
    }
    p_file->close();
}

template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::SetWriteOutput(bool writeOutput)
{
    if (writeOutput && (mOutputDirectory==""))
    {
        EXCEPTION("Can't write output if no output directory was given in constructor");
    }
    mWriteOutput = writeOutput;
}


template<unsigned DIM>
std::vector<double>& AbstractContinuumMechanicsSolver<DIM>::rGetPressures()
{
    mPressureSolution.clear();
    mPressureSolution.resize(mrQuadMesh.GetNumVertices());

    for (unsigned i=0; i<mrQuadMesh.GetNumVertices(); i++)
    {
        mPressureSolution[i] = mCurrentSolution[DIM*mrQuadMesh.GetNumNodes() + i];
    }
    return mPressureSolution;
}

/*
 * This method applies the appropriate BCs to the residual and/or linear system
 *
 * For the latter, and second input parameter==true, the BCs are imposed in such a way as
 * to ensure that a symmetric linear system remains symmetric. For each row with boundary condition
 * applied, both the row and column are zero'd and the RHS vector modified to take into account the
 * zero'd column.
 *
 * Suppose we have a matrix
 * [a b c] [x] = [ b1 ]
 * [d e f] [y]   [ b2 ]
 * [g h i] [z]   [ b3 ]
 * and we want to apply the boundary condition x=v without losing symmetry if the matrix is
 * symmetric. We apply the boundary condition
 * [1 0 0] [x] = [ v  ]
 * [d e f] [y]   [ b2 ]
 * [g h i] [z]   [ b3 ]
 * and then zero the column as well, adding a term to the RHS to take account for the
 * zero-matrix components
 * [1 0 0] [x] = [ v  ] - v[ 0 ]
 * [0 e f] [y]   [ b2 ]    [ d ]
 * [0 h i] [z]   [ b3 ]    [ g ]
 * Note the last term is the first column of the matrix, with one component zeroed, and
 * multiplied by the boundary condition value. This last term is then stored in
 * mDirichletBoundaryConditionsVector, and in general is equal to:
 * SUM_{d=1..D} v_d a'_d
 * where v_d is the boundary value of boundary condition d (d an index into the matrix),
 * and a'_d is the dth-column of the matrix but with the d-th component zeroed, and where
 * there are D boundary conditions
 */
template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::ApplyDirichletBoundaryConditions(ApplyDirichletBcsType type, bool symmetricProblem)
{
    std::vector<unsigned> rows;
    std::vector<double> values;

    // Whether to apply symmetrically, ie alter columns as well as rows (see comment above)
    bool applySymmetrically = (type!=NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY) && symmetricProblem;

    if (applySymmetrically)
    {
        if(mDirichletBoundaryConditionsVector==NULL)
        {
            VecDuplicate(mResidualVector, &mDirichletBoundaryConditionsVector);
        }

        PetscVecTools::Zero(mDirichletBoundaryConditionsVector);
        PetscMatTools::Finalise(mSystemLhsMatrix);
    }

    ///////////////////////////////////////
    // collect the entries to be altered
    ///////////////////////////////////////

    for (unsigned i=0; i<mrProblemDefinition.rGetDirichletNodes().size(); i++)
    {
        unsigned node_index = mrProblemDefinition.rGetDirichletNodes()[i];

        for (unsigned j=0; j<DIM; j++)
        {
            double dirichlet_val = mrProblemDefinition.rGetDirichletNodeValues()[i](j); // problem defn returns DISPLACEMENTS here

            if(dirichlet_val != ContinuumMechanicsProblemDefinition<DIM>::FREE)
            {
                double val;
                unsigned dof_index = DIM*node_index+j;

                if(type == LINEAR_PROBLEM)
                {
                    val = dirichlet_val;
                }
                else
                {
                    // The boundary conditions on the NONLINEAR SYSTEM are x=boundary_values
                    // on the boundary nodes. However:
                    // The boundary conditions on the LINEAR SYSTEM at each Newton step (Ju=f,
                    // where J is the Jacobian, u the negative update vector and f is the residual) is:
                    // u=current_soln-boundary_values on the boundary nodes
                    val = mCurrentSolution[dof_index] - dirichlet_val;
                }
                rows.push_back(dof_index);
                values.push_back(val);
            }
        }
    }

    ///////////////////////////////////////
    // do the alterations
    ///////////////////////////////////////

    if (applySymmetrically)
    {
        // Modify the matrix columns
        for (unsigned i=0; i<rows.size(); i++)
        {
            unsigned col = rows[i];
            double minus_value = -values[i];

            // Get a vector which will store the column of the matrix (column d, where d is
            // the index of the row (and column) to be altered for the boundary condition.
            // Since the matrix is symmetric when get row number "col" and treat it as a column.
            // PETSc uses compressed row format and therefore getting rows is far more efficient
            // than getting columns.
            Vec matrix_col = PetscMatTools::GetMatrixRowDistributed(mSystemLhsMatrix,col);

            // Zero the correct entry of the column
            PetscVecTools::SetElement(matrix_col, col, 0.0);

            // Set up the RHS Dirichlet boundary conditions vector
            // Eg, for a boundary node at the zeroth node (x_0 = value), this is equal to
            //   -value*[0 a_21 a_31 .. a_N1]
            // and will be added to the RHS.
            PetscVecTools::AddScaledVector(mDirichletBoundaryConditionsVector, matrix_col, minus_value);
            PetscTools::Destroy(matrix_col);
        }
    }

    if (type!=NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY) // ie doing a whole linear system
    {
        // Now zero the appropriate rows and columns of the matrix. If the matrix is symmetric we apply the
        // boundary conditions in a way the symmetry isn't lost (rows and columns). If not only the row is
        // zeroed.
        if (applySymmetrically)
        {
            PetscMatTools::ZeroRowsAndColumnsWithValueOnDiagonal(mSystemLhsMatrix, rows, 1.0);
            PetscMatTools::ZeroRowsAndColumnsWithValueOnDiagonal(mPreconditionMatrix, rows, 1.0);

            // Apply the RHS boundary conditions modification if required.
            PetscVecTools::AddScaledVector(mLinearSystemRhsVector, mDirichletBoundaryConditionsVector, 1.0);
        }
        else
        {
            PetscMatTools::ZeroRowsWithValueOnDiagonal(mSystemLhsMatrix, rows, 1.0);
            PetscMatTools::ZeroRowsWithValueOnDiagonal(mPreconditionMatrix, rows, 1.0);
        }
    }

    if (type!=LINEAR_PROBLEM)
    {
        for (unsigned i=0; i<rows.size(); i++)
        {
            PetscVecTools::SetElement(mResidualVector, rows[i], values[i]);
        }
    }

    if (type!=NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY)
    {
        for (unsigned i=0; i<rows.size(); i++)
        {
            PetscVecTools::SetElement(mLinearSystemRhsVector, rows[i], values[i]);
        }
    }
}



template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::AllocateMatrixMemory()
{
    ///////////////////////////
    // three vectors
    ///////////////////////////
    mResidualVector = PetscTools::CreateVec(mNumDofs);
    VecDuplicate(mResidualVector, &mLinearSystemRhsVector);
    // the one is only allocated if it will be needed (in ApplyDirichletBoundaryConditions),
    // depending on whether the matrix is kept symmetric.
    mDirichletBoundaryConditionsVector = NULL;

    ///////////////////////////
    // two matrices
    ///////////////////////////

    if (DIM==2)
    {
        // 2D: N elements around a point => 7N+3 non-zeros in that row? Assume N<=10 (structured mesh would have N_max=6) => 73.
        unsigned num_non_zeros = std::min(75u, mNumDofs);

        PetscTools::SetupMat(mSystemLhsMatrix, mNumDofs, mNumDofs, num_non_zeros, PETSC_DECIDE, PETSC_DECIDE);
        PetscTools::SetupMat(mPreconditionMatrix, mNumDofs, mNumDofs, num_non_zeros, PETSC_DECIDE, PETSC_DECIDE);
    }
    else
    {
        assert(DIM==3);

        // in 3d we get the number of containing elements for each node and use that to obtain an upper bound
        // for the number of non-zeros for each DOF associated with that node.

        int* num_non_zeros_each_row = new int[mNumDofs];
        for (unsigned i=0; i<mNumDofs; i++)
        {
            num_non_zeros_each_row[i] = 0;
        }

        for (unsigned i=0; i<mrQuadMesh.GetNumNodes(); i++)
        {
            // this upper bound neglects the fact that two containing elements will share the same nodes..
            // 4 = max num dofs associated with this node
            // 30 = 3*9+3 = 3 dimensions x 9 other nodes on this element   +  3 vertices with a pressure unknown
            unsigned num_non_zeros_upper_bound = 4 + 30*mrQuadMesh.GetNode(i)->GetNumContainingElements();

            num_non_zeros_upper_bound = std::min(num_non_zeros_upper_bound, mNumDofs);

            num_non_zeros_each_row[DIM*i + 0] = num_non_zeros_upper_bound;
            num_non_zeros_each_row[DIM*i + 1] = num_non_zeros_upper_bound;
            num_non_zeros_each_row[DIM*i + 2] = num_non_zeros_upper_bound;

            if (mCompressibilityType==INCOMPRESSIBLE)
            {
                //Could do !mrQuadMesh.GetNode(i)->IsInternal()
                if (i<mrQuadMesh.GetNumVertices()) // then this is a vertex
                {
                    num_non_zeros_each_row[DIM*mrQuadMesh.GetNumNodes() + i] = num_non_zeros_upper_bound;
                }
            }
        }

        // NOTE: PetscTools::SetupMat() or the below creates a MATAIJ matrix, which means the matrix will
        // be of type MATSEQAIJ if num_procs=1 and MATMPIAIJ otherwise. In the former case
        // MatSeqAIJSetPreallocation MUST be called [MatMPIAIJSetPreallocation will have
        // no effect (silently)], and vice versa in the latter case

        /// We want to allocate different numbers of non-zeros per row, which means
        /// PetscTools::SetupMat isn't that useful. We could call
        //PetscTools::SetupMat(mSystemLhsMatrix, mNumDofs, mNumDofs, 0, PETSC_DECIDE, PETSC_DECIDE);
        //PetscTools::SetupMat(mPreconditionMatrix, mNumDofs, mNumDofs, 0, PETSC_DECIDE, PETSC_DECIDE);
        /// but we would get warnings due to the lack allocation

        // possible todo: create a PetscTools::SetupMatNoAllocation()

        // In the future, when parallelising, remember to think about MAT_IGNORE_OFF_PROC_ENTRIES (see #1682)

#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
        MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,mNumDofs,mNumDofs,&mSystemLhsMatrix);
        MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,mNumDofs,mNumDofs,&mPreconditionMatrix);
#else //New API
        MatCreate(PETSC_COMM_WORLD,&mSystemLhsMatrix);
        MatCreate(PETSC_COMM_WORLD,&mPreconditionMatrix);
        MatSetSizes(mSystemLhsMatrix,PETSC_DECIDE,PETSC_DECIDE,mNumDofs,mNumDofs);
        MatSetSizes(mPreconditionMatrix,PETSC_DECIDE,PETSC_DECIDE,mNumDofs,mNumDofs);
#endif

        if (PetscTools::IsSequential())
        {
            MatSetType(mSystemLhsMatrix, MATSEQAIJ);
            MatSetType(mPreconditionMatrix, MATSEQAIJ);
            MatSeqAIJSetPreallocation(mSystemLhsMatrix,    PETSC_NULL, num_non_zeros_each_row);
            MatSeqAIJSetPreallocation(mPreconditionMatrix, PETSC_NULL, num_non_zeros_each_row);
        }
        else
        {
            PetscInt lo, hi;
            VecGetOwnershipRange(mResidualVector, &lo, &hi);

            int* num_non_zeros_each_row_this_proc = new int[hi-lo];
            int* zero = new int[hi-lo];
            for (unsigned i=0; i<unsigned(hi-lo); i++)
            {
                num_non_zeros_each_row_this_proc[i] = num_non_zeros_each_row[lo+i];
                zero[i] = 0;
            }

            MatSetType(mSystemLhsMatrix, MATMPIAIJ);
            MatSetType(mPreconditionMatrix, MATMPIAIJ);
            MatMPIAIJSetPreallocation(mSystemLhsMatrix,    PETSC_NULL, num_non_zeros_each_row_this_proc, PETSC_NULL, num_non_zeros_each_row_this_proc);
            MatMPIAIJSetPreallocation(mPreconditionMatrix, PETSC_NULL, num_non_zeros_each_row_this_proc, PETSC_NULL, num_non_zeros_each_row_this_proc);
        }

        MatSetFromOptions(mSystemLhsMatrix);
        MatSetFromOptions(mPreconditionMatrix);

        //unsigned total_non_zeros = 0;
        //for (unsigned i=0; i<mNumDofs; i++)
        //{
        //   total_non_zeros += num_non_zeros_each_row[i];
        //}
        //std::cout << total_non_zeros << " versus " << 500*mNumDofs << "\n" << std::flush;

        delete [] num_non_zeros_each_row;
    }
}
#endif // ABSTRACTCONTINUUMMECHANICSSOLVER_HPP_
