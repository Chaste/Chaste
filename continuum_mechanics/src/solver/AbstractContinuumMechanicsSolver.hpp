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
 *  General base class for continuum mechanics solvers
 */
template<unsigned DIM>//, bool MIXED_PROBLEM>
class AbstractContinuumMechanicsSolver
{
protected:
    /**
     * The mesh to be solved on. Requires 6 nodes per triangle (or 10 per tetrahedron)
     * as quadratic bases are used.
     */
    QuadraticMesh<DIM>& mrQuadMesh;

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
     * Allocates memory for the Jacobian and preconditioner matrices (larger number of
     * non-zeros per row than with say linear problems)
     */
    void AllocateMatrixMemory();

    /**
     * Residual vector for the full nonlinear system, also the RHS vector in the linear
     * system used to solve the nonlinear problem using Newton's method.
     */
    Vec mResidualVector;

    /**
     * The RHS side in the linear system that is solved each Newton iteration. Since Newton's method
     * is Ju = f, where J is the Jacobian, u the (negative of the) update and f the residual, it might seem necessary
     * to store this as well as the residual. However, when applying Dirichlet boundary conditions in
     * the compressible case, we alter the rows of the matrix, but also alter the columns in order to
     * maintain symmetry. This requires making further changes to the right-hand vector, meaning that
     * it no longer properly represents the residual. Hence, we have to use two vectors.
     *
     * Overall, this can be represents as
     *  - compute residual f
     *  - compute Jacobian J
     *  - apply BCs to f.
     *  - alter the linear system from Ju=f to (J*)u=f* which enforces the dirichlet boundary conditions but enforces them symmetrically.
     *
     * mLinearSystemRhsVector represents f*.
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
     *
     * In the incompressible case:
     * the preconditioner is the petsc LU factorisation of
     *
     * Jp = [A B] in displacement-pressure block form,
     *      [C M]
     *
     * where the A, B and C are the matrices in the normal jacobian,
     * i.e.
     *
     * J  = [A B]
     *      [C 0]
     *
     * and M is the MASS MATRIX (ie integral phi_i phi_j dV, where phi_i are the
     * pressure basis functions).
     */
    Mat mPreconditionMatrix;

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
     * @param counterToAppend append a counter to the file name
     * @param extension to append at end.
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


template<unsigned DIM>//, bool MIXED_PROBLEM>
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

    this->mCurrentSolution.resize(mNumDofs, 0.0);
}


template<unsigned DIM>//, bool MIXED_PROBLEM>
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
        if (mCompressibilityType==COMPRESSIBLE)
        {
            PetscTools::Destroy(mDirichletBoundaryConditionsVector);
        }
    }
}

template<unsigned DIM>//, bool MIXED_PROBLEM>
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

template<unsigned DIM>//, bool MIXED_PROBLEM>
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

    std::vector<double>& r_pressure = this->rGetPressures();
    for (unsigned i=0; i<r_pressure.size(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            *p_file << this->mrQuadMesh.GetNode(i)->rGetLocation()[j] << " ";
        }

        *p_file << r_pressure[i] << "\n";
    }
    p_file->close();
}

template<unsigned DIM>//, bool MIXED_PROBLEM>
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
    mPressureSolution.resize(this->mrQuadMesh.GetNumVertices());

    for (unsigned i=0; i<this->mrQuadMesh.GetNumVertices(); i++)
    {
        mPressureSolution[i] = this->mCurrentSolution[DIM*this->mrQuadMesh.GetNumNodes() + i];
    }
    return mPressureSolution;
}


template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::AllocateMatrixMemory()
{
    ///////////////////////////
    // three vectors
    ///////////////////////////
    mResidualVector = PetscTools::CreateVec(mNumDofs);
    VecDuplicate(mResidualVector, &mLinearSystemRhsVector);
    if (mCompressibilityType == COMPRESSIBLE)
    {
        VecDuplicate(mResidualVector, &mDirichletBoundaryConditionsVector);
    }

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
