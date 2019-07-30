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

#ifndef ABSTRACTBIDOMAINSOLVER_HPP_
#define ABSTRACTBIDOMAINSOLVER_HPP_

#include "AbstractDynamicLinearPdeSolver.hpp"
#include "BidomainTissue.hpp"
#include "HeartConfig.hpp"

/**
 *  Abstract Bidomain class containing some common functionality
 *  Inherits from AbstractDynamicLinearPdeSolver so child classes must
 *  implement SetupLinearSystem()
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractBidomainSolver : public AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,2>
{
protected:
    /** Whether the simulation involves a perfusing bath */
    bool mBathSimulation;

    /** The PDE to be solved. */
    BidomainTissue<SPACE_DIM>* mpBidomainTissue;

    /** Boundary conditions */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* mpBoundaryConditions;

    /** Used when intialising null-space solver to resolve singularity*/
    bool mNullSpaceCreated;

    /** Local cache of the configuration singleton instance*/
    HeartConfig* mpConfig;

    /** Used when pinning nodes to resolve singularity.
     * This vector indicates the global indices of the nodes to be pinned
     */
    std::vector<unsigned> mFixedExtracellularPotentialNodes;

    /** Used when removing a single row to resolve singularity and
     * replacing it with a constraint on the average phi_e being zero.
     * This number indicates the row of the matrix to be replaced.  This is INT_MAX if unset.
     * It is set from the problem class.
     */
    unsigned mRowForAverageOfPhiZeroed;

    /**
     * Create the linear system object if it hasn't been already.
     * Can use an initial solution as PETSc template, or base it on the mesh size.
     *
     * @param initialSolution an initial guess
     */
    void InitialiseForSolve(Vec initialSolution);

    /**
     *  Checks whether the linear system will have a solution (if so, infinite solutions) instead of
     *  zero solutions. The condition is, if the linear system is Ax=b, that sum b_i over for all the PHI_E
     *  components (ie i=1,3,5,..) is zero.
     *
     *  This check is not made if running in parallel, or in debug mode.
     *
     *  The reason why the sum must be zero: the Fredholm alternative states that a singular system Ax=b has
     *  a solution if and only if v.b=0 for all v in ker(A) (ie all v such that Av=b). The nullspace ker(A)
     *  is one dimensional with basis vector v = (0,1,0,1....,0,1), so v.b = sum_{i=1,3,5..} b_i.
     */
    virtual void CheckCompatibilityCondition();

    /**
     *  PrepareForSetupLinearSystem
     *
     *  Called at the beginning of SetupLinearSystem(). Here, used to integrate cell
     *  model odes.
     *  @param existingSolution is the voltage to feed into the cell models
     */
    void PrepareForSetupLinearSystem(Vec existingSolution);

    /**
     *  FinaliseAssembleSystem
     *
     *  Called at the end of SetupLinearSystem(), before the system is solved.
     *
     *  If no dirichlet boundary conditions
     *  (i) Check compatibility condition to check we are solving
     *      a linear system that can be solved
     *  Then either:
     *  (a) If not setting average(phi)=0, we are solving a singular system,
     *      so set up a null space.
     *  (b) Apply average(phi)=0 constraint by altering the last row, to
     *      get a non-singular system
     *
     *  @param existingSolution Solution at current time
     */
    virtual void FinaliseLinearSystem(Vec existingSolution);

    /**
     *  @return null basis vector
     *
     *  Called by FinaliseAssembleSystem to get the null basis to use for the particular
     *  formulation of the bidomain equations used.
     */
    virtual Vec GenerateNullBasis() const;

    /**
     *  Apply any changes needed to the linear system for problems that
     *  include a bath. Checks the voltage-voltage block of the matrix
     *  at bath-nodes is zero, and puts a 1 on the diagonal.
     *
     *  Precondition: This method requires the system matrix to be in assembled
     *  state. Call FinaliseLhsMatrix() on your linear system if required.
     *
     *  @param computeMatrix Whether the LHS matrix of the linear system has
     *   been computed
     *  @param computeVector Whether the RHS vector of the linear system has
     *   been computed
     */
    void FinaliseForBath(bool computeMatrix, bool computeVector);

public:
    /**
     * Constructor
     *
     * @param bathSimulation  whether the simulation has a perfusing bath
     * @param pMesh  pointer to the mesh
     * @param pTissue  pointer to the tissue
     * @param pBoundaryConditions  pointer to the boundary conditions container
     */
    AbstractBidomainSolver(bool bathSimulation,
                           AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                           BidomainTissue<SPACE_DIM>* pTissue,
                           BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions);

    /**
     *  Destructor
     */
    virtual ~AbstractBidomainSolver();

    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to
     *  zero. This does not necessarily have to be called. If it is not, phi_e
     *  is only defined up to a constant.
     *
     *  @param fixedExtracellularPotentialNodes the nodes to be fixed.
     *
     *  @note currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> fixedExtracellularPotentialNodes);

    /**
     * Used when removing a single row to resolve singularity and
     * replacing it with a constraint on the average phi_e being zero.
     * It is set from the problem class.
     * @param rowMeanPhiEZero  indicates the row of the matrix to be replaced.  Ought to be an odd number...
     */
     void SetRowForAverageOfPhiZeroed(unsigned rowMeanPhiEZero);

    /**
     *  @return the boundary conditions being used
     */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* GetBoundaryConditions()
    {
        return mpBoundaryConditions;
    }

    /**
     *  Reset the boundary conditions being used. The caller should deal with deleting
     *  the old bcc pointer.
     *  @param pBcc The new boundary conditions container.
     */
    void ResetBoundaryConditionsContainer(BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBcc)
    {
        assert(pBcc);
        mpBoundaryConditions = pBcc;
        // Note, in SetupLinearSystem()
        //   neumann_assembler->ResetBoundaryConditionsContainer(mpBcc)
        // MUST be called every timestep, in case the bcc was reset.
    }
};

#endif /*ABSTRACTBIDOMAINSOLVER_HPP_*/
