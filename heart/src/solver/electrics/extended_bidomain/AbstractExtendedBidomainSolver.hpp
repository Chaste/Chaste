/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef ABSTRACTEXTENDEDBIDOMAINSOLVER_HPP_
#define ABSTRACTEXTENDEDBIDOMAINSOLVER_HPP_

#include "AbstractDynamicLinearPdeSolver.hpp"
#include "ExtendedBidomainTissue.hpp"
#include "HeartConfig.hpp"
#include "ExtendedBidomainAssembler.hpp"

/**
 * Abstract class with shared functionalities for different
 * extended bidomain solvers.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractExtendedBidomainSolver
    : public AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,3>
{

protected:

    /** Whether the simulation involves a perfusing bath */
    bool mBathSimulation;


    /** The tissue object over which to solve the equation */
    ExtendedBidomainTissue<SPACE_DIM>* mpExtendedBidomainTissue;

    /** Boundary conditions */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,3>* mpBoundaryConditions;

    /**
     *  Number of quadrature points per dimension (only saved so it can be
     *  passed to the assembler)
     */
    unsigned mNumQuadPoints;

    /**
     *  The extended bidomain assembler, used to set up the LHS matrix,
     *  (but not used to set uo the RHS)
     */
    ExtendedBidomainAssembler<ELEMENT_DIM,SPACE_DIM>* mpExtendedBidomainAssembler;

    /** Local cache of the configuration singleton instance*/
    HeartConfig* mpConfig;

    /** Used when intialising null-space solver to resolve singularity*/
    bool mNullSpaceCreated;

    /**
     * Used when pinning nodes to resolve singularity.
     * This vector indicates the global indices of the nodes to be pinned
     */
    std::vector<unsigned> mFixedExtracellularPotentialNodes;

    /**
     * Used when removing a single row to resolve singularity and
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
     *  components (ie i=2,5,8,..) is zero.
     *
     *  This check is not made if running in parallel, or in debug mode.
     *
     *  The reason why the sum must be zero: the Fredholm alternative states that a singular system Ax=b has
     *  a solution if and only if v.b=0 for all v in ker(A) (ie all v such that Av=b). The nullspace ker(A)
     *  is one dimensional with basis vector v = (0,0,1,0,0,1....0,0,1), so v.b = sum_{i=2,5,8..} b_i.
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
     *  Called at the end of SetupLinearSystem(), before the system is solver
     *  Here, used to avoid problems with phi_e drifting by one of 3 methods:
     *  pinning nodes, using a null space, or using an "average phi_e = 0" row.
     *  @param existingSolution Solution at current time
     */
    virtual void FinaliseLinearSystem(Vec existingSolution);


    /**
     *  GenerateNullBasis
     *
     *  Called by FinaliseAssembleSystem to get the null basis to use for the particular
     *  formulation of the extended idomain equations used.
     */
    virtual Vec GenerateNullBasis() const;

public:

    /**
     * Constructor stores the mesh and pde and sets up boundary conditions.
     *
     * @param bathSimulation Whether the simulation has a perfusing bath
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the PDE
     * @param pBcc pointer to the boundary conditions container
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    AbstractExtendedBidomainSolver(bool bathSimulation,
                         AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                         ExtendedBidomainTissue<SPACE_DIM>* pTissue,
                         BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 3>* pBcc,
                         unsigned numQuadPoints = 2);

    /**
     * Destructor.
     */
    virtual ~AbstractExtendedBidomainSolver();

    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to
     *  zero. This does not necessarily have to be called. If it is not, phi_e
     *  is only defined up to a constant.
     *
     *  @param fixedExtracellularPotentialNodes the nodes to be fixed.
     *
     *  NOTE: currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> fixedExtracellularPotentialNodes);

    /** Used when removing a single row to resolve singularity and
     * replacing it with a constraint on the average phi_e being zero.
     * It is set from the problem class.
     * @param  rowMeanPhiEZero  indicates the row of the matrix to be replaced.
     */
     void SetRowForAverageOfPhiZeroed(unsigned rowMeanPhiEZero);

     /**
      *  Get the boundary conditions being used
      */
     BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,3>* GetBoundaryConditions()
     {
         return mpBoundaryConditions;
     }

};

#endif /*ABSTRACTEXTENDEDBIDOMAINSOLVER_HPP_*/
