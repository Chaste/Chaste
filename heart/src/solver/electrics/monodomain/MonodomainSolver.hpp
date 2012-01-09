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

#ifndef MONODOMAINSOLVER_HPP_
#define MONODOMAINSOLVER_HPP_

#include "AbstractDynamicLinearPdeSolver.hpp"
#include "MassMatrixAssembler.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"
#include "MonodomainCorrectionTermAssembler.hpp"
#include "MonodomainTissue.hpp"
#include "MonodomainAssembler.hpp"

/**
 *  A monodomain solver, which uses various assemblers to set up the
 *  monodomain linear system.
 *
 *  The discretised monodomain equation leads to the linear system (see FEM
 *  implementations document)
 *
 *  ( (chi*C/dt) M  + K ) V^{n+1} = (chi*C/dt) M V^{n} + M F^{n} + c_surf
 *
 *  where chi is the surface-area to volume ratio, C the capacitance, dt the timestep
 *  M the mass matrix, K the stiffness matrix, V^{n} the vector of voltages at time n,
 *  F^{n} the vector of (chi*Iionic + Istim) at each node, and c_surf a vector
 *  arising from any surface stimuli (usually zero).
 *
 *  This solver uses two assemblers, one to assemble the LHS matrix, (chi*C/dt) M  + K,
 *  and also to compute c_surf, and one to assemble the mass matrix M.
 *
 *  Also allows state variable interpolation (SVI) to be used on elements for which it
 *  will be needed, if the appropriate HeartConfig boolean is set.
 *  See wiki page ChasteGuides/StateVariableInterpolation for more details on this.
 *  In this case the equation is
 *  ( (chi*C/dt) M  + K ) V^{n+1} = (chi*C/dt) M V^{n} + M F^{n} + c_surf + c_correction
 *  and another assembler is used to create the c_correction.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainSolver
  : public AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>
{
private:

    /** Monodomain tissue class (collection of cells, and conductivities) */
    MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* mpMonodomainTissue;

    /**
     *  Number of quadrature points per dimension (only saved so it can be
     *  passed to the assembler)
     */
    unsigned mNumQuadPoints;

    /** Boundary conditions */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* mpBoundaryConditions;

    /** The monodomain assembler, used to set up the LHS matrix */
    MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>* mpMonodomainAssembler;

    /** Assembler for surface integrals coming from any non-zero Neumann boundary conditions */
    NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,1>* mpNeumannSurfaceTermsAssembler;

    /**
     * If using state variable interpolation, points to an assembler to use in
     * computing the correction term to apply to the RHS.
     */
    MonodomainCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM>* mpMonodomainCorrectionTermAssembler;

    /** The mass matrix, used to computing the RHS vector */
    Mat mMassMatrix;

    /** The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b (excluding surface integrals), this vector is z where b=Mz.
     */
    Vec mVecForConstructingRhs;


    /**
     *  Implementation of SetupLinearSystem() which uses the assembler to compute the
     *  LHS matrix, but sets up the RHS vector using the mass-matrix (constructed
     *  using a separate assembler) multiplied by a vector
     *
     *  @param currentSolution  Solution at current time
     *  @param computeMatrix  Whether to compute the matrix of the linear system
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);

public:
    /**
     *  Overloaded PrepareForSetupLinearSystem() methods which
     *  gets the cell models to solve themselves
     *
     *  @param currentSolution solution at current time
     */
    void PrepareForSetupLinearSystem(Vec currentSolution);

    /**
     *  Overloaded InitialiseForSolve
     *
     *  @param initialSolution initial solution
     */
    virtual void InitialiseForSolve(Vec initialSolution);

    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    MonodomainSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                     MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                     BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                     unsigned numQuadPoints = 2);

    /**
     *  Destructor
     */
    ~MonodomainSolver();
};



#endif /*MONODOMAINSOLVER_HPP_*/
