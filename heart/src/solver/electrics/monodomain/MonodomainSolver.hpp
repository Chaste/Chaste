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
 *  This solver uses two assemblers, MonodomainAssembler to assemble the LHS matrix, [(chi*C/dt) M  + K],
 *  and also to compute c_surf, and MassMatrixAssembler to assemble the mass matrix M used on the RHS.
 *  Note that MonodomainAssembler itself calls:
 *  MassMatrixAssembler for M (as this class does directly for RHS),
 *  and MonodomainStiffnessMatrixAssembler for K.
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
     */
    MonodomainSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                     MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                     BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions);

    /**
     *  Destructor
     */
    virtual ~MonodomainSolver();
};

#endif /*MONODOMAINSOLVER_HPP_*/
