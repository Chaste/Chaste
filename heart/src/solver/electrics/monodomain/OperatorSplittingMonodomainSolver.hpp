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


#ifndef OPERATORSPLITTINGMONODOMAINSOLVER_HPP_
#define OPERATORSPLITTINGMONODOMAINSOLVER_HPP_


#include "MonodomainTissue.hpp"
#include "MonodomainAssembler.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "MassMatrixAssembler.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"

/**
 *  A monodomain solver that uses Strang operator splitting of the diffusion (conductivity) term and the reaction
 *  (ionic current) term, instead of solving the full reaction-diffusion PDE. This does NOT refer to operator splitting
 *  of the two PDEs in the bidomain equations. For details see for example Sundnes et al "Computing the Electrical
 *  Activity of the Heart".
 *
 *  The algorithm is, for solving from t=T to T+dt.
 *
 *  (i)   Solve ODEs   dV/dt = Iionic             for t=T to T+dt/2     [giving updated V (internally, and in solution vector) and updated state variables]
 *  (ii)  Solve PDE    dV/dt = div sigma grad V   for t=T to dt         [using V from step i, --> updated V]
 *  (iii) Solve ODEs   dV/dt = Iionic             for t=T+dt/2 to T+dt  [using V from step ii, --> final V]
 *
 *  Notes
 *   (a)  Stages (iii) and (i) can normally be solved together in one go, except just before/after printing the voltage to file.
 *        However for simplicity of code this has not been implemented
 *   (b)  Therefore, the effective ODE timestep will be:  min(ode_dt, pde_dt/2), where ode_dt and pde_dt are those
 *        given via HeartConfig.
 *   (c)  This solver is FOR COMPARING ACCURACY, NOT PERFORMANCE. It has not been optimised and may or may not
 *        perform well in parallel.
 *   (d)  We don't implement the simpler form of operator splitting, Godunov splitting, where the ODEs are
 *        solved for one timestep and the PDEs are solved for one timestep, since this is formally equivalent
 *        to the default implementation where the ionic current is interpolated from the nodal values
 *        (ie ICI - see ICI/SVI discussion in documentation)
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class OperatorSplittingMonodomainSolver : public AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>
{
private:

    /** Boundary conditions */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* mpBoundaryConditions;

    /** Monodomain tissue class (collection of cells, and conductivities) */
    MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* mpMonodomainTissue;

    /** The monodomain assembler, used to set up the LHS matrix */
    MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>* mpMonodomainAssembler;

    /** Assembler for surface integrals coming from any non-zero Neumann boundary conditions */
    NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,1>* mpNeumannSurfaceTermsAssembler;

    /** The mass matrix, used to computing the RHS vector*/
    Mat mMassMatrix;

    /**
     *  The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b, this vector is z where b=Mz.
     *
     *  In the normal solver this is equal to alpha*V + beta*Iionic + gamma*Istim, here there
     *  is no Iionic term.
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

    /**
     *  Called before setting up the linear system, used to solve the cell models for first half timestep (step (i) above)
     *  @param currentSolution the latest solution vector
     */
    void PrepareForSetupLinearSystem(Vec currentSolution);

    /**
     *  Called after solving the linear system, used to solve the cell models for second half timestep (step (iii) above)
     *  @param currentSolution the latest solution vector (ie the solution of the linear system).
     */
    void FollowingSolveLinearSystem(Vec currentSolution);

public:

    /** Overloaded InitialiseForSolve() which calls base version but also
     *  initialises #mMassMatrix and #mVecForConstructingRhs.
     *
     *  @param initialSolution  initial solution
     */
    void InitialiseForSolve(Vec initialSolution);


    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param pBoundaryConditions pointer to the boundary conditions
     */
    OperatorSplittingMonodomainSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                      MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                                      BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions);

    /**
     *  Destructor
     */
    ~OperatorSplittingMonodomainSolver();
};

#endif /* OPERATORSPLITTINGMONODOMAINSOLVER_HPP_ */
