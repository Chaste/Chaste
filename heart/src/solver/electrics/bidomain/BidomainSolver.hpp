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


#ifndef BIDOMAINSOLVER_HPP_
#define BIDOMAINSOLVER_HPP_


#include "UblasIncludes.hpp"

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "AbstractBidomainSolver.hpp"
#include "HeartConfig.hpp"
#include "BidomainAssembler.hpp"
#include "BidomainMassMatrixAssembler.hpp"
#include "BidomainCorrectionTermAssembler.hpp"
#include "BidomainNeumannSurfaceTermAssembler.hpp"

/**
 *  A bidomain solver, which uses various assemblers to set up the bidomain
 *  FEM linear system.
 *
 *  The discretised bidomain equation leads to the linear system (see FEM
 *  implementations document)
 *
 *  [ (chi*C/dt) M + K1    K1   ] [ V^{n+1}   ]  =  [  (chi*C/dt) M V^{n} + M F^{n} + c1_surf ]
 *  [        K1            K2   ] [ PhiE^{n+1}]     [              c2_surf                    ]
 *
 *  where chi is the surface-area to volume ratio, C the capacitance, dt the timestep
 *  M the mass matrix, K1 and K2 stiffness matrices, V^{n} and PhiE^{n} the vector of
 *  voltages and phi_e at time n, F^{n} the vector of (chi*Iionic + Istim) at each node,
 *  and c1_surf and c2_surf vectors arising from any surface stimuli (usually zero).
 *
 *  This solver uses two assemblers, one to assemble the whole LHS matrix,
 *  and also to compute c1_surf and c2_surf, and one to assemble the mass matrix M.
 *
 *  Also allows state variable interpolation (SVI) to be used on elements for which it
 *  will be needed, if the appropriate HeartConfig boolean is set.
 *  See wiki page ChasteGuides/StateVariableInterpolation for more details on this. In this
 *  case the vector [c_correction, 0] is added to the above, and another assembler is
 *  used to create the c_correction.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainSolver : public AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>
{
private:
    /** Mass matrix, used to computing the RHS vector (actually: mass-matrix in
     *  voltage-voltage block, zero elsewhere)
     */
    Mat mMassMatrix;

    /**
     *  The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b, this vector is z where b=Mz.
     */
    Vec mVecForConstructingRhs;

    /** The bidomain assembler, used to set up the LHS matrix */
    BidomainAssembler<ELEMENT_DIM,SPACE_DIM>* mpBidomainAssembler;

    /** Assembler for surface integrals coming from any non-zero Neumann boundary conditions */
    BidomainNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM>* mpBidomainNeumannSurfaceTermAssembler;

    /**
     * If using state variable interpolation, points to an assembler to use in
     * computing the correction term to apply to the RHS.
     */
    BidomainCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM>* mpBidomainCorrectionTermAssembler;




    /**
     *  Implementation of SetupLinearSystem() which uses the assembler to compute the
     *  LHS matrix, but sets up the RHS vector using the mass-matrix (constructed
     *  using a separate assembler) multiplied by a vector
     *
     *  @param currentSolution Solution at current time
     *  @param computeMatrix Whether to compute the matrix of the linear system
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);

public:
    /**
     * Constructor
     *
     * @param bathSimulation Whether the simulation involves a perfusing bath
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param pBoundaryConditions pointer to the boundary conditions
     */
    BidomainSolver(bool bathSimulation,
                   AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                   BidomainTissue<SPACE_DIM>* pTissue,
                   BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions);

    ~BidomainSolver();

    /** Overloaded InitialiseForSolve() which calls base version but also
     *  initialises mMassMatrix and mVecForConstructingRhs
     *
     *  @param initialSolution initial solution
     */
    void InitialiseForSolve(Vec initialSolution);
};


#endif /*BIDOMAINSOLVER_HPP_*/

