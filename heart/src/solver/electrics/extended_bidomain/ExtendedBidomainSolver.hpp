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


#ifndef EXTENDEDBIDOMAINSOLVER_HPP_
#define EXTENDEDBIDOMAINSOLVER_HPP_


#include "UblasIncludes.hpp"

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "AbstractCardiacFeVolumeIntegralAssembler.hpp"
#include "AbstractExtendedBidomainSolver.hpp"
#include "HeartConfig.hpp"
#include "BidomainAssembler.hpp"
#include "BidomainMassMatrixAssembler.hpp"
#include "ExtendedBidomainNeumannSurfaceTermAssembler.hpp"

/**
 *  An extended bidomain solver, which computes the right-hand-side (RHS) vector of the linear
 *  system to be solved using matrix-vector products.
 */
template<unsigned ELEM_DIM, unsigned SPACE_DIM>
class ExtendedBidomainSolver : public AbstractExtendedBidomainSolver<ELEM_DIM,SPACE_DIM>
{
private:
    /**
     * Mass matrix, used to computing the RHS vector (actually: mass-matrix in
     *  voltage-voltage block, zero elsewhere)
     */
    Mat mMassMatrix;

    /**
     *  The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b, this vector is z where b=Mz.
     */
    Vec mVecForConstructingRhs;

    /** The bidomain assembler, used to set up the LHS matrix */
    ExtendedBidomainAssembler<ELEM_DIM,SPACE_DIM>* mpExtendedBidomainAssembler;

    /** Assembler for surface integrals coming from any non-zero Neumann boundary conditions */
    ExtendedBidomainNeumannSurfaceTermAssembler<ELEM_DIM,SPACE_DIM>* mpExtendedBidomainNeumannSurfaceTermAssembler;


    /** Overloaded InitialiseForSolve() which calls base version but also
     *  initialises mMassMatrix and mVecForConstructingRhs
     *
     *  @param initialSolution initial solution
     */
    void InitialiseForSolve(Vec initialSolution);

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
     * @param pTissue pointer to the PDE
     * @param pBoundaryConditions pointer to the boundary conditions
     */
    ExtendedBidomainSolver(bool bathSimulation,
                              AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                              ExtendedBidomainTissue<SPACE_DIM>* pTissue,
                              BoundaryConditionsContainer<ELEM_DIM,SPACE_DIM,3>* pBoundaryConditions);

    /**
     * Destructor
     */
    ~ExtendedBidomainSolver();
};


#endif /*EXTENDEDBIDOMAINSOLVER_HPP_*/

