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
#ifndef COMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_
#define COMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_

/*
 * NOTE ON COMPILATION ERRORS:
 *
 * (The following applies to IncompressibleNonlinearElasticityAssembler; possibly/probably holds for this class too).
 *
 * This file won't compile with Intel icpc version 9.1.039, with error message:
 * "Terminate with:
  (0): internal error: backend signals"
 *
 * Try recompiling with icpc version 10.0.025.
 */

#include "AbstractNonlinearElasticitySolver.hpp"
#include "AbstractCompressibleMaterialLaw.hpp"
#include "QuadraticMesh.hpp"

class TestStressRecoveror;

/**
 * Finite elasticity solver. Solves static *compressible* nonlinear elasticity
 * problems with arbitrary (compressible) material laws and a body force.
 *
 * Uses quadratic basis functions for displacement, and is therefore outside the
 * other assembler or solver hierarchy.
 */
template<size_t DIM>
class CompressibleNonlinearElasticitySolver : public AbstractNonlinearElasticitySolver<DIM>
{
    friend class TestCompressibleNonlinearElasticitySolver;
    friend class TestStressRecoveror;

protected:

    /** Number of nodes per element. */
    static const size_t NUM_NODES_PER_ELEMENT    = AbstractNonlinearElasticitySolver<DIM>::NUM_NODES_PER_ELEMENT;

    /** Number of vertices per element. */
    static const size_t NUM_VERTICES_PER_ELEMENT = AbstractNonlinearElasticitySolver<DIM>::NUM_VERTICES_PER_ELEMENT;

    /** Number of nodes per boundary element. */
    static const size_t NUM_NODES_PER_BOUNDARY_ELEMENT = AbstractNonlinearElasticitySolver<DIM>::NUM_NODES_PER_BOUNDARY_ELEMENT;

    /**
     * Stencil size - number of unknowns per element (DIM*NUM_NODES_PER_ELEMENT displacement unknowns,
     * no pressure unknowns.
     */
    static const size_t STENCIL_SIZE = DIM*NUM_NODES_PER_ELEMENT;

    /** Boundary stencil size. */
    static const size_t BOUNDARY_STENCIL_SIZE = DIM*NUM_NODES_PER_BOUNDARY_ELEMENT;


    /**
     * Assemble residual or Jacobian on an element, using the current solution
     * stored in mCurrrentSolution. The ordering assumed is (in 2d)
     * rBElem = [u0 v0 u1 v1 .. u5 v5].
     *
     * @param rElement The element to assemble on.
     * @param rAElem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     * @param rAElemPrecond The element's contribution to the matrix passed to PetSC
     *     in creating a preconditioner
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     */
    virtual void AssembleOnElement(Element<DIM, DIM>& rElement,
                                   c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
                                   c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElemPrecond,
                                   c_vector<double, STENCIL_SIZE>& rBElem,
                                   bool assembleResidual,
                                   bool assembleJacobian);

    /**
     * Assemble the residual vector (using the current solution stored
     * in mCurrentSolution, output going to mpLinearSystem->rGetRhsVector),
     * or Jacobian matrix (using the current solution stored in
     * mCurrentSolution, output going to mpLinearSystem->rGetLhsMatrix).
     *
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     */
    void AssembleSystem(bool assembleResidual, bool assembleJacobian);

public:

    /**
     * Constructor
     *
     * @param rQuadMesh The quadratic mesh to solve on
     * @param rProblemDefinition an object defining in particular the body force and boundary conditions
     * @param outputDirectory The output directory
     */
    CompressibleNonlinearElasticitySolver(AbstractTetrahedralMesh<DIM,DIM>& rQuadMesh,
                                          SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                          std::string outputDirectory);

    /** Destructor. */
    virtual ~CompressibleNonlinearElasticitySolver();
};

#endif /*COMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_*/
