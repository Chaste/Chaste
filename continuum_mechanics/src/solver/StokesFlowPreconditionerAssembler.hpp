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


#ifndef STOKESFLOWPRECONDITIONERASSEMBLER_HPP_
#define STOKESFLOWPRECONDITIONERASSEMBLER_HPP_

#include "StokesFlowAssembler.hpp"


/**
 *  Assembler for setting up (volume-integral parts of) the preconditioner matrix used in
 *  the FEM solver Stokes' Flow.
 *
 *  The system matrix has the block form (except see comment below)
 *  [A   B]
 *  [B^T 0]
 *  In contrast, the preconditioner is:
 *  [A   B]
 *  [B^T M]
 *
 *  The class therefore just needs to inherit from StokesFlowAssembler, which will assemble
 *  the A,B,B^T terms, and it just has to overload the pressure-pressure block method.
 *
 *  NOTE: The elemental matrix and vector is as above. The full matrix and vector uses a completely
 *  different ordering: for parallelisation reasons the pressure variables are interleaved with the
 *  spatial variables and dummy pressure variables are used for internal nodes. For example, in 2d,
 *  the ordering is
 *  [U1 V1 P1 , .. , Un Vn, Pn]
 *  where n is the total number of nodes.
 */
template<unsigned DIM>
class StokesFlowPreconditionerAssembler : public StokesFlowAssembler<DIM>
{
private:
    /** Number of vertices per element  */
    static const unsigned NUM_VERTICES_PER_ELEMENT = DIM+1;

    /**
     * Size of the pressure block, per element, equal num_vertices times 1 (as p solved
     * for at each vertex.
     */
    static const unsigned PRESSURE_BLOCK_SIZE_ELEMENTAL = NUM_VERTICES_PER_ELEMENT;

    /**
     *  For a continuum mechanics problem in mixed form (displacement-pressure or velocity-pressure), the matrix
     *  has the form (except see comments about ordering above)
     *  [A     B1]
     *  [B2^T  C ]
     *  The function is related to the pressure-pressure block, i.e. C
     *
     *  @return C=M, the mass matrix.
     *
     *  @param rLinearPhi  All the linear basis functions on this element, evaluated at the current quad point
     *  @param rGradLinearPhi  Gradients of all the linear basis functions on this element, evaluated at the current quad point
     *  @param rX Current location (physical position)
     *  @param pElement Current element
     */
    c_matrix<double,PRESSURE_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputePressurePressureMatrixTerm(
        c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
        c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        return outer_prod(rLinearPhi,rLinearPhi);
    }

public:
    /**
     * Constructor
     * @param pMesh mesh
     * @param pProblemDefinition problem definition
     */
    StokesFlowPreconditionerAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh,
                                      StokesFlowProblemDefinition<DIM>* pProblemDefinition)
        : StokesFlowAssembler<DIM>(pMesh,pProblemDefinition)
    {
    }
};

#endif // STOKESFLOWPRECONDITIONERASSEMBLER_HPP_
