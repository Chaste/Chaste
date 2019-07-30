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

#ifndef STOKESFLOWASSEMBLER_HPP_
#define STOKESFLOWASSEMBLER_HPP_

#include "AbstractContinuumMechanicsAssembler.hpp"
#include "StokesFlowProblemDefinition.hpp"


/**
 *  Assembler for setting up (volume-integral parts of) the matrix and vector used in
 *  the FEM discretisation of Stokes' Flow.
 *
 *  The matrix has the block-form (except see comment below)
 *  [A   B]
 *  [B^T 0]
 *  and the vector has the block form (except see comment below)
 *  [b]
 *  [0]
 *
 *  NOTE: The elemental matrix and vector is as above. The full matrix and vector uses a completely
 *  different ordering: for parallelisation reasons the pressure variables are interleaved with the
 *  spatial variables and dummy pressure variables are used for internal nodes. For example, in 2d,
 *  the ordering is
 *  [U1 V1 P1 , .. , Un Vn, Pn]
 *  where n is the total number of nodes.
 */
template<unsigned DIM>
class StokesFlowAssembler : public AbstractContinuumMechanicsAssembler<DIM,true,true>
{
friend class TestStokesFlowAssembler;

private:
    /** Number of vertices per element  */
    static const unsigned NUM_VERTICES_PER_ELEMENT = DIM+1;

    /** Number of nodes per element. */
    static const unsigned NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2; // assuming quadratic

    /**
     * Size of the spatial block, per element, equal num_nodes times DIM (as say in 2d (u,v) solved
     * for at each node
     */
    static const unsigned SPATIAL_BLOCK_SIZE_ELEMENTAL = DIM*NUM_NODES_PER_ELEMENT;

    /**
     * Size of the pressure block, per element, equal num_vertices times 1 (as p solved
     * for at each vertex.
     */
    static const unsigned PRESSURE_BLOCK_SIZE_ELEMENTAL = NUM_VERTICES_PER_ELEMENT;

    /** Stokes' flow problem definition */
    StokesFlowProblemDefinition<DIM>* mpProblemDefinition;

    /** This variable is initialised to 1.0 and almost never changed, and is used in the spatial-spatial matrix term.
     *  One test (see TestStokesFlowAssembler)
     *  sets it to 0.0 before assembling the matrix. Basically, a different weak form USED to be implemented here (corresponding to
     *  one kind of Neumann boundary condition), and the matrix for that weak form was compared against exact solutions in this test.
     *  mScaleFactor = 0.0 corresponds to old weak form, mScaleFactor = 1 corresponds to new weak form as documented in fem implementations
     *  pdf.
     */
    double mScaleFactor;


    /**
     *  The matrix has the form (except see comments about ordering above)
     *  [A    B ]
     *  [B^T  0 ]
     *  The function is related to the spatial-spatial block, ie matrix A.
     *
     *  For the (volume-integral) contribution to A from a given element, this method returns the
     *  INTEGRAND in the definition of A.
     *
     *  @param rQuadPhi  All the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rGradQuadPhi  Gradients of all the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rX Current location (physical position corresponding to quad point)
     *  @param pElement Current element
     *  @return stencil matrix
     */
    c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL> ComputeSpatialSpatialMatrixTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL> ret = zero_matrix<double>(SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL);

        double mu = mpProblemDefinition->GetViscosity();

        for (unsigned index1=0; index1<NUM_NODES_PER_ELEMENT*DIM; index1++)
        {
            unsigned spatial_dim1 = index1%DIM;
            unsigned node_index1 = (index1-spatial_dim1)/DIM;

            for (unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
            {
                unsigned spatial_dim2 = index2%DIM;
                unsigned node_index2 = (index2-spatial_dim2)/DIM;

                ret(index1,index2) +=   mu
                                      * mScaleFactor // virtually always 1, see doxygen for this variable
                                      * rGradQuadPhi(spatial_dim1, node_index2)
                                      * rGradQuadPhi(spatial_dim2, node_index1);

                for (unsigned k=0; k<DIM; k++)
                {
                    ret(index1,index2) +=   mu
                                          * (spatial_dim1==spatial_dim2)
                                          * rGradQuadPhi(k, node_index1)
                                          * rGradQuadPhi(k, node_index2);
                }
            }
        }
        return ret;

    }

    /**
     *  The matrix has the form (except see comments about ordering above)
     *  [A    B ]
     *  [B^T  0 ]
     *  The function is related to the spatial-pressure block, ie matrix B.
     *
     *  For the (volume-integral) contribution to B from a given element, this method returns the
     *  INTEGRAND in the definition of B.
     *
     *  @param rQuadPhi  All the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rGradQuadPhi  Gradients of all the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rLinearPhi  All the linear basis functions on this element, evaluated at the current quad point
     *  @param rGradLinearPhi  Gradients of all the linear basis functions on this element, evaluated at the current quad point
     *  @param rX Current location (physical position corresponding to quad point)
     *  @param pElement Current element
     *  @return stencil matrix
     */
    c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputeSpatialPressureMatrixTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
        c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ret = zero_matrix<double>(SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL);

        for (unsigned index1=0; index1<NUM_NODES_PER_ELEMENT*DIM; index1++)
        {
            unsigned spatial_dim1 = index1%DIM;
            unsigned node_index1 = (index1-spatial_dim1)/DIM;

            for (unsigned index2=0; index2<NUM_VERTICES_PER_ELEMENT; index2++)
            {
                ret(index1,index2) += -rGradQuadPhi(spatial_dim1, node_index1) * rLinearPhi(index2);
            }
        }

        return ret;
    }

    // We don't implement this method - so it is a zero block
    //c_matrix<double,PRESSURE_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputePressurePressureMatrixTerm(
    //    c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
    //    c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
    //    c_vector<double,DIM>& rX,
    //    Element<DIM,DIM>* pElement)


    /**
     *  The matrix has the form (except see comments about ordering above)
     *  [A    B ]
     *  [B^T  0 ]
     *  and the vector has the form
     *  [b1]
     *  [b2]
     *  The function is related to the spatial-block in the vector, ie b1.
     *
     *  For the contribution to b1 from a given element, this method should return the INTEGRAND in the definition of b1.
     *
     *  @param rQuadPhi  All the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rGradQuadPhi  Gradients of all the quadratic basis functions on this element, evaluated at the current quad point
     *  @param rX Current location (physical position)
     *  @param pElement Current element
     *  @return stencil vector
     */
    c_vector<double,SPATIAL_BLOCK_SIZE_ELEMENTAL> ComputeSpatialVectorTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        c_vector<double,SPATIAL_BLOCK_SIZE_ELEMENTAL> ret = zero_vector<double>(SPATIAL_BLOCK_SIZE_ELEMENTAL);

        c_vector<double,DIM> body_force = mpProblemDefinition->GetBodyForce(rX, 0.0);

        for (unsigned index=0; index<NUM_NODES_PER_ELEMENT*DIM; index++)
        {
            unsigned spatial_dim = index%DIM;
            unsigned node_index = (index-spatial_dim)/DIM;

            ret(index) += body_force(spatial_dim) * rQuadPhi(node_index);
        }

        return ret;
    }

    // We don't implement this method - so it is a zero block of the vector:
    //c_vector<double,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputePressureVectorTerm(
    //        c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
    //        c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
    //        c_vector<double,DIM>& rX,
    //        Element<DIM,DIM>* pElement)

public:
    /**
     * Constructor
     * @param pMesh
     * @param pProblemDefinition
     */
    StokesFlowAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh,
                        StokesFlowProblemDefinition<DIM>* pProblemDefinition)
        : AbstractContinuumMechanicsAssembler<DIM,true,true>(pMesh),
          mpProblemDefinition(pProblemDefinition),
          mScaleFactor(1.0)
    {
    }
};

#endif // STOKESFLOWASSEMBLER_HPP_
