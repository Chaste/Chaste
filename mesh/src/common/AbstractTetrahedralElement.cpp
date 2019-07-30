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

#include "UblasCustomFunctions.hpp"
#include "AbstractTetrahedralElement.hpp"
#include "Exception.hpp"

#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::RefreshJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian)
{
    if (this->mIsDeleted)
    {
        EXCEPTION("Attempting to Refresh a deleted element");
    }
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        for (unsigned j=0; j!=ELEMENT_DIM; j++) // Does a j<ELEMENT_DIM without ever having to test j<0U (#186: pointless comparison of unsigned integer with zero)
        {
            rJacobian(i,j) = this->GetNodeLocation(j+1,i) - this->GetNodeLocation(0,i);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : AbstractElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes)
{
    // Sanity checking
    unsigned num_vectices = ELEMENT_DIM+1;

    double det;

    if (SPACE_DIM == ELEMENT_DIM)
    {
        // This is so we know it's the first time of asking
        // Create Jacobian
        c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacobian;
        try
        {
            CalculateJacobian(jacobian, det);
        }
        catch (Exception)
        {
            // if the Jacobian is negative the orientation of the element is probably
            // wrong, so swap the last two nodes around.

            this->mNodes[num_vectices-1] = rNodes[num_vectices-2];
            this->mNodes[num_vectices-2] = rNodes[num_vectices-1];

            CalculateJacobian(jacobian, det);
            // If determinant < 0 then element nodes are listed clockwise.
            // We want them anticlockwise.
            assert(det > 0.0);
        }
    }
    else
    {
        //This is not a full-dimensional element
        c_vector<double, SPACE_DIM> weighted_direction;
        CalculateWeightedDirection(weighted_direction, det);
    }

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralElement(unsigned index)
    : AbstractElement<ELEMENT_DIM,SPACE_DIM>(index)
{}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double& rJacobianDeterminant)
{

    assert(ELEMENT_DIM <= SPACE_DIM);
    RefreshJacobian(rJacobian);

    {
        rJacobianDeterminant = Determinant(rJacobian);
        if (rJacobianDeterminant <= DBL_EPSILON)
        {
            std::stringstream message;
            message << "Jacobian determinant is non-positive: "
                    << "determinant = " << rJacobianDeterminant
                    << " for element " << this->mIndex << " (" << ELEMENT_DIM
                    << "D element in " << SPACE_DIM << "D space). Nodes are at:" << std::endl;

            for (unsigned local_node_index=0u; local_node_index != ELEMENT_DIM+1; local_node_index++)
            {
                c_vector<double, SPACE_DIM> location = this->GetNodeLocation(local_node_index);
                message << "Node " << this->GetNodeGlobalIndex(local_node_index) << ":\t";

                for (unsigned i=0; i<SPACE_DIM; i++)
                {
                    message << location[i];
                    if (i==SPACE_DIM - 1u)
                    {
                        message << std::endl;
                    }
                    else
                    {
                        message << "\t";
                    }
                }
            }
            EXCEPTION(message.str());
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateWeightedDirection(c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant)
{

    if (ELEMENT_DIM >= SPACE_DIM)
    {
        assert(ELEMENT_DIM == SPACE_DIM);
        EXCEPTION("WeightedDirection undefined for fully dimensional element");
    }

    c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacobian;
    RefreshJacobian(jacobian);

    /*
     * At this point we're only dealing with subspace (ELEMENT_DIM < SPACE_DIM) elem.
     * We assume that the rWeightedDirection vector and rJacobianDeterminant (length
     * of vector) are the values from a previous call.
     */

    // This code is only used when ELEMENT_DIM<SPACE_DIM
    switch (ELEMENT_DIM)
    {
        case 0:
            // See specialised template for ELEMENT_DIM==0
            NEVER_REACHED;
            break;
        case 1:
            // Linear edge in a 2D plane or in 3D
            rWeightedDirection=matrix_column<c_matrix<double,SPACE_DIM,ELEMENT_DIM> >(jacobian, 0);
            break;
        case 2:
            // Surface triangle in a 3d mesh
            assert(SPACE_DIM == 3);
            rWeightedDirection(0) = -SubDeterminant(jacobian, 0, 2);
            rWeightedDirection(1) =  SubDeterminant(jacobian, 1, 2);
            rWeightedDirection(2) = -SubDeterminant(jacobian, 2, 2);
            break;
        default:
           ; // Not going to happen
    }
    rJacobianDeterminant = norm_2(rWeightedDirection);

    if (rJacobianDeterminant < DBL_EPSILON)
    {
        EXCEPTION("Jacobian determinant is zero");
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateNormal()
{
    if (ELEMENT_DIM == 1 && SPACE_DIM == 3)
    {
        EXCEPTION("Don't have enough information to calculate a normal vector");
    }
    c_vector<double, SPACE_DIM> normal;
    double determinant;
    CalculateWeightedDirection(normal, determinant);
    normal /= determinant;
    if (ELEMENT_DIM == 1)
    {
        // Need to rotate so tangent becomes normal
        double x = normal[0];
        normal[0] = normal[1];
        normal[1] = -x;
    }
    return normal;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateCentroid() const
{
    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);
    for (unsigned i=0; i<=ELEMENT_DIM; i++)
    {
        centroid += this->mNodes[i]->rGetLocation();
    }
    return centroid/((double)(ELEMENT_DIM + 1));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateInverseJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double& rJacobianDeterminant, c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian)
{
    assert(ELEMENT_DIM <= SPACE_DIM);     // LCOV_EXCL_LINE
    CalculateJacobian(rJacobian, rJacobianDeterminant);

    // CalculateJacobian should make sure that the determinant is not close to zero (or, in fact, negative)
    assert(rJacobianDeterminant > 0.0);
    rInverseJacobian = Inverse(rJacobian);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::GetVolume(double determinant) const
{

    if (this->mIsDeleted)
    {
        return 0.0;
    }

    double scale_factor = 1.0;

    if (ELEMENT_DIM == 2)
    {
        scale_factor = 2.0; // both the volume of the canonical triangle is 1/2
    }
    else if (ELEMENT_DIM == 3)
    {
        scale_factor = 6.0; // both the volume of the canonical triangle is 1/6
    }
    return determinant/scale_factor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::GetStiffnessMatrixGlobalIndices(unsigned problemDim, unsigned* pIndices) const
{
    for (unsigned local_index=0; local_index<ELEMENT_DIM+1; local_index++)
    {
        unsigned node = this->GetNodeGlobalIndex(local_index);

        for (unsigned problem_index=0; problem_index<problemDim; problem_index++)
        {
            //std::cout << local_index*problemDim + problem_index << std::endl;
            pIndices[local_index*problemDim + problem_index] = node*problemDim + problem_index;
        }
    }
}


//////////////////////////////////////////////////////////////////////
//                  Specialization for 0d elements                  //
//////////////////////////////////////////////////////////////////////

/**
 * Specialization for 0d elements so we don't get errors from Boost on some
 * compilers.
 */
template<unsigned SPACE_DIM>
AbstractTetrahedralElement<0, SPACE_DIM>::AbstractTetrahedralElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : AbstractElement<0, SPACE_DIM>(index, rNodes)
{
    // Sanity checking
    //unsigned total_nodes = 1;
    assert(this->mNodes.size() == 1);
    assert(SPACE_DIM > 0);     // LCOV_EXCL_LINE

    // This is so we know it's the first time of asking
    // Create Jacobian
    c_vector<double, SPACE_DIM> weighted_direction;
    double det;

    CalculateWeightedDirection(weighted_direction, det);

    // If determinant < 0 then element nodes are listed clockwise.
    // We want them anticlockwise.
    assert(det > 0.0);
}

template<unsigned SPACE_DIM>
AbstractTetrahedralElement<0, SPACE_DIM>::AbstractTetrahedralElement(unsigned index)
    : AbstractElement<0, SPACE_DIM>(index)
{
}

template<unsigned SPACE_DIM>
void AbstractTetrahedralElement<0, SPACE_DIM>::CalculateWeightedDirection(
        c_vector<double, SPACE_DIM>& rWeightedDirection,
        double& rJacobianDeterminant)
{
    assert(SPACE_DIM > 0);     // LCOV_EXCL_LINE

    // End point of a line
    rWeightedDirection = zero_vector<double>(SPACE_DIM);
    rWeightedDirection(0) = 1.0;

    rJacobianDeterminant = 1.0;
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractTetrahedralElement<0, SPACE_DIM>::CalculateNormal()
{
    assert(SPACE_DIM > 0);     // LCOV_EXCL_LINE

    // End point of a line
    c_vector<double, SPACE_DIM> normal = zero_vector<double>(SPACE_DIM);
    ///\todo should throw?
    return normal;
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractTetrahedralElement<0, SPACE_DIM>::CalculateCentroid() const
{
    c_vector<double, SPACE_DIM> centroid = this->mNodes[0]->rGetLocation();
    return centroid;
}

template<unsigned SPACE_DIM>
void AbstractTetrahedralElement<0, SPACE_DIM>::GetStiffnessMatrixGlobalIndices(unsigned problemDim, unsigned* pIndices) const
{
    for (unsigned local_index=0; local_index<1; local_index++)
    {
        unsigned node = this->GetNodeGlobalIndex(local_index);

        for (unsigned problem_index=0; problem_index<problemDim; problem_index++)
        {
            //std::cout << local_index*problemDim + problem_index << std::endl;
            pIndices[local_index*problemDim + problem_index] = node*problemDim + problem_index;
        }
    }
}

// Explicit instantiation
template class AbstractTetrahedralElement<0,1>;
template class AbstractTetrahedralElement<1,1>;
template class AbstractTetrahedralElement<0,2>;
template class AbstractTetrahedralElement<1,2>;
template class AbstractTetrahedralElement<2,2>;
template class AbstractTetrahedralElement<0,3>;
template class AbstractTetrahedralElement<1,3>;
template class AbstractTetrahedralElement<2,3>;
template class AbstractTetrahedralElement<3,3>;
