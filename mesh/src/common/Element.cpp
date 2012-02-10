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

#include "Element.hpp"

#include <cfloat>
#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>::Element(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes, bool registerWithNodes)
    : AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes)
{
    if (registerWithNodes)
    {
        RegisterWithNodes();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>::Element(const Element& rElement, const unsigned index)
{
    *this = rElement;
    this->mIndex = index;

    RegisterWithNodes();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void Element<ELEMENT_DIM, SPACE_DIM>::RegisterWithNodes()
{
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->AddElement(this->mIndex);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void Element<ELEMENT_DIM, SPACE_DIM>::MarkAsDeleted()
{
    this->mIsDeleted = true;
    // Update nodes in this element so they know they are not contained by us
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveElement(this->mIndex);
    }
}

/** Update node at the given index
 *  @param rIndex is an local index to which node to change
 *  @param pNode is a pointer to the replacement node
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void Element<ELEMENT_DIM, SPACE_DIM>::UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
{
    assert(rIndex < this->mNodes.size());

    // Remove it from the node at this location
    this->mNodes[rIndex]->RemoveElement(this->mIndex);

    // Update the node at this location
    this->mNodes[rIndex] = pNode;

    // Add element to this node
    this->mNodes[rIndex]->AddElement(this->mIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void Element<ELEMENT_DIM, SPACE_DIM>::ResetIndex(unsigned index)
{
    //std::cout << "ResetIndex - removing nodes.\n" << std::flush;
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
       //std::cout << "Node " << this->mNodes[i]->GetIndex() << " element "<< this->mIndex << std::flush;
       this->mNodes[i]->RemoveElement(this->mIndex);
    }
    //std::cout << "\nResetIndex - done.\n" << std::flush;
    this->mIndex=index;
    RegisterWithNodes();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,SPACE_DIM+1> Element<ELEMENT_DIM, SPACE_DIM>::CalculateCircumsphere(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian)
{
    /*Assuming that x0,y0.. is at the origin then we need to solve
     *
     * ( 2x1 2y1 2z1  ) (x)    (x1^2+y1^2+z1^2)
     * ( 2x2 2y2 2z2  ) (y)    (x2^2+y2^2+z2^2)
     * ( 2x3 2y3 2z3  ) (z)    (x3^2+y3^2+z3^2)
     * where (x,y,z) is the circumcentre
     *
     */

    assert(ELEMENT_DIM == SPACE_DIM);
    c_vector<double, ELEMENT_DIM> rhs;

    for (unsigned j=0; j<ELEMENT_DIM; j++)
    {
        double squared_location=0.0;
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            //mJacobian(i,j) is the i-th component of j-th vertex (relative to vertex 0)
            squared_location += rJacobian(i,j)*rJacobian(i,j);
        }
        rhs[j]=squared_location/2.0;
    }

    c_vector<double, SPACE_DIM> centre;
    centre = prod(rhs, rInverseJacobian);
    c_vector<double, SPACE_DIM+1> circum;
    double squared_radius = 0.0;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        circum[i] = centre[i] + this->GetNodeLocation(0,i);
        squared_radius += centre[i]*centre[i];
    }
    circum[SPACE_DIM] = squared_radius;

    return circum;
}

/**
 * The quality of a triangle/tetrahedron is the ratio between the
 * volume of the shape and the volume of its circumsphere.
 * This is normalised by dividing through by the Platonic ratio.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double Element<ELEMENT_DIM, SPACE_DIM>::CalculateQuality()
{
    assert(SPACE_DIM == ELEMENT_DIM);
    if (SPACE_DIM == 1)
    {
        return 1.0;
    }

    c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacobian;
    c_matrix<double, ELEMENT_DIM, SPACE_DIM> jacobian_inverse;
    double jacobian_determinant;

    this->CalculateInverseJacobian(jacobian, jacobian_determinant, jacobian_inverse);

    c_vector<double, SPACE_DIM+1> circum=CalculateCircumsphere(jacobian, jacobian_inverse);
    if (SPACE_DIM == 2)
    {
        /* Want Q=(Area_Tri / Area_Cir) / (Area_Equilateral_Tri / Area_Equilateral_Cir)
         * Area_Tri = |Jacobian| /2
         * Area_Cir = Pi * r^2
         * Area_Eq_Tri = (3*sqrt(3)/4)*R^2
         * Area_Eq_Tri = Pi * R^2
         * Q= (2*|Jacobian|)/3*sqrt(3)*r^2)
         */
        return 2.0*jacobian_determinant/(3.0*sqrt(3.0)*circum[SPACE_DIM]);
    }
    assert(SPACE_DIM == 3);
    /* Want Q=(Vol_Tet / Vol_CirS) / (Vol_Plat_Tet / Vol_Plat_CirS)
      *  Vol_Tet  = |Jacobian| /6
      *  Vol_CirS = 4*Pi*r^3/3
      *  Vol_Plat_Tet  = 8*sqrt(3)*R^3/27
      *  Vol_Plat_CirS = 4*Pi*R^3/3
     * Q= 3*sqrt(3)*|Jacobian|/ (16*r^3)
      */

    return (3.0*sqrt(3.0)*jacobian_determinant)
           /(16.0*circum[SPACE_DIM]*sqrt(circum[SPACE_DIM]));
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM+1> Element<ELEMENT_DIM, SPACE_DIM>::CalculateInterpolationWeights(const ChastePoint<SPACE_DIM>& rTestPoint)
{
    //Can only test if it's a tetrahedal mesh in 3d, triangles in 2d...
    assert(ELEMENT_DIM == SPACE_DIM);

    c_vector<double, SPACE_DIM+1> weights;

    c_vector<double, SPACE_DIM> xi=CalculateXi(rTestPoint);

    //Copy 3 weights and compute the fourth weight
    weights[0]=1.0;
    for (unsigned i=1; i<=SPACE_DIM; i++)
    {
        weights[0] -= xi[i-1];
        weights[i] = xi[i-1];
    }
    return weights;
}

/**
 * Calculate the interpolation weights, but if we are not within
 * the element (one or more negative weights), we project onto the
 * element, rather than extrapolating from it.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM+1> Element<ELEMENT_DIM, SPACE_DIM>::CalculateInterpolationWeightsWithProjection(const ChastePoint<SPACE_DIM>& rTestPoint)
{
    //Can only test if it's a tetrahedal mesh in 3d, triangles in 2d...
    assert(ELEMENT_DIM == SPACE_DIM);

    c_vector<double, SPACE_DIM+1> weights = CalculateInterpolationWeights(rTestPoint);

    // Check for negative weights and set them to zero.
    bool negative_weight = false;

    for (unsigned i=0; i<=SPACE_DIM; i++)
    {
        if (weights[i] < 0.0)
        {
            weights[i] = 0.0;

            negative_weight = true;
        }
    }

    if (negative_weight == false)
    {
        // If there are no negative weights, there is nothing to do.
        return weights;
    }

    // Renormalise so that all weights add to 1.0.

    // Note that all elements of weights are now non-negative and so the l1-norm (sum of magnitudes) is equivalent to the sum of the elements of the vector
    double sum = norm_1 (weights);

    //l1-norm ought to be above 1 (because we scrubbed negative weights)
    //However, if we scrubbed weights that were the size of the machine precision then we might be close to one (even less than 1).
    assert( sum + DBL_EPSILON >= 1.0);

    //We might skip this division when sum ~= 1
    weights = weights/sum;

    return weights;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> Element<ELEMENT_DIM, SPACE_DIM>::CalculateXi(const ChastePoint<SPACE_DIM>& rTestPoint)
{
    //Can only test if it's a tetrahedal mesh in 3d, triangles in 2d...
    assert(ELEMENT_DIM == SPACE_DIM);

    // Find the location with respect to node 0
///\todo: #1361 ComputeContainingElements and related methods, and methods called by that down to
/// here, should really take in const c_vector& rather than ChastePoints.
    ChastePoint<SPACE_DIM> copy = rTestPoint; // as rGetLocation in line below is not a const method
    c_vector<double, SPACE_DIM> test_location=copy.rGetLocation()-this->GetNodeLocation(0);

    //Multiply by inverse Jacobian
    c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacobian;
    c_matrix<double, ELEMENT_DIM, SPACE_DIM> inverse_jacobian;
    double jacobian_determinant;

    ///\todo #1326 This method shouldn't need a new Jacobian inverse for every Xi
    this->CalculateInverseJacobian(jacobian, jacobian_determinant, inverse_jacobian);

    return prod(inverse_jacobian, test_location);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool Element<ELEMENT_DIM, SPACE_DIM>::IncludesPoint(const ChastePoint<SPACE_DIM>& rTestPoint, bool strict)
{
    //Can only test if it's a tetrahedal mesh in 3d, triangles in 2d...
    assert(ELEMENT_DIM == SPACE_DIM);

    c_vector<double, SPACE_DIM+1> weights=CalculateInterpolationWeights(rTestPoint);

    //If the point is in the simplex then all the weights should be positive

    for (unsigned i=0; i<=SPACE_DIM; i++)
    {
        if (strict)
        {
            //Points can't be close to a face
            if (weights[i] <= 2*DBL_EPSILON)
            {
                return false;
            }
        }
        else
        {
            //Allow point to be close to a face
            if (weights[i] < -2*DBL_EPSILON)
            {
                return false;
            }
        }
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class Element<1,1>;
template class Element<1,2>;
template class Element<1,3>;
template class Element<2,2>;
template class Element<2,3>;
template class Element<3,3>;
