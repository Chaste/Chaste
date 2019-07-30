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


#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

#include "UblasMatrixInclude.hpp"
#include "UblasVectorInclude.hpp"

#include <vector>

#include "AbstractTetrahedralElement.hpp"
#include "Node.hpp"
#include "ChastePoint.hpp"

/**
 * A concrete element class which inherits from AbstractTetrahedralElement.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class Element : public AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index  the index of the element in the mesh
     * @param rNodes  the nodes owned by the element
     * @param registerWithNodes  whether to tell the nodes that they are contained in this element
     */
    Element(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes, bool registerWithNodes=true);

    /**
     * "Copy" constructor which allows a new index to be specified.
     *
     * @param rElement  an element to copy
     * @param index the index of the new element
     */
    Element(const Element& rElement, const unsigned index);

    /**
     * Inform all nodes forming this element that they are in this element.
     */
    void RegisterWithNodes();

    /**
     * Mark the element as having been removed from the mesh.
     * Also notify nodes in the element that it has been removed.
     */
    void MarkAsDeleted();

    /**
     * Update node at the given index.
     *
     * @param rIndex is an local index to which node to change
     * @param pNode is a pointer to the replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode);

    /**
     * Reset the index of this element in the mesh.
     *
     * @param index the new index of the element
     */
    void ResetIndex(unsigned index);

    /**
     * Calculate the circumsphere/circumcircle of this element.
     *
     * After reconstructing a cylindrical 2d mesh, the Jacobian data of the periodic elements is not valid anymore.
     * We want to use the jacobians computed before swapping the nodes.

     * @return a vector containing x_centre, y_centre,...,radius^2
     *
     * @param rJacobian  the Jacobian matrix
     * @param rInverseJacobian  the inverse Jacobian matrix
     */
    c_vector<double,SPACE_DIM+1> CalculateCircumsphere(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian,
                                                       c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian);

    /**
     * @return the volume of the circumsphere, or area of the circumcircle, of this element.
     */
    double CalculateCircumsphereVolume();

    /**
     * @return the quality of a triangle/tetrahedron is the ratio between the
     * volume of the shape and the volume of its circumsphere.
     * This is normalised by dividing through by the Platonic ratio.
     */
    double CalculateQuality();

    /**
     * @return calculated maximum and minimum edge lengths
     */
    c_vector <double, 2> CalculateMinMaxEdgeLengths();

    /**
     * @return calculated interpolation weights: the vector
     *  (1-xi(0)-xi(1)-xi(2), xi(0), xi(1), xi(2))
     * (in the 3D case) for a given point. (see CalculateXi() documentation)
     *
     * @param rTestPoint reference to the point
     */
    c_vector<double, SPACE_DIM+1> CalculateInterpolationWeights(const ChastePoint<SPACE_DIM>& rTestPoint);

    /**
     * @return calculated interpolation weights (see CalculateInterpolationWeights() documentation),
     * but if we are not within the element (one or more negative weights), we project onto the
     * element, rather than extrapolating from it.
     *
     * @param rTestPoint reference to the point
     */
    c_vector<double, SPACE_DIM+1> CalculateInterpolationWeightsWithProjection(const ChastePoint<SPACE_DIM>& rTestPoint);

    /**
     * @return calculated xi at a given point. These are the values in the canonical element, using the
     *  the canonical element coordinate system (relative to node 0), corresponding to the test point in this element.
     *  For example, if the test point is node 0, xi=(0,0,0); if node 2, then xi=(0,1,0); if the
     *  test point is halfway between nodes 1 and 2 on the edge between then, then xi=(0.5,0.5,0);
     *  if the test point is the interior, then xi=(a,b,c), where a,b,c>0 and 1-a-b-c > 0.
     *
     * @param rTestPoint reference to the point
     */
    c_vector<double, SPACE_DIM> CalculateXi(const ChastePoint<SPACE_DIM>& rTestPoint);

    /**
     * @return whether a given point lies inside this element.
     *
     * @param rTestPoint reference to the point
     * @param strict whether the point must not be too close to an edge/face (defaults to false)
     */
    bool IncludesPoint(const ChastePoint<SPACE_DIM>& rTestPoint, bool strict=false);
};

#endif //_ELEMENT_HPP_
