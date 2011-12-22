/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef DEFORMEDBOUNDARYELEMENT_HPP_
#define DEFORMEDBOUNDARYELEMENT_HPP_

#include "Element.hpp"
#include "BoundaryElement.hpp"



/**
 *  The mechanics solvers in some cases require performing integration over
 *  deformed surfaces, for which this class is a helper class. It takes in
 *  an base boundary element, and a set of displacements for each node of the mesh,
 *  and sets up the corresponding deformed boundary element, on which
 *  CalculateNormal(), ComputeJacobianDeterminant() etc can be called.
 *
 *  Note: just use the vertices of the elements, ignores internal nodes, hence this
 *  class is treats the element as a linear element (ie ignores the possible curved
 *  deformed edges). Hence the outward normal, jacobian determinant etc is constant
 *  across the element.
 */
template<unsigned ELEM_DIM, unsigned SPACE_DIM>
class DeformedBoundaryElement : public BoundaryElement<ELEM_DIM,SPACE_DIM>
{
private:
    /** Number of nodes (or actually, vertices) per boundary element. */
    static const size_t NUM_NODES = ELEM_DIM+1;

public:
    /**
     *  Constructor.
     *
     *  Sets up the nodes of the mesh. Gives them the locations of the canonical
     *  element (although this will be overwritten when the undeformed element is
     *  given in ApplyUndeformedElementAndDisplacement()
     */
    DeformedBoundaryElement();

    /**
     *  Destructor deletes the nodes in the element (as they were created in the constructor)
     */
    ~DeformedBoundaryElement();

    /**
     *  Apply an undeformed boundary element and corresponding displacements of each node
     *  (Note: just displacements of the vertices of the element, not any internal nodes).
     *
     *  Over-writes the nodal locations with the given data.
     *
     *  @param pUndeformedElement The undeformed element
     *  @param rDisplacement Displacement of each node (should be of size ELEM_DIM+1, ie num vertices
     *   in the element).
     */
    void ApplyUndeformedElementAndDisplacement( BoundaryElement<ELEM_DIM,SPACE_DIM>* pUndeformedElement,
                                                std::vector<c_vector<double,SPACE_DIM> >& rDisplacement);

};

#endif /* DEFORMEDBOUNDARYELEMENT_HPP_ */
