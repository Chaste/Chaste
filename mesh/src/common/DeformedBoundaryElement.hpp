/*

Copyright (c) 2005-2012, University of Oxford.
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
