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

#ifndef QUADRATICMESHHELPER_HPP_
#define QUADRATICMESHHELPER_HPP_

#include "AbstractTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

/**
 * A helper class for QuadraticMesh and DistributedQuadratic mesh, containing utility methods
 * that are common to both classes (and only require public access).
 */
template<unsigned DIM>
class QuadraticMeshHelper
{
public:
    /**
     * Top level method for adding internal nodes to the mesh's elements.
     *
     * Note that the internal nodes must already exist in the mesh; this method just associates them with
     * the relevant elements.
     *
     * @param pMesh  the mesh to modify
     * @param pMeshReader  pointer to the reader for accessing the on-disk mesh data
     */
    static void AddInternalNodesToElements(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                           AbstractMeshReader<DIM,DIM>* pMeshReader);

    /**
     * Top level method for adding internal nodes to the mesh's boundary elements.
     *
     * Note that the internal nodes must already exist in the mesh; this method just associates them with
     * the relevant boundary elements.  If the mesh on disk has internal node indices in the boundary
     * elements file, then we just look them up.  Otherwise we figure out which normal element has each
     * boundary element as a face, and extract the internal nodes from that.
     *
     * @param pMesh  the mesh to modify
     * @param pMeshReader  pointer to the reader for accessing the on-disk mesh data
     */
    static void AddInternalNodesToBoundaryElements(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                                   AbstractMeshReader<DIM,DIM>* pMeshReader);


    /**
     * Used by AddInternalNodesToBoundaryElements in the case where the on-disk mesh doesn't have
     * the associations already, and by the linear -> quadratic mesh conversion routines.
     *
     * If there is an on-mesh disk that has containing element information in its edge/face file,
     * then this can be supplied to speed up finding the full element containing each boundary element.
     *
     * @param pMesh  the mesh to modify
     * @param pMeshReader  pointer to the reader for accessing the on-disk mesh data, if any; NULL otherwise
     */
    static void AddNodesToBoundaryElements(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                            AbstractMeshReader<DIM,DIM>* pMeshReader);

    /**
     * Check whether all the boundary elements in the given mesh have the expected number of nodes.
     *
     * @param pMesh  the mesh to check
     */
    static void CheckBoundaryElements(AbstractTetrahedralMesh<DIM, DIM>* pMesh);


    /**
     * This method adds the given node (defined by an element and a node index)
     * to the given boundary element, and also sets the node as a boundary
     * node and adds it to the mesh's boundary nodes.
     *
     * @param pMesh  the (quadratic) mesh to modify
     * @param pBoundaryElement  pointer to a boundary element in the mesh
     * @param pElement  pointer to an element in the mesh
     * @param internalNode  index of a node in the mesh
     */
    static void AddNodeToBoundaryElement(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                         BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                         Element<DIM,DIM>* pElement,
                                         unsigned internalNode);

    /**
     * This method adds the given node to the given boundary element, and also sets
     * the node as a boundary node and adds it to the mesh's boundary nodes.
     *
     * @param pMesh  the (quadratic) mesh to modify
     * @param pBoundaryElement  pointer to a boundary element in the mesh
     * @param pNode  pointer to a node in the mesh
     */
    static void AddNodeToBoundaryElement(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                         BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                         Node<DIM>* pNode);

    /**
     * Given a face in an element (defined by giving an element and the opposite
     * node number to the face) that corresponds to a given boundary element,
     * this method adds in the face's internal nodes to the boundary element
     * (in the correct order).
     *
     * @param pMesh  the (quadratic) mesh to modify
     * @param pBoundaryElement  pointer to a boundary element in the mesh
     * @param pElement  pointer to an element in the mesh
     * @param nodeIndexOppositeToFace  index of a node in the mesh
     */
    static void AddExtraBoundaryNodes(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                      BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                      Element<DIM,DIM>* pElement,
                                      unsigned nodeIndexOppositeToFace);

    /**
     * Nasty helper method for AddNodeToBoundaryElement() in 3D.
     *
     * This method takes in the three vertices of a face which match the given boundary
     * element, and figure out if the order of the nodes in the face is reversed in
     * the boundary element (returned in the bool 'rReverse'). Also, the offset between
     * the first node in the face (as given to this method) and the first node in
     * the boundary element is computed (returned in the variable 'rOffset'). Offset
     * should then be applied before reverse to match the face nodes to the boundary
     * element nodes.
     *
     * \todo document these parameters
     *
     * @param boundaryElemNode0
     * @param boundaryElemNode1
     * @param pElement
     * @param node0
     * @param node1
     * @param node2
     * @param rOffset
     * @param rReverse
     */
    static void HelperMethod1(unsigned boundaryElemNode0, unsigned boundaryElemNode1,
                              Element<DIM,DIM>* pElement,
                              unsigned node0, unsigned node1, unsigned node2,
                              unsigned& rOffset,
                              bool& rReverse);

    /**
     * Nasty helper method for AddNodeToBoundaryElement() in 3D.
     *
     * This method takes the three internal nodes for some face in some element,
     * applies the given offset and reverse (see HelperMethod1) to them, to get
     * the ordered internal nodes which should given to the boundary element.
     * It then calls AddNodeToBoundaryElement with each of the three internal nodes.
     *
     * \todo document these parameters
     *
     * @param pMesh
     * @param pBoundaryElement
     * @param pElement
     * @param internalNode0
     * @param internalNode1
     * @param internalNode2
     * @param offset
     * @param reverse
     */
    static void HelperMethod2(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                              BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                              Element<DIM,DIM>* pElement,
                              unsigned internalNode0, unsigned internalNode1, unsigned internalNode2,
                              unsigned offset,
                              bool reverse);
};

#endif // QUADRATICMESHHELPER_HPP_
