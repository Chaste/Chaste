/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef MONOLAYERVERTEXMESHCUSTOMFUNCTIONS_HPP_
#define MONOLAYERVERTEXMESHCUSTOMFUNCTIONS_HPP_

#include <set>
#include <vector>

// Forward declaration prevents circular include chain
template <unsigned DIM>
class Node;
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexElement;
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MutableVertexMesh;

/////////////////////////////////////////////////////////
///      Some function that are relevant for all      ///
/////////////////////////////////////////////////////////

/**
 * A simple function to check if an element or a face has a particular node.
 * @param pElement  the element that might have the node
 * @param nodeIndex  the global index of that particular node
 * @return boolean, if node is in the element
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ElementHasNode(const VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, const unsigned nodeIndex);

/**
 * A simple function to get the element index/indices which own(s) both nodes.
 * @param pNodeA  node number 1
 * @param pNodeB  node number 2
 * @return  index of elements that contain both nodes
 */
std::set<unsigned> GetSharedElementIndices(const Node<3>* pNodeA, const Node<3>* pNodeB);

/**
 * A simple function to get the face index/indices which own(s) both nodes.
 * @param pNodeA  node number 1
 * @param pNodeB  node number 2
 * @return  index of faces that contain both nodes
 */
std::set<unsigned> GetSharedFaceIndices(const Node<3>* pNodeA, const Node<3>* pNodeB);

/**
 * A simple function to get the lateral face which owns both nodes.
 * @param pMesh  pointer of the mesh (this is necessary because nodes only have the indices, could be easier with weak pointer)
 * @param pNodeA  node number 1
 * @param pNodeB  node number 2
 * @return  index of faces that contain both nodes
 */
VertexElement<2, 3>* GetSharedLateralFace(const MutableVertexMesh<3, 3>* pMesh,
                                          const Node<3>* pNodeA, const Node<3>* pNodeB);

/**
 * Output element, its faces and nodes in terminal
 * (adapted from MonolayerVertexMeshGenerator::PrintMesh()).
 * @param pElement pointer of the element which will be inspected.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PrintElement(const VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement);

/**
 * Prints out the elements, nodes and faces in verbose mode.
 * Mainly for debugging purpose.
 *
 * @param printDeletedObjects  whether to include deleted objects (false by default)
 */
void PrintMesh(const MutableVertexMesh<3, 3>* pMesh, const bool printDeletedObjects = false);

/**
 * A face is on boundary when it contains only one element.
 * @param pFace  pointer of the face
 * @return if the face is on boundary
 */
bool IsFaceOnBoundary(const VertexElement<2, 3>* pFace);

///////////////////////////////////////////////////////////////////////////////////
///                       Functions for monolayer classes                       ///
///////////////////////////////////////////////////////////////////////////////////

/**
 * Define some variables within a namespace that will be used by monolayer
 * classes for better readability and prevent variable name conflicts.
 *
 * Note that setter values are different from reference value because nodes
 * and elements store attributes as double, but comparison of double is slightly
 * more troublesome than the comparison of integer types.
 */
namespace Monolayer
{
typedef unsigned v_type;

const double SetBasalValue = 1.1;
const double SetApicalValue = 2.1;
const double SetLateralValue = 3.1;
const double SetElementValue = 4.1;

const v_type BasalValue = 1;
const v_type ApicalValue = 2;
const v_type LateralValue = 3;
const v_type ElementValue = 4;
}

///////////////////////////////////
///     Functions for nodes     ///
///////////////////////////////////

/**
 * Set a node as an apical node.
 * @param pNode  pointer of a new apical node
 */
void SetNodeAsApical(Node<3>* pNode);

/**
 * Set a node as a basal node.
 * @param pNode  pointer of a new basal node
 */
void SetNodeAsBasal(Node<3>* pNode);

/**
 * Getter function for node type. For more details on node type,
 * refer to `namespace Monolayer`.
 * @param pNode  pointer of the interested node
 * @return  the type of node (apical/basal)
 */
Monolayer::v_type GetNodeType(const Node<3>* pNode);

/**
 * @param pNode  pointer of a node
 * @return  whether this node is an apical node
 */
bool IsApicalNode(const Node<3>* pNode);

/**
 * @param pNode  pointer of a node
 * @return  whether this node is a basal node
 */
bool IsBasalNode(const Node<3>* pNode);

//////////////////////////////////
///     Functions for face     ///
//////////////////////////////////

/**
 * Set a face as an apical face.
 * @param pFace  pointer of a new apical face
 */
void SetFaceAsApical(VertexElement<2, 3>* pFace);

/**
 * Set a face as a basal face.
 * @param pFace  pointer of a new apical face
 */
void SetFaceAsBasal(VertexElement<2, 3>* pFace);

/**
 * Set a face as a lateral face.
 * @param pFace  pointer of a new lateral face
 */
void SetFaceAsLateral(VertexElement<2, 3>* pFace);

/**
 * Getter function for face type. For more details on face type,
 * refer to `namespace Monolayer`.
 * @param pFace  pointer of the face of interest
 * @return  the type of face (apical/basal/lateral)
 */
Monolayer::v_type GetFaceType(const VertexElement<2, 3>* pFace);

/**
 * @param pFace  pointer of a face
 * @return  whether this face is an apical face
 */
bool IsApicalFace(const VertexElement<2, 3>* pFace);

/**
 * @param pFace  pointer of a face
 * @return  whether this face is a basal face
 */
bool IsBasalFace(const VertexElement<2, 3>* pFace);

/**
 * @param pFace  pointer of a face
 * @return  whether this face is a lateral face
 */
bool IsLateralFace(const VertexElement<2, 3>* pFace);

////////////////////////////////////////////
///     Functions for Vertex Element     ///
////////////////////////////////////////////

/**
 * Set an element as a monolayer element. This function will at the same time assign
 * its faces their respective types if it hasn't been done. Besides, it will call
 * `VertexElement::MonolayerElementRearrangeFacesNodes` so that other operations can
 * be done more efficiently.
 * @param pElement  pointer of a new monolayer element
 */
void SetElementAsMonolayer(VertexElement<3, 3>* pElement);

/**
 * @param pElement  pointer of an element
 * @return  whether this element is within a monolayer tissue
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool IsMonolayerElement(const VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement);

/**
 * @param pElement  pointer of an element
 * @return  pointer of the apical face of the input element
 */
VertexElement<2, 3>* GetApicalFace(const VertexElement<3, 3>* pElement);

/**
 * @param pElement  pointer of an element
 * @return  pointer of the basal face of the input element
 */
VertexElement<2, 3>* GetBasalFace(const VertexElement<3, 3>* pElement);

/**
 * Get half the number of nodes of a monolayer element (which is much useful than
 * the actual number of nodes in many cases)
 * @param pElement  pointer of an element
 */
unsigned MonolayerGetHalfNumNodes(const VertexElement<3, 3>* pElement);

/**
 * Get the lateral face that is contained by this element and contains the two
 * given nodes. pElement is required by this function as the nodes don't have the
 * pointer of the faces but only its indices.
 * (Function which work with the mesh would be even better)
 * @param pElement  pointer of the element
 * @param nodeIndexA  global index of node A
 * @param nodeIndexB  global index of node B
 * @return  if the nodes doesn't share lateral face, it will return a vector with one
 *          element with UINT_MAX
 *          else, the return vector = {face_global_index, face_orientation, face_local_index}
 */
std::vector<unsigned> GetLateralFace(const VertexElement<3, 3>* pElement, const unsigned nodeIndexA,
                                     const unsigned nodeIndexB);

//void AddPairNode(VertexElement<3, 3>* pElement, const unsigned index, Node<3>* pBasalNode, Node<3>* pApicalNode);

#endif /* MONOLAYERVERTEXMESHCUSTOMFUNCTIONS_HPP_ */
