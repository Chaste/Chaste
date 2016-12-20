/*

Copyright (c) 2005-2016, University of Oxford.
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

// Forward declaration prevents circular include chain
template<unsigned DIM>
class Node;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexElement;

#include <set>
#include <vector>


/// ===============================================================
/// Some function that can be added into trunk and relevant for all

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ElementHasNode(const VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, const unsigned nodeIndex);

std::set<unsigned> GetSharedElementIndices(const Node<3>* pNodeA, const Node<3>* pNodeB);

std::set<unsigned> GetSharedFaceIndices(const Node<3>* pNodeA, const Node<3>* pNodeB);

/**
 * Face is only considered a boundary face when all nodes are boundary nodes.
 * @return
 */
bool IsFaceOnBoundary(VertexElement<2, 3>* pFace);

/// ===============================================================
/// Functions for monolayer

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

void SetNodeAsApical(Node<3>* pNode);

void SetNodeAsBasal(Node<3>* pNode);

Monolayer::v_type GetNodeType(const Node<3>* pNode);

bool IsApicalNode(const Node<3>* pNode);

bool IsBasalNode(const Node<3>* pNode);


// VertexElement
void SetFaceAsApical(VertexElement<2, 3>* pFace);

void SetFaceAsBasal(VertexElement<2, 3>* pFace);

void SetFaceAsLateral(VertexElement<2, 3>* pFace);

Monolayer::v_type GetFaceType(const VertexElement<2, 3>* pFace);

bool IsApicalFace(const VertexElement<2, 3>* pFace);

bool IsBasalFace(const VertexElement<2, 3>* pFace);

bool IsLateralFace(const VertexElement<2, 3>* pFace);


void SetElementAsMonolayer(VertexElement<3, 3>* pElement);

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool IsMonolayerElement(const VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement);


VertexElement<2, 3>* GetApicalFace(const VertexElement<3, 3>* pElement);

VertexElement<2, 3>* GetBasalFace(const VertexElement<3, 3>* pElement);

/**
 * Get half the number of nodes of a monolayer vertex element
 * to eliminate /2 everywhere in the code.
 * @param pElement
 */
unsigned MonolayerGetNumNodes(const VertexElement<3, 3>* pElement);

/**
 * Get the lateral face that is contained by this element and contains the two
 * given nodes. pElement is required as node only contains the index, but not the
 * pointer of the elements.
 *
 * @param pElement
 * @param nodeIndexA
 * @param nodeIndexB
 * @return
 */
std::vector<unsigned> GetLateralFace(const VertexElement<3, 3>* pElement, const unsigned nodeIndexA, const unsigned nodeIndexB);

bool GetFaceOrientation(const VertexElement<3, 3>* pElement, const unsigned faceIndex);


//void AddPairNode(VertexElement<3, 3>* pElement, const unsigned index, Node<3>* pBasalNode, Node<3>* pApicalNode);


#endif /* MONOLAYERVERTEXMESHCUSTOMFUNCTIONS_HPP_ */
