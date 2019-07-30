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


#include "AirwayTreeWalker.hpp"

AirwayTreeWalker::AirwayTreeWalker(AbstractTetrahedralMesh<1,3>& rAirwaysMesh,
                                    unsigned rootIndex=0u) :
                                    mMesh(rAirwaysMesh),
                                    mOutletNodeIndex(rootIndex),
                                    mNodesAreGraphOrdered(true)
{
    if (mMesh.GetNode(mOutletNodeIndex)->IsBoundaryNode() == false)
    {
        EXCEPTION("Outlet node is not a boundary node");
    }
    Node<3>* p_node = mMesh.GetNode(mOutletNodeIndex);
    assert(p_node->GetNumContainingElements() == 1u);

    //Get the head element & process
    Element<1,3>* p_root_element = mMesh.GetElement(*(p_node->ContainingElementsBegin()));
    mOutletElementIndex = p_root_element->GetIndex();
    ProcessElement(p_root_element, p_node);
    CalculateElementProperties(p_root_element);
}

Element<1,3>* AirwayTreeWalker::GetParentElement(Element<1,3>* pElement)
{
    if (mParentElementMap.count(pElement->GetIndex()) == 0)
    {
        return nullptr;
    }

    return mMesh.GetElement(mParentElementMap[pElement->GetIndex()]);
}

unsigned AirwayTreeWalker::GetParentElementIndex(Element<1,3>* pElement)
{
    return mParentElementMap[pElement->GetIndex()];
}

Element<1,3>* AirwayTreeWalker::GetParentElement(unsigned index)
{
    if (mParentElementMap.count(index) == 0)
    {
        return nullptr;
    }

    return mMesh.GetElement(mParentElementMap[index]);
}

unsigned AirwayTreeWalker::GetParentElementIndex(unsigned index)
{
    return mParentElementMap[index];
}


unsigned AirwayTreeWalker::GetNumberOfChildElements(Element<1,3>* pElement)
{
    return mChildElementsMap[pElement->GetIndex()].size();
}

unsigned AirwayTreeWalker::GetNumberOfChildElements(unsigned index)
{
    return mChildElementsMap[index].size();
}

std::vector<Element<1,3>* > AirwayTreeWalker::GetChildElements(Element<1,3>* pElement)
{
    std::vector<Element<1,3>* > child_elements;

    for (unsigned i = 0; i < mChildElementsMap[pElement->GetIndex()].size(); ++i)
    {
        child_elements.push_back(mMesh.GetElement(mChildElementsMap[pElement->GetIndex()][i]));
    }

    return child_elements;
}

std::vector<unsigned> AirwayTreeWalker::GetChildElementIndices(Element<1,3>* pElement)
{
    return mChildElementsMap[pElement->GetIndex()];
}


Node<3>* AirwayTreeWalker::GetDistalNode(Element<1,3>* pElement)
{
    return mMesh.GetNode(mDistalNodeMap[pElement->GetIndex()]);
}

unsigned AirwayTreeWalker::GetDistalNodeIndex(Element<1,3>* pElement)
{
    return mDistalNodeMap[pElement->GetIndex()];
}

void AirwayTreeWalker::ProcessElement(Element<1,3>* pElement, Node<3>* pParentNode)
{
    for (unsigned node_index = 0; node_index < pElement->GetNumNodes(); ++node_index)
    {
        Node<3>* p_current_node = pElement->GetNode(node_index);

        if (p_current_node != pParentNode) //Prevent moving back up the tree
        {
            mDistalNodeMap[pElement->GetIndex()] = p_current_node->GetIndex();
            if (pParentNode->GetIndex() > p_current_node->GetIndex() )
            {
                mNodesAreGraphOrdered = false;
            }
            for (Node<3>::ContainingElementIterator ele_iter = p_current_node->ContainingElementsBegin();
                ele_iter != p_current_node->ContainingElementsEnd();
                ++ele_iter)
            {
                Element<1,3>* p_child_element = mMesh.GetElement(*ele_iter);

                if (p_child_element != pElement) //Prevent moving back up the tree
                {
                    mParentElementMap[p_child_element->GetIndex()] = pElement->GetIndex();
                    mChildElementsMap[pElement->GetIndex()].push_back(p_child_element->GetIndex());

                    ProcessElement(p_child_element, p_current_node);
                }
            }
        }
    }
}


unsigned AirwayTreeWalker::GetElementGeneration(Element<1,3>* pElement)
{
    return mElementGenerations[pElement->GetIndex()];
}

unsigned AirwayTreeWalker::GetElementGeneration(unsigned element_index)
{
    return mElementGenerations[element_index];
}

//
// Auxiliary predicate to simplify finding the maximum entry in a generation
//
bool less_than_predicate(const std::pair<unsigned, unsigned>& lhs, const std::pair<unsigned, unsigned>& rhs)
{
    return lhs.second < rhs.second;
}

unsigned AirwayTreeWalker::GetMaxElementGeneration()
{
    return std::max_element(mElementGenerations.begin(), mElementGenerations.end(), less_than_predicate)->second;
}

unsigned AirwayTreeWalker::GetElementHorsfieldOrder(Element<1,3>* pElement)
{
    return mElementHorsfieldOrder[pElement->GetIndex()];
}

unsigned AirwayTreeWalker::GetElementHorsfieldOrder(unsigned element_index)
{
    return mElementHorsfieldOrder[element_index];
}

unsigned AirwayTreeWalker::GetMaxElementHorsfieldOrder()
{
    return std::max_element(mElementHorsfieldOrder.begin(), mElementHorsfieldOrder.end(), less_than_predicate)->second;
}


unsigned AirwayTreeWalker::GetElementStrahlerOrder(Element<1,3>* pElement)
{
    return mElementStrahlerOrder[pElement->GetIndex()];
}

unsigned AirwayTreeWalker::GetElementStrahlerOrder(unsigned element_index)
{
    return mElementStrahlerOrder[element_index];
}

unsigned AirwayTreeWalker::GetMaxElementStrahlerOrder()
{
    return std::max_element(mElementStrahlerOrder.begin(), mElementStrahlerOrder.end(), less_than_predicate)->second;
}

void AirwayTreeWalker::CalculateElementProperties(Element<1,3>* pElement)
{
    Element<1,3>* p_parent_element = GetParentElement(pElement);

    //pre-order traversal calculations
    if (p_parent_element == nullptr)
    {
        mElementGenerations[pElement->GetIndex()] = 0u;
    }
    else if (GetNumberOfChildElements(p_parent_element) == 1u)
    {
        mElementGenerations[pElement->GetIndex()] = mElementGenerations[p_parent_element->GetIndex()];
    }
    else //if (GetNumberOfChildElements(p_parent_element) > 1u)
    {
        mElementGenerations[pElement->GetIndex()] = mElementGenerations[p_parent_element->GetIndex()] + 1;
    }

    //Process children
    std::vector<Element<1,3>* > child_eles = GetChildElements(pElement);
    for (unsigned i = 0; i < child_eles.size(); ++i)
    {
        CalculateElementProperties(child_eles[i]);
    }

    //post-order traversal calculations
    if (GetNumberOfChildElements(pElement) == 0u)
    {
        mElementHorsfieldOrder[pElement->GetIndex()] = 1u;
        mElementStrahlerOrder[pElement->GetIndex()] = 1u;
    }
    else if (GetNumberOfChildElements(pElement) == 1u)
    {
        mElementHorsfieldOrder[pElement->GetIndex()] = mElementHorsfieldOrder[child_eles[0]->GetIndex()];
        mElementStrahlerOrder[pElement->GetIndex()] = mElementStrahlerOrder[child_eles[0]->GetIndex()];
    }
    else //if (GetNumberOfChildElements(pElement) > 1u)
    {
        {
            unsigned new_horsfield_order = 0;
            for (unsigned i = 0; i < child_eles.size(); ++i)
            {
                if (mElementHorsfieldOrder[child_eles[i]->GetIndex()] > new_horsfield_order)
                {
                    new_horsfield_order = mElementHorsfieldOrder[child_eles[i]->GetIndex()];
                }
            }
            new_horsfield_order++;
            mElementHorsfieldOrder[pElement->GetIndex()] = new_horsfield_order;
        }

        {
            unsigned new_strahler_order = 0;
            for (unsigned i = 0; i < child_eles.size(); ++i)
            {
                if (mElementStrahlerOrder[child_eles[i]->GetIndex()] > new_strahler_order)
                {
                    new_strahler_order = mElementStrahlerOrder[child_eles[i]->GetIndex()];
                }
                else if (mElementStrahlerOrder[child_eles[i]->GetIndex()] == new_strahler_order)
                {
                    new_strahler_order++;
                }
            }
            mElementStrahlerOrder[pElement->GetIndex()] = new_strahler_order;
        }
    }
}
