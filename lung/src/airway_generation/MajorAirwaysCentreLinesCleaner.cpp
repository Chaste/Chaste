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

#include "MajorAirwaysCentreLinesCleaner.hpp"


MajorAirwaysCentreLinesCleaner::MajorAirwaysCentreLinesCleaner(MutableMesh<1,3>& rMesh,
                                                               unsigned rootIndex) : mrMesh(rMesh),
                                                                                     mOutletNodeIndex(rootIndex),
                                                                                     mWalker(rMesh, rootIndex),
                                                                                     mCalculator(rMesh, rootIndex)
{
}

void MajorAirwaysCentreLinesCleaner::CleanUsingHorsfieldOrder(unsigned order)
{
    mMaxOrder = order;

    Node<3>* p_node = mrMesh.GetNode(mOutletNodeIndex);
    Element<1,3>* p_element = mrMesh.GetElement(*(p_node->ContainingElementsBegin()));
    CleanElementUsingHorsfieldOrder(p_element, false);
}

void  MajorAirwaysCentreLinesCleaner::CleanElementUsingHorsfieldOrder(Element<1,3>* pElement, bool delete_me)
{
    std::vector<Element<1,3>* > child_eles = mWalker.GetChildElements(pElement);

    // If either of the children of this element are below the order limit delete BOTH of them.
    bool delete_children = delete_me; //Always delete children if this element is to be deleted.
    for (unsigned i = 0; i < child_eles.size(); ++i)
    {
        if (mWalker.GetElementHorsfieldOrder(child_eles[i]) <= mMaxOrder)
        {
            delete_children = true;
        }
    }

    // Recursively process children
    for (unsigned i = 0; i < child_eles.size(); ++i)
    {
        CleanElementUsingHorsfieldOrder(child_eles[i], delete_children);
    }

    // Remove this element if necessary
    if (delete_me)
    {
        mrMesh.DeleteElement(pElement->GetIndex());
    }
}

void MajorAirwaysCentreLinesCleaner::CleanTerminalsHueristic()
{
    std::vector<AirwayBranch*> p_branches = mCalculator.GetBranches();

    for (std::vector<AirwayBranch*>::iterator branch_iter = p_branches.begin();
         branch_iter != p_branches.end();
         ++branch_iter)
    {
         std::list<Element<1,3>* > eles = (*branch_iter)->GetElements();

         if (mWalker.GetElementHorsfieldOrder(eles.front()) == 1u) //Only clean terminals
         {
             double parent_length = (*branch_iter)->GetParent()->GetLength();
             double length = (*branch_iter)->GetLength();
             double mean_radius = (*branch_iter)->GetAverageRadius();

             // Smooth out radii
             for (std::list<Element<1,3>* >::iterator ele_iter = eles.begin();
                  ele_iter != eles.end();
                  ++ele_iter)
             {
                 (*ele_iter)->GetNode(1)->rGetNodeAttributes()[0] = mean_radius;
             }

             // If it's too long then shorten it
             if ((length / parent_length) > 0.9)
             {
                 Node<3>* p_start_node = eles.front()->GetNode(0);
                 const c_vector<double, 3>& start_location = p_start_node->rGetLocation();
                 c_vector<double, 3> branch_direction = (*branch_iter)->GetDirection();

                 for (std::list<Element<1,3>* >::iterator ele_iter = eles.begin();
                     ele_iter != eles.end();
                     ++ele_iter)
                 {
                     mrMesh.DeleteElement((*ele_iter)->GetIndex());
                 }

                 c_vector<double,3> new_terminal_point = start_location + 0.8*parent_length*branch_direction;
                 Node<3>* p_new_node = new Node<3>(0u, new_terminal_point);
                 mrMesh.AddNode(p_new_node);

                 p_new_node->AddNodeAttribute(mean_radius);

                 std::vector<Node<3>*> nodes;
                 nodes.push_back(p_start_node);
                 nodes.push_back(p_new_node);

                 Element<1,3>* p_element = new Element<1,3>(UINT_MAX, nodes);
                 mrMesh.AddElement(p_element);
             }
             else if ((length / parent_length) < 0.6) // If it's too short then lengthen it
             {
                 Node<3>* p_start_node = eles.front()->GetNode(0);
                 const c_vector<double, 3>& start_location = p_start_node->rGetLocation();
                 c_vector<double, 3> branch_direction = (*branch_iter)->GetDirection();

                 for (std::list<Element<1,3>* >::iterator ele_iter = eles.begin();
                      ele_iter != eles.end();
                      ++ele_iter)
                 {
                     mrMesh.DeleteElement((*ele_iter)->GetIndex());
                 }

                 c_vector<double,3> new_terminal_point = start_location + 0.7*parent_length*branch_direction;
                 Node<3>* p_new_node = new Node<3>(0u, new_terminal_point);
                 mrMesh.AddNode(p_new_node);

                 p_new_node->AddNodeAttribute(mean_radius);

                 std::vector<Node<3>*> nodes;
                 nodes.push_back(p_start_node);
                 nodes.push_back(p_new_node);

                 Element<1,3>* p_element = new Element<1,3>(UINT_MAX, nodes);
                 mrMesh.AddElement(p_element);
             }
         }
    }
}

void MajorAirwaysCentreLinesCleaner::CleanIsolatedNodes()
{
    for (AbstractTetrahedralMesh<1,3>::NodeIterator iter = mrMesh.GetNodeIteratorBegin();
         iter != mrMesh.GetNodeIteratorEnd();
         ++iter)
    {
        if (iter->ContainingElementsBegin() == iter->ContainingElementsEnd()) //There are no containing elements
        {
            mrMesh.DeleteNodePriorToReMesh(iter->GetIndex());
        }
    }
}
