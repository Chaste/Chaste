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

#include "AirwayBranch.hpp"
#include "Exception.hpp"
#include "UblasCustomFunctions.hpp"

AirwayBranch::AirwayBranch(bool radiusOnEdge) : mpChildOne(nullptr),
                                                mpChildTwo(nullptr),
                                                mpParent(nullptr),
                                                mpSibling(nullptr),
                                                mIndex(UINT_MAX),
                                                mRadiusOnEdge(radiusOnEdge)
{}

void AirwayBranch::AddElement(Element<1,3>* pElement)
{
    mElements.push_back(pElement);
}

std::list<Element<1,3>* > AirwayBranch::GetElements()
{
    return mElements;
}

double AirwayBranch::GetLength()
{
    double length = 0.0;

    c_matrix<double, 3, 1> jacobian; //not used
    double element_length = 0.0;

    for (std::list<Element<1,3>*>::iterator iter = mElements.begin();
         iter != mElements.end();
         ++iter)
    {
        (*iter)->CalculateJacobian(jacobian, element_length);
        length += element_length;
    }

    return length;
}

double AirwayBranch::GetAverageRadius()
{
    double length = 0.0;
    double radius = 0.0;

    c_matrix<double, 3, 1> jacobian; //not used
    double element_length = 0.0;

    for (std::list<Element<1,3>* >::iterator iter = mElements.begin();
         iter != mElements.end();
         ++iter)
    {
        (*iter)->CalculateJacobian(jacobian, element_length);
        length += element_length;

        if (mRadiusOnEdge)
        {
            radius += element_length*((*iter)->GetAttribute());
        }
        else
        {
            radius += element_length*((*iter)->GetNode(0)->rGetNodeAttributes()[0] + (*iter)->GetNode(1)->rGetNodeAttributes()[0])/2.0;
        }
    }

    return radius/length;
}

double AirwayBranch::GetPoiseuilleResistance()
{
    double resistance = 0.0;

    c_matrix<double, 3, 1> jacobian; //not used
    double element_length = 0.0;

    for (std::list<Element<1,3>*>::iterator iter = mElements.begin();
         iter != mElements.end();
         ++iter)
    {
        (*iter)->CalculateJacobian(jacobian, element_length);

        double radius = 0.0;
        if (mRadiusOnEdge)
        {
            radius = (*iter)->GetAttribute();
        }
        else
        {
            radius = ((*iter)->GetNode(0)->rGetNodeAttributes()[0] + (*iter)->GetNode(1)->rGetNodeAttributes()[0])/2.0;
        }

        resistance += element_length/SmallPow(radius, 4);
    }

    return resistance;
}

c_vector<double, 3> AirwayBranch::GetDirection()
{
    c_vector<double, 3> direction = this->GetDistalNode()->rGetLocation() - this->GetProximalNode()->rGetLocation();
    direction = direction/norm_2(direction);

    return direction;
}

bool AirwayBranch::IsMajor()
{
    if (this->GetSibling() == nullptr)
    {
        return true;
    }

    return (GetAverageRadius() > GetSibling()->GetAverageRadius());
}

double AirwayBranch::GetBranchAngle()
{
    if (this->GetParent() == nullptr)
    {
        EXCEPTION("Insufficient airway tree structure to calculate branch angle.");
    }

    c_vector<double, 3> dir = GetDirection();
    c_vector<double, 3> parent_dir = this->GetParent()->GetDirection();

    return std::acos(inner_prod(dir,parent_dir)/(norm_2(dir)*norm_2(parent_dir)));
}

double AirwayBranch::GetRotationAngle()
{
    if (this->GetParent() == nullptr || this->GetParent()->GetSibling() == nullptr || this->GetSibling() == nullptr)
    {
        EXCEPTION("Insufficient airway tree structure to calculate rotation angle.");
    }

    c_vector<double, 3> n1 = VectorProduct(GetDirection(), GetSibling()->GetDirection());
    c_vector<double, 3> n2 = VectorProduct(GetParent()->GetDirection(), GetParent()->GetSibling()->GetDirection());

    double rotation_factor = inner_prod(n1,n2)/(norm_2(n1)*norm_2(n2));

    // Sometimes the bifurcations are co-planar, which leads to an undefined angle (0.0 or pi radians are valid).
    // For our purposes we consider this angle to be zero.
    if (fabs(rotation_factor) == 1.0)
    {
        return 0.0;
    }
    else
    {
        return std::acos(rotation_factor);
    }
}

std::vector<AirwayBranch*> AirwayBranch::GetAllChildren()
{
    return mAllChildren;
}

void AirwayBranch::AddChild(AirwayBranch* pChild)
{
    mAllChildren.push_back(pChild);
}

AirwayBranch* AirwayBranch::GetChildOne()
{
    return mpChildOne;
}

AirwayBranch* AirwayBranch::GetChildTwo()
{
    return mpChildTwo;
}

void AirwayBranch::SetChildOne(AirwayBranch* pChildOne)
{
    mpChildOne = pChildOne;
}

void AirwayBranch::SetChildTwo(AirwayBranch* pChildTwo)
{
    mpChildTwo = pChildTwo;
}

AirwayBranch* AirwayBranch::GetSibling()
{
    return mpSibling;
}

AirwayBranch* AirwayBranch::GetParent()
{
    return mpParent;
}

void AirwayBranch::SetSibling(AirwayBranch* pSibling)
{
    mpSibling = pSibling;
}

void AirwayBranch::SetParent(AirwayBranch* pParent)
{
    mpParent = pParent;
}

unsigned AirwayBranch::GetIndex()
{
    return mIndex;
}

void AirwayBranch::SetIndex(unsigned index)
{
    mIndex = index;
}

Node<3>* AirwayBranch::GetProximalNode()
{
    if (mElements.size() < 2) // if we only have one element
    {
        // Assume for now nodes are ordered correctly.  ///\todo: make this more robust
        return mElements.front()->GetNode(0);
    }
    else // we have more than one element
    {
        /**
         * We have this situation:
         *
         *  x  <-- only in first element
         *  |
         *  x  <-- in both first and second element
         *  |
         *  x  <-- only in second element
         *
         *  Find the global index of first element node 1.  Check against both nodes in second element.
         *  If either matches, first element node 0 must be proximal.  Else first element node 1 must be.
         */

        // Get pointers to both elements
        Element<1,3>* p_first_element = mElements.front();
        Element<1,3>* p_second_element = *(++(mElements.begin()));

        assert(p_first_element->GetNumNodes()==2);
        assert(p_second_element->GetNumNodes()==2);

        // Get necessary global indices
        unsigned first_elem_node_1_global_idx = p_first_element->GetNode(1)->GetIndex();
        UNUSED_OPT(first_elem_node_1_global_idx);
        unsigned second_elem_node_0_global_idx = p_second_element->GetNode(0)->GetIndex();
        UNUSED_OPT(second_elem_node_0_global_idx);

        // In case of failure, there's a problem with node ordering in the mesh: look at commented code below.
        assert(first_elem_node_1_global_idx == second_elem_node_0_global_idx);
        return p_first_element->GetNode(0);
//        unsigned second_elem_node_1_global_idx = p_second_element->GetNode(1)->GetIndex();
//
//        // Do logic check as described above
//        if (first_elem_node_1_global_idx == second_elem_node_0_global_idx)
//        {
//            return p_first_element->GetNode(0);
//        }
//        else if (first_elem_node_1_global_idx == second_elem_node_1_global_idx)
//        {
//            return p_first_element->GetNode(0);
//        }
//        else
//        {
//            return p_first_element->GetNode(1);
//        }
    }
}

Node<3>* AirwayBranch::GetDistalNode()
{
    if (mElements.size() < 2)    // if we only have one element
    {
        // Assume for now nodes are ordered correctly.  \todo: make this more robust
        return mElements.front()->GetNode(1);
    }
    else // we have more than one element
    {

        /**
         * We have this situation:
         *
         *  x  <-- only in penultimate element
         *  |
         *  x  <-- in both penultimate and last element
         *  |
         *  x  <-- only in last element
         *
         *  Find the global index of last element node 0.  Check against both nodes in penultimate element.
         *  If either matches, last element node 1 must be distal.  Else last element node 0 must be.
         */

        // Get pointers to both elements
        Element<1,3>* p_last_element = mElements.back();
        Element<1,3>* p_penultimate_element = *(--(--(mElements.end())));

        assert(p_last_element->GetNumNodes()==2);
        assert(p_penultimate_element->GetNumNodes()==2);

        // Get necessary global indices
        unsigned last_elem_node_0_global_idx = p_last_element->GetNode(0)->GetIndex();
        UNUSED_OPT(last_elem_node_0_global_idx);
        unsigned penultimate_elem_node_1_global_idx = p_penultimate_element->GetNode(1)->GetIndex();
        UNUSED_OPT(penultimate_elem_node_1_global_idx);

        // In case of failure, there's a problem with node ordering in the mesh: look at commented code below.
        assert(last_elem_node_0_global_idx == penultimate_elem_node_1_global_idx);
        return p_last_element->GetNode(1);

//        unsigned penultimate_elem_node_0_global_idx = p_penultimate_element->GetNode(0)->GetIndex();
//
//        // Do logic check as described above
//        if (last_elem_node_0_global_idx == penultimate_elem_node_1_global_idx)
//        {
//            return p_last_element->GetNode(1);
//        }
//        else if (last_elem_node_0_global_idx == penultimate_elem_node_0_global_idx)
//        {
//            return p_last_element->GetNode(1);
//        }
//        else
//        {
//            return p_last_element->GetNode(0);
//        }
    }
}

double AirwayBranch::GetBranchVolume()
{
    assert(!mRadiusOnEdge);

    double volume = 0.0;

    for (std::list<Element<1,3>* >::iterator iter = mElements.begin();
         iter != mElements.end();
         ++iter)
    {
        Element<1,3>* current_elem = *iter;

        double r1 = current_elem->GetNode(0)->rGetNodeAttributes()[0];
        double r2 = current_elem->GetNode(1)->rGetNodeAttributes()[0];

        double elem_length = norm_2(current_elem->GetNodeLocation(0) - current_elem->GetNodeLocation(1));

        // Volume of a truncated cone
        volume += (M_PI * (r1*r1 + r1*r2 + r2*r2) * elem_length / 3);
    }

    return volume;
}

double AirwayBranch::GetBranchLateralSurfaceArea()
{
    assert(!mRadiusOnEdge);

    double lateralSurfaceArea = 0.0;

    for (std::list<Element<1,3>* >::iterator iter = mElements.begin();
         iter != mElements.end();
         ++iter)
    {
        Element<1,3>* current_elem = *iter;

        double r1 = current_elem->GetNode(0)->rGetNodeAttributes()[0];
        double r2 = current_elem->GetNode(1)->rGetNodeAttributes()[0];

        double elem_length = norm_2(current_elem->GetNodeLocation(0) - current_elem->GetNodeLocation(1));

        // Volume of a truncated cone
        lateralSurfaceArea += (M_PI * (r1 + r2) * sqrt((r1-r2)*(r1-r2) + elem_length * elem_length));
    }

    return lateralSurfaceArea;
}

c_vector<double, 3> AirwayBranch::GetBranchCentroid()
{
    double volume = 0.0;

    c_vector<double, 3> centroid;
    centroid.clear();

    for (std::list<Element<1,3>* >::iterator iter = mElements.begin();
         iter != mElements.end();
         ++iter)
    {
        Element<1,3>* current_elem = *iter;

        double r1 = current_elem->GetNode(0)->rGetNodeAttributes()[0];
        double r2 = current_elem->GetNode(1)->rGetNodeAttributes()[0];

        c_vector<double, 3> along_element = current_elem->GetNodeLocation(1) - current_elem->GetNodeLocation(0);
        double elem_length = norm_2(along_element);

        // Volume of a truncated cone
        double local_volume = (M_PI * (r1*r1 + r1*r2 + r2*r2) * elem_length / 3);
        volume += local_volume;

        // Centroid height of truncated cone
        double centroid_offset = elem_length * (r1*r1 + 2*r1*r2 + 3*r2*r2) / (4 * (r1*r1 + r1*r2 + r2*r2));

        // Increment centroid by average node position, weighted by local volume
        centroid += local_volume * ( current_elem->GetNodeLocation(0) + centroid_offset * (along_element / elem_length) );
    }

    return centroid / volume;
}

bool AirwayBranch::IsTerminal()
{
    Node<3>* p_distal_node = GetDistalNode();

    return p_distal_node->IsBoundaryNode();
}
