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

#include <map>
#include "Cylindrical2dNodesOnlyMesh.hpp"

Cylindrical2dNodesOnlyMesh::Cylindrical2dNodesOnlyMesh(double width)
        : NodesOnlyMesh<2>(),
          mWidth(width)
{
    assert(width > 0.0);
}

void Cylindrical2dNodesOnlyMesh::SetUpBoxCollection(double cutOffLength, c_vector<double, 2*2> domainSize, int numLocalRows, bool isPeriodic)
{
    if (cutOffLength < mWidth)
    {
        EXCEPTION("Need to specify a cut off length larger than the width with Cylindrical2dNodesOnlyMeshes.");
    }

    NodesOnlyMesh<2>::SetUpBoxCollection(cutOffLength, domainSize, PETSC_DECIDE, true);    // Only difference is that this "true" makes the boxes periodic.

    this->AddNodesToBoxes();
}

double Cylindrical2dNodesOnlyMesh::GetWidth(const unsigned& rDimension) const
{
    double width = 0.0;
    assert(rDimension==0 || rDimension==1);
    if (rDimension==0)
    {
        width = mWidth;
    }
    else
    {
        width = MutableMesh<2,2>::GetWidth(rDimension);
    }
    return width;
}


void Cylindrical2dNodesOnlyMesh::SetNode(unsigned nodeIndex, ChastePoint<2> point, bool concreteMove)
{
    // concreteMove should always be false for NodesOnlyMesh as no elements to check
    assert(!concreteMove);

    double x_coord = point.rGetLocation()[0];

    // Perform a periodic movement if necessary
    if (x_coord >= mWidth)
    {
        // Move point to the left
        point.SetCoordinate(0, x_coord - mWidth);
    }
    else if (x_coord < 0.0)
    {
        // Move point to the right
        point.SetCoordinate(0, x_coord + mWidth);
    }

    // Update the node's location
    this->GetNode(nodeIndex)->SetPoint(point);
}

unsigned Cylindrical2dNodesOnlyMesh::AddNode(Node<2>* pNewNode)
{
    // Call method on parent class
    unsigned node_index = NodesOnlyMesh<2>::AddNode(pNewNode);

    // If necessary move it to be back on the cylinder
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point);

    return node_index;

}

c_vector<double, 2> Cylindrical2dNodesOnlyMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    c_vector<double, 2> vector = rLocation2 - rLocation1;
    vector[0] = fmod(vector[0], mWidth);

    /*
     * Handle the cylindrical condition here: if the points are more
     * than halfway around the cylinder apart, measure the other way.
     */
    if (vector[0] > 0.5*mWidth)
    {
        vector[0] -= mWidth;
    }
    else if (vector[0] < -0.5*mWidth)
    {
        vector[0] += mWidth;
    }
    return vector;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Cylindrical2dNodesOnlyMesh)
