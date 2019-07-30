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

#include "QuadraticMesh.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "Warnings.hpp"
#include "QuadraticMeshHelper.hpp"

//Jonathan Shewchuk's Triangle and Hang Si's TetGen
#define REAL double
#define VOID void
#include "triangle.h"
#include "tetgen.h"
#undef REAL
#undef VOID

template<unsigned DIM>
void QuadraticMesh<DIM>::CountVertices()
{
    mNumVertices = 0;
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        bool is_internal = this->GetNode(i)->IsInternal();
        if (is_internal==false)
        {
            mNumVertices++;
        }
    }
}

template<unsigned DIM>
QuadraticMesh<DIM>::QuadraticMesh(double spaceStep, double width, double height, double depth)
{
    this->ConstructRegularSlabMesh(spaceStep, width, height, depth);
}

//////////////////////////////////////////////////////////
// Badly-named (name inherited from parent class),
// 'linear' here refers to the fact it creates a 1d mesh
// on a line
//////////////////////////////////////////////////////////
template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructLinearMesh(unsigned numElemX)
{
    assert(DIM==1); // LCOV_EXCL_LINE

    AbstractTetrahedralMesh<DIM,DIM>::ConstructLinearMesh(numElemX);
    assert (this->mNodes.size() == numElemX+1);
    mNumVertices = numElemX+1;
    c_vector<double, DIM> top;
    top[0] = numElemX;

    unsigned mid_node_index=mNumVertices;
    for (unsigned element_index=0; element_index<numElemX; element_index++)
    {
        c_vector<double, DIM> x_value_mid_node;
        x_value_mid_node[0] = element_index+0.5;

        Node<DIM>* p_mid_node = MakeNewInternalNode(mid_node_index, x_value_mid_node, top);

        //Put in element and cross-reference
        this->mElements[element_index]->AddNode(p_mid_node);
        p_mid_node->AddElement(element_index);
    }

    this->RefreshMesh();
}

template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructRectangularMesh(unsigned numElemX, unsigned numElemY, bool stagger)
{
    assert(DIM==2); // LCOV_EXCL_LINE
    assert(numElemX > 0);
    assert(numElemY > 0);

    AbstractTetrahedralMesh<DIM,DIM>::ConstructRectangularMesh(numElemX, numElemY, stagger);

    this->mMeshIsLinear=false;
    //Make the internal nodes in y-order.  This is important for the distributed case, since we want the top and bottom
    //layers to have predictable numbers
    std::map<std::pair<unsigned, unsigned>, unsigned> edge_to_internal_map;

    unsigned node_index = this->GetNumNodes();
    c_vector<double, DIM> top;
    top[0]=numElemX;
    top[1]=numElemY;
    c_vector<double, DIM> node_pos;

    for (unsigned j=0; j<numElemY+1; j++)
    {
        node_pos[1]=j;
        //Add mid-way nodes to horizontal edges in this slice
        for (unsigned i=0; i<numElemX; i++)
        {
            unsigned left_index = j*(numElemX+1) + i;
            std::pair<unsigned,unsigned> edge(left_index, left_index+1 );
            edge_to_internal_map[edge] = node_index;
            node_pos[0]=i+0.5;
            MakeNewInternalNode(node_index, node_pos, top);
        }

        //Add the vertical and diagonal nodes to the mid-way above the last set of horizontal edges
        node_pos[1] = j+0.5;
        for (unsigned i=0; i<numElemX+1; i++)
        {
            node_pos[0] = i;
            unsigned left_index = j*(numElemX+1) + i;
            std::pair<unsigned,unsigned> edge(left_index, left_index+(numElemX+1) );
            edge_to_internal_map[edge] = node_index;
            MakeNewInternalNode(node_index, node_pos, top);
            unsigned parity=(i+(numElemY-j))%2;
            if (stagger==false || parity==1) //Default when no stagger
            {
                //backslash
                std::pair<unsigned,unsigned> back_edge(left_index+1, left_index+(numElemX+1) );
                edge_to_internal_map[back_edge] = node_index;
            }
            else
            {
                //foward slash
                std::pair<unsigned,unsigned> forward_edge(left_index, left_index+(numElemX+1)+1 );
                edge_to_internal_map[forward_edge] = node_index;
            }
            node_pos[0] = i+0.5;
            MakeNewInternalNode(node_index, node_pos, top);
        }
    }
    CountVertices();

//    assert(edge_to_internal_map.size() == this->GetNumNodes()-this->GetNumVertices());
    for (typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        unsigned local_index1=0;
        for (unsigned index=0; index<=DIM; index++)
        {
            local_index1 = (local_index1+1)%(DIM+1);
            unsigned local_index2 = (local_index1+1)%(DIM+1);
            unsigned global_index1 =  iter->GetNodeGlobalIndex(local_index1);
            unsigned global_index2 =  iter->GetNodeGlobalIndex(local_index2);
            unsigned new_node_index = LookupInternalNode(global_index1, global_index2, edge_to_internal_map);
            iter->AddNode(this->mNodes[new_node_index]);
            this->mNodes[new_node_index]->AddElement(iter->GetIndex());
        }
    }

    for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter = this->GetBoundaryElementIteratorBegin();
         iter != this->GetBoundaryElementIteratorEnd();
         ++iter)
    {
        unsigned global_index1 =  (*iter)->GetNodeGlobalIndex(0);
        unsigned global_index2 =  (*iter)->GetNodeGlobalIndex(1);
        unsigned new_node_index = LookupInternalNode(global_index1, global_index2, edge_to_internal_map);
        (*iter)->AddNode(this->mNodes[new_node_index]);
        this->mNodes[new_node_index]->AddBoundaryElement((*iter)->GetIndex());
    }

    this->RefreshMesh();
}

template<unsigned DIM>
Node<DIM>* QuadraticMesh<DIM>::MakeNewInternalNode(unsigned& rIndex, c_vector<double, DIM>& rLocation, c_vector<double, DIM>& rTop)
{
    bool boundary = false;
    for (unsigned dim=0; dim<DIM; dim++)
    {
        if (rLocation[dim] > rTop[dim])
        {
            //Outside the box so don't do anything
            return nullptr;
        }
        if ((rLocation[dim] == 0.0) || (rLocation[dim] == rTop[dim]))
        {
            boundary = true;
        }
    }
    //The caller needs to know that rIndex is in sync with what's in the mesh
    assert(rIndex == this->mNodes.size());
    Node<DIM>* p_node   = new Node<DIM>(rIndex++, rLocation, boundary);
    p_node->MarkAsInternal();
    //Put in mesh
    this->mNodes.push_back(p_node);
    if (boundary)
    {
        this->mBoundaryNodes.push_back(p_node);
    }
    return p_node;
}

template<unsigned DIM>
unsigned QuadraticMesh<DIM>::LookupInternalNode(unsigned globalIndex1, unsigned globalIndex2, std::map<std::pair<unsigned, unsigned>, unsigned>& rEdgeMap)
{
    unsigned node_index = 0u;
    assert(globalIndex1 != globalIndex2);
    if (globalIndex1 < globalIndex2)
    {
        node_index = rEdgeMap[std::pair<unsigned,unsigned>(globalIndex1, globalIndex2)];
    }
    else
    {
        node_index = rEdgeMap[std::pair<unsigned,unsigned>(globalIndex2, globalIndex1)];
    }
    //A failure to find the key would result in a new zero entry in the map.  Note that no *internal* node will have global index zero.
    assert(node_index != 0u);
    return node_index;
}

template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructCuboid(unsigned numElemX, unsigned numElemY, unsigned numElemZ)
{
    assert(DIM==3); // LCOV_EXCL_LINE

    assert(numElemX > 0);
    assert(numElemY > 0);
    assert(numElemZ > 0);

    AbstractTetrahedralMesh<DIM,DIM>::ConstructCuboid(numElemX, numElemY, numElemZ);
    c_vector<double, DIM> top;
    top[0]=numElemX;
    top[1]=numElemY;
    top[2]=numElemZ;
    c_vector<double, DIM> node_pos;
    this->mMeshIsLinear=false;
    //Make the internal nodes in z-order.  This is important for the distributed case, since we want the top and bottom
    //layers to have predictable numbers
    std::map<std::pair<unsigned, unsigned>, unsigned> edge_to_internal_map;
    unsigned node_index = this->GetNumNodes();
    for (unsigned k=0; k<numElemZ+1; k++)
    {
        //Add a slice of the mid-points to the edges and faces at this z=level
        node_pos[2] = k;
        for (unsigned j=0; j<numElemY+1; j++)
        {
            unsigned lo_z_lo_y = (numElemX+1)*((numElemY+1)*k + j);
            unsigned lo_z_hi_y = (numElemX+1)*((numElemY+1)*k + j + 1);

            node_pos[1] = j;

            //The midpoints along the horizontal (y fixed) edges
            for (unsigned i=0; i<numElemX+1; i++)
            {
                // i+.5, j, k
                std::pair<unsigned,unsigned> edge(lo_z_lo_y+i, lo_z_lo_y+i+1);
                edge_to_internal_map[edge] = node_index;
                node_pos[0] = i+0.5;
                MakeNewInternalNode(node_index, node_pos, top);
            }
            //The midpoints and face centres between two horizontal (y-fixed) strips
            node_pos[1] = j+0.5;
            for (unsigned i=0; i<numElemX+1; i++)
            {
                // i, j+0.5, k
                std::pair<unsigned,unsigned> edge(lo_z_lo_y+i, lo_z_hi_y+i);
                edge_to_internal_map[edge] = node_index;
                node_pos[0] = i;
                MakeNewInternalNode(node_index, node_pos, top);
                //Centre of face node
                // i+0.5, j+0.5, k
                std::pair<unsigned,unsigned> edge2(lo_z_lo_y+i, lo_z_hi_y+i+1);
                edge_to_internal_map[edge2] = node_index;
                node_pos[0] = i+0.5;
                MakeNewInternalNode(node_index, node_pos, top);
            }
        }
        //Add a slice of the mid-points to the edges and faces mid-way up the cube z=level
        node_pos[2] = k+0.5;
        for (unsigned j=0; j<numElemY+1; j++)
        {
            node_pos[1] = j;
            unsigned lo_z_lo_y = (numElemX+1)*((numElemY+1)*k + j);
            unsigned hi_z_lo_y = (numElemX+1)*((numElemY+1)*(k+1) + j);
            unsigned hi_z_hi_y = (numElemX+1)*((numElemY+1)*(k+1) + j + 1);

            //The midpoints along the horizontal (y fixed) edges
            for (unsigned i=0; i<numElemX+1; i++)
            {
                // i, j, k+0.5
                std::pair<unsigned,unsigned> edge(lo_z_lo_y+i, hi_z_lo_y+i);
                edge_to_internal_map[edge] = node_index;
                node_pos[0] = i;
                MakeNewInternalNode(node_index, node_pos, top);

                // i+0.5, j, k+0.5
                node_pos[0] = i+0.5;
                std::pair<unsigned,unsigned> edge2(lo_z_lo_y+i, hi_z_lo_y+i+1);
                edge_to_internal_map[edge2] = node_index;
                MakeNewInternalNode(node_index, node_pos, top);
            }
            //The midpoints and face centres between two horizontal (y-fixed) strips
            node_pos[1] = j+0.5;
            for (unsigned i=0; i<numElemX+1; i++)
            {
                // i, j+0.5, k+0.5
                std::pair<unsigned,unsigned> edge(lo_z_lo_y+i, hi_z_hi_y+i);
                edge_to_internal_map[edge] = node_index;
                node_pos[0] = i;
                MakeNewInternalNode(node_index, node_pos, top);
                //Centre of face node on the main diagonal
                // i+0.5, j+0.5, k+0.5
                std::pair<unsigned,unsigned> edge2(lo_z_lo_y+i, hi_z_hi_y+i+1);
                edge_to_internal_map[edge2] = node_index;
                node_pos[0] = i+0.5;
                MakeNewInternalNode(node_index, node_pos, top);
            }
        }
    }
    CountVertices();
    for (typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        /* The standard tetgen ordering of the internal nodes 4,5,6..9 (using the
         * zero-based numbering scheme) is
         * 4 (0,1), 5 (1,2), 6 (0,2) 7 (0,3), 8 (1,3), 9 (2,3)
         * i.e. internal node with local index 4 is half-way between vertex nodes
         * with local indices 0 and 1.
         */
         unsigned v0 = iter->GetNodeGlobalIndex(0);
         unsigned v1 = iter->GetNodeGlobalIndex(1);
         unsigned v2 = iter->GetNodeGlobalIndex(2);
         unsigned v3 = iter->GetNodeGlobalIndex(3);
         unsigned internal_index;

         //4
         internal_index=LookupInternalNode(v0, v1, edge_to_internal_map);
         iter->AddNode(this->mNodes[internal_index]);
         this->mNodes[internal_index]->AddElement(iter->GetIndex());
         //5
         internal_index=LookupInternalNode(v1, v2, edge_to_internal_map);
         iter->AddNode(this->mNodes[internal_index]);
         this->mNodes[internal_index]->AddElement(iter->GetIndex());
         //6
         internal_index=LookupInternalNode(v0, v2, edge_to_internal_map);
         iter->AddNode(this->mNodes[internal_index]);
         this->mNodes[internal_index]->AddElement(iter->GetIndex());
         //7
         internal_index=LookupInternalNode(v0, v3, edge_to_internal_map);
         iter->AddNode(this->mNodes[internal_index]);
         this->mNodes[internal_index]->AddElement(iter->GetIndex());
         //8
         internal_index=LookupInternalNode(v1, v3, edge_to_internal_map);
         iter->AddNode(this->mNodes[internal_index]);
         this->mNodes[internal_index]->AddElement(iter->GetIndex());
         //9
         internal_index=LookupInternalNode(v2, v3, edge_to_internal_map);
         iter->AddNode(this->mNodes[internal_index]);
         this->mNodes[internal_index]->AddElement(iter->GetIndex());

    }
    for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter = this->GetBoundaryElementIteratorBegin();
         iter != this->GetBoundaryElementIteratorEnd();
         ++iter)
    {
        unsigned local_index1=0;
        for (unsigned index=0; index<DIM; index++)
        {
            local_index1 = (local_index1+1)%(DIM);
            unsigned local_index2 = (local_index1+1)%(DIM);
            unsigned global_index1 =  (*iter)->GetNodeGlobalIndex(local_index1);
            unsigned global_index2 =  (*iter)->GetNodeGlobalIndex(local_index2);
            unsigned new_node_index = LookupInternalNode(global_index1, global_index2, edge_to_internal_map);
            (*iter)->AddNode(this->mNodes[new_node_index]);
            this->mNodes[new_node_index]->AddBoundaryElement((*iter)->GetIndex());
        }
    }
    this->RefreshMesh();
}

template<unsigned DIM>
unsigned QuadraticMesh<DIM>::GetNumVertices() const
{
    return mNumVertices;
}

template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructFromLinearMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader)
{
    assert(DIM != 1); // LCOV_EXCL_LINE

    //Make a linear mesh
    TetrahedralMesh<DIM,DIM>::ConstructFromMeshReader(rMeshReader);

    NodeMap unused_map(this->GetNumNodes());

    if (DIM==2)  // In 2D, remesh using triangle via library calls
    {
        struct triangulateio mesher_input, mesher_output;
        this->InitialiseTriangulateIo(mesher_input);
        this->InitialiseTriangulateIo(mesher_output);

        mesher_input.numberoftriangles = this->GetNumElements();
        mesher_input.trianglelist = (int *) malloc(this->GetNumElements() * (DIM+1) * sizeof(int));
        this->ExportToMesher(unused_map, mesher_input, mesher_input.trianglelist);

        // Library call
        triangulate((char*)"Qzero2", &mesher_input, &mesher_output, nullptr);

        this->ImportFromMesher(mesher_output, mesher_output.numberoftriangles, mesher_output.trianglelist, mesher_output.numberofedges, mesher_output.edgelist, mesher_output.edgemarkerlist);
        CountVertices();
        QuadraticMeshHelper<DIM>::AddNodesToBoundaryElements(this, nullptr);

        //Tidy up triangle
        this->FreeTriangulateIo(mesher_input);
        this->FreeTriangulateIo(mesher_output);
    }
    else // in 3D, remesh using tetgen
    {

        class tetgen::tetgenio mesher_input, mesher_output;

        mesher_input.numberoftetrahedra = this->GetNumElements();
        mesher_input.tetrahedronlist = new int[this->GetNumElements() * (DIM+1)];
        this->ExportToMesher(unused_map, mesher_input, mesher_input.tetrahedronlist);

        // Library call
        tetgen::tetrahedralize((char*)"Qzro2", &mesher_input, &mesher_output);

        this->ImportFromMesher(mesher_output, mesher_output.numberoftetrahedra, mesher_output.tetrahedronlist, mesher_output.numberoftrifaces, mesher_output.trifacelist, nullptr);
        CountVertices();
        QuadraticMeshHelper<DIM>::AddNodesToBoundaryElements(this, nullptr);
    }
}

template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rAbsMeshReader)
{
    //Some mesh readers will let you read with non-linear elements
    unsigned order_of_elements = rAbsMeshReader.GetOrderOfElements();

    // If it is a linear mesh reader
    if (order_of_elements == 1)
    {
        WARNING("Reading a (linear) tetrahedral mesh and converting it to a QuadraticMesh.  This involves making an external library call to Triangle/Tetgen in order to compute internal nodes");
        ConstructFromLinearMeshReader(rAbsMeshReader);
        return;
    }

    TetrahedralMesh<DIM,DIM>::ConstructFromMeshReader(rAbsMeshReader);
    assert(this->GetNumBoundaryElements() > 0);

    QuadraticMeshHelper<DIM>::AddInternalNodesToElements(this, &rAbsMeshReader);
    CountVertices();
    QuadraticMeshHelper<DIM>::AddInternalNodesToBoundaryElements(this, &rAbsMeshReader);
    QuadraticMeshHelper<DIM>::CheckBoundaryElements(this);
}

///////// Explicit instantiation///////


template class QuadraticMesh<1>;
template class QuadraticMesh<2>;
template class QuadraticMesh<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(QuadraticMesh)
