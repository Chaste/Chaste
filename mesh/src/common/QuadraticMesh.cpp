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

#include <map>
//#include <pair>

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

template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructLinearMesh(unsigned numElemX)
{
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
void QuadraticMesh<DIM>::ConstructRectangularMesh(unsigned numElemX, unsigned numElemY, bool unused)
{
    assert(DIM==2);

    assert(numElemX > 0);
    assert(numElemY > 0);
    assert(unused);
    ///\todo #2224 The call looks like it is going to apply a stagger, but it does not.
    AbstractTetrahedralMesh<DIM,DIM>::ConstructRectangularMesh(numElemX, numElemY, false);

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
        bool boundary = (j==0) || (j==numElemY);
        //Add mid-way nodes to horizontal edges in this slice
        for (unsigned i=0; i<numElemX; i++)
        {
            unsigned left_index = j*(numElemX+1) + i;
            std::pair<unsigned,unsigned> edge(left_index, left_index+1 ) ;
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
            std::pair<unsigned,unsigned> edge(left_index, left_index+(numElemX+1) ) ;
            edge_to_internal_map[edge] = node_index;
            boundary = (i==0) || (i==numElemX);
            MakeNewInternalNode(node_index, node_pos, top);
//                    unsigned parity=(i+(numElemY-j))%2;
//                    if (parity==1)
            {
                //backslash
                std::pair<unsigned,unsigned> edge(left_index+1, left_index+(numElemX+1) ) ;
                edge_to_internal_map[edge] = node_index;
            }
//                    else
//                    {
//                        //foward slash
//                        std::pair<unsigned,unsigned> edge(left_index, left_index+(numElemX+1)+1 ) ;
//                        edge_to_internal_map[edge] = node_index;
//                    }
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
            if (global_index1 > global_index2)
            {
                unsigned tmp=global_index1;
                global_index1=global_index2;
                global_index2=tmp;
            }
            unsigned node_index = edge_to_internal_map[std::pair<unsigned,unsigned>(global_index1, global_index2)];
            iter->AddNode(this->mNodes[node_index]);
            this->mNodes[node_index]->AddElement(iter->GetIndex());
        }        
    }

    for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter = this->GetBoundaryElementIteratorBegin();
         iter != this->GetBoundaryElementIteratorEnd();
         ++iter)
    {
        unsigned global_index1 =  (*iter)->GetNodeGlobalIndex(0);
        unsigned global_index2 =  (*iter)->GetNodeGlobalIndex(1);
        if (global_index1 > global_index2)
        {
            unsigned tmp=global_index1;
            global_index1=global_index2;
            global_index2=tmp;
        }
        unsigned node_index = edge_to_internal_map[std::pair<unsigned,unsigned>(global_index1, global_index2)];
        (*iter)->AddNode(this->mNodes[node_index]);
        this->mNodes[node_index]->AddBoundaryElement((*iter)->GetIndex());
    }
    
            

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
            return NULL;
        }
        if ( (rLocation[dim] == 0.0) || (rLocation[dim] == rTop[dim]) )
        {
            boundary = true;
        }
    }
    if (DIM == 2)
    {
//        PRINT_3_VARIABLES(rIndex, rLocation[0], rLocation[1]);
    }
    if (DIM == 3)
    {
//        PRINT_4_VARIABLES(rIndex, rLocation[0], rLocation[1], rLocation[2]);
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
void QuadraticMesh<DIM>::ConstructCuboidNewImp(unsigned numElemX, unsigned numElemY, unsigned numElemZ)
{
    assert(DIM==3);

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
        //
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
            if (global_index1 > global_index2)
            {
                unsigned tmp=global_index1;
                global_index1=global_index2;
                global_index2=tmp;
            }
            unsigned node_index = edge_to_internal_map[std::pair<unsigned,unsigned>(global_index1, global_index2)];
            assert(node_index != 0); //Failure here means that we can't find something in the map
            (*iter)->AddNode(this->mNodes[node_index]);
            this->mNodes[node_index]->AddElement((*iter)->GetIndex());
        }
        
    }
}

template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructCuboid(unsigned numElemX, unsigned numElemY, unsigned numElemZ)
{
    assert(DIM==3);

    assert(numElemX > 0);
    assert(numElemY > 0);
    assert(numElemZ > 0);
    this->mMeshIsLinear = false;
    unsigned num_nodes = (numElemX+1)*(numElemY+1)*(numElemZ+1);

    struct tetgen::tetgenio mesher_input;
    mesher_input.pointlist = new double[num_nodes * DIM];
    mesher_input.numberofpoints = num_nodes;
    unsigned new_index = 0;
    for (unsigned k=0; k<=numElemZ; k++)
    {
        double z = k;
        for (unsigned j=0; j<=numElemY; j++)
        {
            double y = j;
            for (unsigned i=0; i<=numElemX; i++)
            {
                double x = i;
                mesher_input.pointlist[DIM*new_index] = x;
                mesher_input.pointlist[DIM*new_index + 1] = y;
                mesher_input.pointlist[DIM*new_index + 2] = z;
                new_index++;
            }
        }
    }

    // Library call
    struct tetgen::tetgenio mesher_output;
    tetgen::tetrahedralize((char*)"Qo2", &mesher_input, &mesher_output, NULL);

    assert(mesher_output.numberofcorners == (DIM+1)*(DIM+2)/2);//Nodes per element (including internals, one per edge)

    this->ImportFromMesher(mesher_output, mesher_output.numberoftetrahedra, mesher_output.tetrahedronlist, mesher_output.numberoftrifaces, mesher_output.trifacelist, NULL);

    CountVertices();
    QuadraticMeshHelper<DIM>::AddNodesToBoundaryElements(this, NULL);
}


template<unsigned DIM>
unsigned QuadraticMesh<DIM>::GetNumVertices() const
{
    return mNumVertices;
}


template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructFromLinearMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader)
{
    assert(DIM != 1);

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
        triangulate((char*)"Qzero2", &mesher_input, &mesher_output, NULL);

        this->ImportFromMesher(mesher_output, mesher_output.numberoftriangles, mesher_output.trianglelist, mesher_output.numberofedges, mesher_output.edgelist, mesher_output.edgemarkerlist);
        CountVertices();
        QuadraticMeshHelper<DIM>::AddNodesToBoundaryElements(this, NULL);

        //Tidy up triangle
        this->FreeTriangulateIo(mesher_input);
        this->FreeTriangulateIo(mesher_output);
    }
    else // in 3D, remesh using tetgen
    {

        struct tetgen::tetgenio mesher_input, mesher_output;

        mesher_input.numberoftetrahedra = this->GetNumElements();
        mesher_input.tetrahedronlist = new int[this->GetNumElements() * (DIM+1)];
        this->ExportToMesher(unused_map, mesher_input, mesher_input.tetrahedronlist);

        // Library call
        tetgen::tetrahedralize((char*)"Qzro2", &mesher_input, &mesher_output);

        this->ImportFromMesher(mesher_output, mesher_output.numberoftetrahedra, mesher_output.tetrahedronlist, mesher_output.numberoftrifaces, mesher_output.trifacelist, NULL);
        CountVertices();
        QuadraticMeshHelper<DIM>::AddNodesToBoundaryElements(this, NULL);
    }
}


template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rAbsMeshReader)
{
    TrianglesMeshReader<DIM, DIM>* p_mesh_reader=dynamic_cast<TrianglesMeshReader<DIM, DIM>*>(&rAbsMeshReader);

    unsigned order_of_elements = 1;
    if (p_mesh_reader)
    {
        //A triangles mesh reader will let you read with non-linear elements
        order_of_elements = p_mesh_reader->GetOrderOfElements();
    }

    // If it is a linear TrianglesMeshReader or any other reader (which are all linear)
    if (order_of_elements == 1)
    {
        WARNING("Reading a (linear) tetrahedral mesh and converting it to a QuadraticMesh.  This involves making an external library call to Triangle/Tetgen in order to compute internal nodes");
        ConstructFromLinearMeshReader(rAbsMeshReader);
        return;
    }

    TetrahedralMesh<DIM,DIM>::ConstructFromMeshReader(*p_mesh_reader);
    assert(this->GetNumBoundaryElements() > 0);

    QuadraticMeshHelper<DIM>::AddInternalNodesToElements(this, p_mesh_reader);
    CountVertices();
    QuadraticMeshHelper<DIM>::AddInternalNodesToBoundaryElements(this, p_mesh_reader);
    QuadraticMeshHelper<DIM>::CheckBoundaryElements(this);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class QuadraticMesh<1>;
template class QuadraticMesh<2>;
template class QuadraticMesh<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(QuadraticMesh)
