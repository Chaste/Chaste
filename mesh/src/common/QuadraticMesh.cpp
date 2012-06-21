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

#include "QuadraticMesh.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "Warnings.hpp"
#include "QuadraticMeshHelper.hpp"

//Jonathan Shewchuk's triangle and Hang Si's tetgen
#define REAL double
#define VOID void
#include "triangle.h"
#include "tetgen.h"
#undef REAL
#undef VOID

template<unsigned DIM>
void QuadraticMesh<DIM>::CountAndCheckVertices()
{
    // count the number of vertices, and also check all vertices come before the
    // rest of the nodes (as this is assumed in
    // AbstractNonlinearElasticitySolver<DIM>::AllocateMatrixMemory() )
    //
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
    this->mNodes.resize(2*numElemX+1);
    mNumVertices = numElemX+1;

    // create the left-most node
    Node<DIM>* p_edge_node = new Node<DIM>(0, true, 0.0);
    this->mNodes[0] = p_edge_node; // create nodes
    this->mBoundaryNodes.push_back(p_edge_node);
    this->mBoundaryElements.push_back(new BoundaryElement<DIM-1,DIM>(0, p_edge_node) );

    for (unsigned element_index=0; element_index<numElemX; element_index++)
    {
        unsigned right_node_index = element_index+1;
        unsigned mid_node_index = mNumVertices + element_index;

        double x_value_right_x = right_node_index;
        double x_value_mid_node = x_value_right_x-0.5;

        bool is_boundary = (element_index+1==numElemX);
        Node<DIM>* p_right_node = new Node<DIM>(right_node_index, is_boundary, x_value_right_x);
        Node<DIM>* p_mid_node   = new Node<DIM>(mid_node_index, false, x_value_mid_node);

        this->mNodes[right_node_index] = p_right_node;
        this->mNodes[mid_node_index] = p_mid_node;

        if (element_index+1==numElemX) // right boundary
        {
            this->mBoundaryNodes.push_back(p_right_node);
            this->mBoundaryElements.push_back(new BoundaryElement<DIM-1,DIM>(1, p_right_node) );
        }

        std::vector<Node<DIM>*> nodes;
        nodes.push_back(this->mNodes[right_node_index-1]);
        nodes.push_back(this->mNodes[right_node_index]);
        nodes.push_back(this->mNodes[mid_node_index]);
        this->mElements.push_back(new Element<DIM,DIM>(element_index, nodes) );
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

    this->mMeshIsLinear=false;
    unsigned num_nodes=(numElemX+1)*(numElemY+1);
    struct triangulateio mesher_input;

    this->InitialiseTriangulateIo(mesher_input);
    mesher_input.pointlist = (double *) malloc( num_nodes * DIM * sizeof(double));
    mesher_input.numberofpoints = num_nodes;

    unsigned new_index = 0;
    for (unsigned j=0; j<=numElemY; j++)
    {
        double y = j;
        for (unsigned i=0; i<=numElemX; i++)
        {
            double x = i;

            mesher_input.pointlist[DIM*new_index] = x;
            mesher_input.pointlist[DIM*new_index + 1] = y;
            new_index++;
        }
    }

    // Make structure for output
    struct triangulateio mesher_output;

    this->InitialiseTriangulateIo(mesher_output);

    // Library call
    triangulate((char*)"Qzeo2", &mesher_input, &mesher_output, NULL);

    assert(mesher_output.numberofcorners == (DIM+1)*(DIM+2)/2);//Nodes per element (including internals, one per edge)

    this->ImportFromMesher(mesher_output, mesher_output.numberoftriangles, mesher_output.trianglelist, mesher_output.numberofedges, mesher_output.edgelist, mesher_output.edgemarkerlist);

    CountAndCheckVertices();
    QuadraticMeshHelper<DIM>::AddNodesToBoundaryElements(this, NULL);

    this->FreeTriangulateIo(mesher_input);
    this->FreeTriangulateIo(mesher_output);
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

    CountAndCheckVertices();
    QuadraticMeshHelper<DIM>::AddNodesToBoundaryElements(this, NULL);
}


template<unsigned DIM>
unsigned QuadraticMesh<DIM>::GetNumVertices()
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
        CountAndCheckVertices();
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
        CountAndCheckVertices();
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
    CountAndCheckVertices();
    QuadraticMeshHelper<DIM>::AddInternalNodesToBoundaryElements(this, p_mesh_reader);
    QuadraticMeshHelper<DIM>::CheckBoundaryElements(this);
}


template<unsigned DIM>
void QuadraticMesh<DIM>::WriteBoundaryElementFile(std::string directory, std::string fileName)
{
    OutputFileHandler handler(directory, false);
    out_stream p_file = handler.OpenOutputFile(fileName);

    unsigned expected_num_nodes;
    assert(DIM > 1);
    if (DIM == 2)
    {
        expected_num_nodes = 3;
    }
    else if (DIM == 3)
    {
        expected_num_nodes = 6;
    }

    unsigned num_elements = 0;

    for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
          = this->GetBoundaryElementIteratorBegin();
          iter != this->GetBoundaryElementIteratorEnd();
          ++iter)
    {
        assert((*iter)->GetNumNodes()==expected_num_nodes);
        num_elements++;
    }

    *p_file << num_elements << " 0\n";

    unsigned counter = 0;
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
          = this->GetBoundaryElementIteratorBegin();
          iter != this->GetBoundaryElementIteratorEnd();
          ++iter)
    {
        *p_file << counter++ << " ";
        for (unsigned i=0; i<(*iter)->GetNumNodes(); i++)
        {
            *p_file << (*iter)->GetNodeGlobalIndex(i) << " ";
        }
        *p_file << "\n";
    }

    p_file->close();
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
