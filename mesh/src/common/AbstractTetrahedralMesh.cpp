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

#include <limits>
#include "AbstractTetrahedralMesh.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships()
{
    unsigned lo=this->GetDistributedVectorFactory()->GetLow();
    unsigned hi=this->GetDistributedVectorFactory()->GetHigh();
    for (unsigned element_index=0; element_index<mElements.size(); element_index++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element = mElements[element_index];
        p_element->SetOwnership(false);
        for (unsigned local_node_index=0; local_node_index< p_element->GetNumNodes(); local_node_index++)
        {
            unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
            if (lo<=global_node_index && global_node_index<hi)
            {
                p_element->SetOwnership(true);
                break;
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralMesh()
    : mMeshIsLinear(true)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~AbstractTetrahedralMesh()
{
    // Iterate over elements and free the memory
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    // Iterate over boundary elements and free the memory
    for (unsigned i=0; i<mBoundaryElements.size(); i++)
    {
        delete mBoundaryElements[i];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalElements() const
{
    return GetNumElements();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
{
    return mBoundaryElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalBoundaryElements() const
{
    return GetNumBoundaryElements();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllBoundaryElements() const
{
    return mBoundaryElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumCableElements() const
{
    return 0u;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumVertices() const
{
    return this->GetNumNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetMaximumNodeIndex()
{
    return this->GetNumAllNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    unsigned local_index = SolveElementMapping(index);
    return mElements[local_index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElement(unsigned index) const
{
    unsigned local_index = SolveBoundaryElementMapping(index);
    return mBoundaryElements[local_index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElementIteratorBegin() const
{
    return mBoundaryElements.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElementIteratorEnd() const
{
    return mBoundaryElements.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetInverseJacobianForElement(
        unsigned elementIndex,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian,
        double& rJacobianDeterminant,
        c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const
{
    mElements[SolveElementMapping(elementIndex)]->CalculateInverseJacobian(rJacobian, rJacobianDeterminant, rInverseJacobian);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForBoundaryElement(
        unsigned elementIndex,
        c_vector<double, SPACE_DIM>& rWeightedDirection,
        double& rJacobianDeterminant) const
{
    mBoundaryElements[SolveBoundaryElementMapping(elementIndex)]->CalculateWeightedDirection(rWeightedDirection, rJacobianDeterminant );
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CheckOutwardNormals()
{
    if (ELEMENT_DIM <= 1)
    {
        //If the ELEMENT_DIM of the mesh is 1 then the boundary will have ELEMENT_DIM = 0
        EXCEPTION("1-D mesh has no boundary normals");
    }
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator face_iter = this->GetBoundaryElementIteratorBegin();
         face_iter != this->GetBoundaryElementIteratorEnd();
         ++face_iter)
    {
        //Form a set for the boundary element node indices
        std::set<unsigned> boundary_element_node_indices;
        for (unsigned i=0; i<ELEMENT_DIM; i++)
        {
            boundary_element_node_indices.insert( (*face_iter)->GetNodeGlobalIndex(i) );
        }

        Node<SPACE_DIM>* p_opposite_node = nullptr;
        Node<SPACE_DIM>* p_representative_node = (*face_iter)->GetNode(0);
          for (typename Node<SPACE_DIM>::ContainingElementIterator element_iter = p_representative_node->ContainingElementsBegin();
             element_iter != p_representative_node->ContainingElementsEnd();
             ++element_iter)
        {
            Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*element_iter);
            //Form a set for the element node indices
            std::set<unsigned> element_node_indices;
            for (unsigned i=0; i<=ELEMENT_DIM; i++)
            {
                element_node_indices.insert( p_element->GetNodeGlobalIndex(i) );
            }

            std::vector<unsigned> difference(ELEMENT_DIM);

            std::vector<unsigned>::iterator set_iter = std::set_difference(
                    element_node_indices.begin(),element_node_indices.end(),
                    boundary_element_node_indices.begin(), boundary_element_node_indices.end(),
                    difference.begin());
            if (set_iter - difference.begin() == 1)
            {
                p_opposite_node = this -> GetNodeOrHaloNode(difference[0]);
                break;
            }
        }
        assert(p_opposite_node != nullptr);

        // Vector from centroid of face to opposite node
        c_vector<double, SPACE_DIM> into_mesh = p_opposite_node->rGetLocation() - (*face_iter)->CalculateCentroid();
        c_vector<double, SPACE_DIM> normal = (*face_iter)->CalculateNormal();

        if (inner_prod(into_mesh, normal) > 0.0)
        {
            EXCEPTION("Inward facing normal in boundary element index "<<(*face_iter)->GetIndex());
        }
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructLinearMesh(unsigned width)
{
    assert(ELEMENT_DIM == 1);     // LCOV_EXCL_LINE

    for (unsigned node_index=0; node_index<=width; node_index++)
    {
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, node_index==0 || node_index==width, node_index);
        this->mNodes.push_back(p_node); // create node
        if (node_index==0) // create left boundary node and boundary element
        {
            this->mBoundaryNodes.push_back(p_node);
            this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(0, p_node) );
        }
        if (node_index==width) // create right boundary node and boundary element
        {
            this->mBoundaryNodes.push_back(p_node);
            this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(1, p_node) );
        }
        if (node_index > 0) // create element
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            nodes.push_back(this->mNodes[node_index-1]);
            nodes.push_back(this->mNodes[node_index]);
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(node_index-1, nodes) );
        }
    }

    this->RefreshMesh();
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRectangularMesh(unsigned width, unsigned height, bool stagger)
{
    assert(SPACE_DIM == 2);     // LCOV_EXCL_LINE
    assert(ELEMENT_DIM == 2);     // LCOV_EXCL_LINE

    //Construct the nodes
    unsigned node_index=0;
    for (unsigned j=0; j<height+1; j++)
    {
        for (unsigned i=0; i<width+1; i++)
        {
            bool is_boundary=false;
            if (i==0 || j==0 || i==width || j==height)
            {
                is_boundary=true;
            }
            //Check in place for parallel
            assert(node_index==(width+1)*(j) + i);
            Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index++, is_boundary, i, j);
            this->mNodes.push_back(p_node);
            if (is_boundary)
            {
                this->mBoundaryNodes.push_back(p_node);
            }
        }
    }

    //Construct the boundary elements
    unsigned belem_index=0;
    //Top
    for (unsigned i=0; i<width; i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[height*(width+1)+i+1]);
        nodes.push_back(this->mNodes[height*(width+1)+i]);
        assert(belem_index==i);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Right
    for (unsigned j=1; j<=height; j++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[(width+1)*j-1]);
        nodes.push_back(this->mNodes[(width+1)*(j+1)-1]);
        assert(belem_index==width+j-1);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Bottom
    for (unsigned i=0; i<width; i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[i]);
        nodes.push_back(this->mNodes[i+1]);
        assert(belem_index==width+height+i);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Left
    for (unsigned j=0; j<height; j++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[(width+1)*(j+1)]);
        nodes.push_back(this->mNodes[(width+1)*(j)]);
        assert(belem_index==2*width+height+j);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }

    //Construct the elements
    unsigned elem_index = 0;
    for (unsigned j=0; j<height; j++)
    {
        for (unsigned i=0; i<width; i++)
        {
            unsigned parity=(i+(height-j))%2;//Note that parity is measured from the top-left (not bottom left) for historical reasons
            unsigned nw=(j+1)*(width+1)+i; //ne=nw+1
            unsigned sw=(j)*(width+1)+i;   //se=sw+1
            std::vector<Node<SPACE_DIM>*> upper_nodes;
            upper_nodes.push_back(this->mNodes[nw]);
            upper_nodes.push_back(this->mNodes[nw+1]);
            if (stagger==false  || parity == 1)
            {
                upper_nodes.push_back(this->mNodes[sw+1]);
            }
            else
            {
                upper_nodes.push_back(this->mNodes[sw]);
            }
            assert(elem_index==2*(j*width+i));
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++,upper_nodes));
            std::vector<Node<SPACE_DIM>*> lower_nodes;
            lower_nodes.push_back(this->mNodes[sw+1]);
            lower_nodes.push_back(this->mNodes[sw]);
            if (stagger==false  ||parity == 1)
            {
                lower_nodes.push_back(this->mNodes[nw]);
            }
            else
            {
                lower_nodes.push_back(this->mNodes[nw+1]);
            }
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++,lower_nodes));
        }
    }

    this->RefreshMesh();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructCuboid(unsigned width,
        unsigned height,
        unsigned depth)
{
    assert(SPACE_DIM == 3);     // LCOV_EXCL_LINE
    assert(ELEMENT_DIM == 3);     // LCOV_EXCL_LINE
    //Construct the nodes

    unsigned node_index = 0;
    for (unsigned k=0; k<depth+1; k++)
    {
        for (unsigned j=0; j<height+1; j++)
        {
            for (unsigned i=0; i<width+1; i++)
            {
                bool is_boundary = false;
                if (i==0 || j==0 || k==0 || i==width || j==height || k==depth)
                {
                    is_boundary = true;
                }
                assert(node_index == (k*(height+1)+j)*(width+1)+i);
                Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index++, is_boundary, i, j, k);

                this->mNodes.push_back(p_node);
                if (is_boundary)
                {
                    this->mBoundaryNodes.push_back(p_node);
                }
            }
        }
    }

    // Construct the elements

    unsigned elem_index = 0;
    unsigned belem_index = 0;
    unsigned element_nodes[6][4] = {{0, 1, 5, 7}, {0, 1, 3, 7},
                                        {0, 2, 3, 7}, {0, 2, 6, 7},
                                        {0, 4, 6, 7}, {0, 4, 5, 7}};
/* Alternative tessellation - (gerardus)
 * Note that our method (above) has a bias that all tetrahedra share a
 * common edge (the diagonal 0 - 7).  In the following method the cube is
 * split along the "face diagonal" 1-2-5-6 into two prisms.  This also has a bias.
 *
    unsigned element_nodes[6][4] = {{ 0, 6, 5, 4},
                                    { 0, 2, 6, 1},
                                    { 0, 1, 6, 5},
                                    { 1, 2, 3, 7},
                                    { 1, 2, 6, 7},
                                    { 1, 6, 7, 5 }};
*/
    std::vector<Node<SPACE_DIM>*> tetrahedra_nodes;

    for (unsigned k=0; k<depth; k++)
    {
        if (k!=0)
        {
            // height*width squares on upper face, k layers of 2*height+2*width square aroun
            assert(belem_index ==   2*(height*width+k*2*(height+width)) );
        }
        for (unsigned j=0; j<height; j++)
        {
            for (unsigned i=0; i<width; i++)
            {
                // Compute the nodes' index
                unsigned global_node_indices[8];
                unsigned local_node_index = 0;

                for (unsigned z = 0; z < 2; z++)
                {
                    for (unsigned y = 0; y < 2; y++)
                    {
                        for (unsigned x = 0; x < 2; x++)
                        {
                            global_node_indices[local_node_index] = i+x+(width+1)*(j+y+(height+1)*(k+z));

                            local_node_index++;
                        }
                    }
                }

                for (unsigned m = 0; m < 6; m++)
                {
                    // Tetrahedra #m

                    tetrahedra_nodes.clear();

                    for (unsigned n = 0; n < 4; n++)
                    {
                        tetrahedra_nodes.push_back(this->mNodes[global_node_indices[element_nodes[m][n]]]);
                    }

                    assert(elem_index == 6 * ((k*height+j)*width+i)+m );
                    this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++, tetrahedra_nodes));
                }

                // Are we at a boundary?
                std::vector<Node<SPACE_DIM>*> triangle_nodes;

                if (i == 0) //low face at x==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[6]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[6]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[4]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (i == width-1) //high face at x=width
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[5]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[3]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (j == 0) //low face at y==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[5]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[1]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[5]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (j == height-1) //high face at y=height
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[3]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[6]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (k == 0) //low face at z==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[3]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[2]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[3]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (k == depth-1) //high face at z=depth
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[5]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[6]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
            }//i
        }//j
    }//k

    this->RefreshMesh();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRegularSlabMesh(double spaceStep, double width, double height, double depth)
{
    assert(spaceStep>0.0);
    assert(width>0.0);
    if (ELEMENT_DIM > 1)
    {
        assert(height>0.0);
    }
    if (ELEMENT_DIM > 2)
    {
        assert(depth>0.0);
    }
    unsigned num_elem_x=(unsigned)((width+0.5*spaceStep)/spaceStep); //0.5*spaceStep is to ensure that rounding down snaps to correct number
    unsigned num_elem_y=(unsigned)((height+0.5*spaceStep)/spaceStep);
    unsigned num_elem_z=(unsigned)((depth+0.5*spaceStep)/spaceStep);

    //Make it obvious that actual_width_x etc. are temporaries used in spotting for exception
    {
        double actual_width_x=num_elem_x*spaceStep;
        double actual_width_y=num_elem_y*spaceStep;
        double actual_width_z=num_elem_z*spaceStep;

        //Note here that in ELEMENT_DIM > 1 cases there may be a zero height or depth - in which case we don't need to use relative comparisons
        // Doing relative comparisons with zero is okay - if we avoid division by zero.
        // However, it's best not to test whether " fabs( 0.0 - 0.0) > DBL_EPSILON*0.0 "
        if (fabs (actual_width_x - width) > DBL_EPSILON*width
            ||( height!= 0.0 &&  fabs (actual_width_y - height) > DBL_EPSILON*height)
            ||( depth != 0.0 &&  fabs (actual_width_z - depth) > DBL_EPSILON*depth ))
        {
            EXCEPTION("Space step does not divide the size of the mesh");
        }
    }
    switch (ELEMENT_DIM)
    {
        case 1:
            this->ConstructLinearMesh(num_elem_x);
            break;
        case 2:
            this->ConstructRectangularMesh(num_elem_x, num_elem_y); // Stagger=default value
            break;
        default:
        case 3:
            this->ConstructCuboid(num_elem_x, num_elem_y, num_elem_z);
    }
    this->Scale(spaceStep, spaceStep, spaceStep);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRegularSlabMeshWithDimensionSplit(
        unsigned dimension, double spaceStep,
        double width, double height, double depth)
{
    assert(ELEMENT_DIM == SPACE_DIM);     // LCOV_EXCL_LINE
    if (dimension >= SPACE_DIM)
    {
        EXCEPTION("Cannot split on non-existent dimension");
    }
    // Rotate the width -> height -> depth around (if dimension is non-default)
    if (SPACE_DIM == 2 && dimension == 0)
    {
        double temp = height ;
        height = width;
        width = temp;
    }
    else if (SPACE_DIM == 3)
    {
        unsigned rotate_perm = SPACE_DIM - 1u - dimension; // How many shuffles to get the user's axis to the top
        for (unsigned i=0; i<rotate_perm; i++)
        {
            double temp = depth;
            depth = height;
            height = width;
            width = temp;
        }
    }
    this->ConstructRegularSlabMesh(spaceStep, width, height, depth);
    if (SPACE_DIM == 2 && dimension == 0)
    {
        // Rotate the positions back again x -> y -> x
        // this->Rotate(M_PI_2);
        c_matrix<double, 2, 2> axis_rotation = zero_matrix<double>(2, 2);
        axis_rotation(0,1)=1.0;
        axis_rotation(1,0)=-1.0;
        this->Rotate(axis_rotation);
        this->Translate(0.0, width); // Formerly known as height, but we rotated it
    }
    else if (SPACE_DIM == 3 && dimension == 0)
    {
        //        this->RotateZ(M_PI_2);
        //        this->RotateY(M_PI_2);
        // RotY * RotZ = [0 0 1; 1 0 0; 0 1 0] x->y->z->x
        //this->Translate(depth /*old width*/, width /*old height*/, 0.0);
        c_matrix<double, 3, 3> axis_permutation = zero_matrix<double>(3, 3);
        axis_permutation(0, 2)=1.0;
        axis_permutation(1, 0)=1.0;
        axis_permutation(2, 1)=1.0;
        this->Rotate(axis_permutation);
    }
    else if (SPACE_DIM == 3 && dimension == 1)
    {
        //        this->RotateY(-M_PI_2);
        //        this->RotateZ(-M_PI_2);
        //        // RotZ' after RotY' = RotZ' * RotY' = [0 1 0; 0 0 1; 1 0 0] x->z->y->x
        //        this->Translate(height /*old width*/, 0.0, width /*old depth*/);
        c_matrix<double, 3, 3> axis_permutation = zero_matrix<double>(3, 3);
        axis_permutation(0, 1)=1.0;
        axis_permutation(1, 2)=1.0;
        axis_permutation(2, 0)=1.0;
        this->Rotate(axis_permutation);
    }

}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateDesignatedOwnershipOfBoundaryElement( unsigned faceIndex )
{
    // This may throw in the distributed parallel case
    unsigned tie_break_index = this->GetBoundaryElement(faceIndex)->GetNodeGlobalIndex(0);

    // if it is in my range
    if (this->GetDistributedVectorFactory()->IsGlobalIndexLocal(tie_break_index))
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateDesignatedOwnershipOfElement( unsigned elementIndex )
{
    // This may throw in the distributed parallel case
    unsigned tie_break_index = this->GetElement(elementIndex)->GetNodeGlobalIndex(0);

    // if it is in my range
    if (this->GetDistributedVectorFactory()->IsGlobalIndexLocal(tie_break_index))
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMaximumNodeConnectivityPerProcess() const
{

    if (this->mNodes.size() == 0u)
    {
// LCOV_EXCL_START
        /*
         * Coverage of this block requires a mesh regular slab mesh with the number of
         * elements in the primary dimension less than (num_procs - 1), e.g. a 1D mesh
         * one element wide with num_procs >=3.
         * Note that if a process owns no nodes, then it will still need to enter the
         * collective call to MatMPIAIJSetPreallocation. In PetscTools::SetupMat, the
         * rowPreallocation parameter uses the special value zero to indicate no preallocation.
         */

        // This process owns no nodes and thus owns none of the mesh
        return (1u);
// LCOV_EXCL_STOP
    }

    unsigned nodes_per_element = this->mElements[0]->GetNumNodes(); //Usually ELEMENT_DIM+1, except in Quadratic case
    if (ELEMENT_DIM <= 2u)
    {
        /*
         *  Note that we start assuming that each internal node is connected to 1 element.
         *  This is to avoid the trivial situation in which there are no internal nodes (a
         *  single line or a single triangle.  We need the minimum connectivity to be 2 (in 1D) or
         *  3 (in 2D) even if there is only one element.
         */
        unsigned max_num = 1u; // See note above.
        unsigned boundary_max_num = 0u;
        for (unsigned local_node_index=0; local_node_index<this->mNodes.size(); local_node_index++)
        {
            unsigned num = this->mNodes[local_node_index]->GetNumContainingElements();
            if (this->mNodes[local_node_index]->IsBoundaryNode()==false && num>max_num)
            {
                max_num = num;
            }
            if (this->mNodes[local_node_index]->IsBoundaryNode() && num>boundary_max_num)
            {
                boundary_max_num = num;
            }
        }
        bool linear = (nodes_per_element == ELEMENT_DIM + 1);
        /*
         * In 1d each containing element is connected to one node (or 2 if quadratic), add to this the node itself
         * and the connectivity is GetNumContainingElements() + 1 or 2*GetNumContainingElements() + 1
         */
        if (ELEMENT_DIM == 1)
        {
            if (linear)
            {
                return max_num+1;
            }
            else
            {
                return 2*max_num+1;
            }
        }
        // Not returned ...else if  (ELEMENT_DIM == 2)
        /*
         * In 2d each containing element adds one connected node (since one node will be shared by a previous element)
         * this leads to a connectivity of GetNumContainingElements() + 1 (or 3*GetNumContainingElements() + 1) in the quadratic case
         *
         * If the node is on a boundary then the one elements will have an unpaired node and the connectivity is
         * GetNumContainingElements() + 1 (or 3*(GetNumContainingElements() + 3 for quadratic)
         */
        if (linear)
        {
            return std::max(max_num+1, boundary_max_num+2);
        }
        else
        {
            return std::max(3*max_num+1, 3*boundary_max_num+3);
        }
    }

    /*
     * In 3d there are many more cases.  In general a non-boundary node has fewer connecting nodes than it has elements.
     * A node on the boundary may have even fewer, unless it is on a corner and has more faces than it has elements.
     * We can, in the linear case estimate an upper bound as max(elements, faces)+2.
     * However for the sake of accuracy we are going for a brute force solution.
     *
     * This may prove to be a bottle-neck...
     */

    std::set<unsigned> forward_star_nodes; // Used to collect each node's neighbours
    unsigned max_connectivity = 0u;
    for (unsigned local_node_index=0; local_node_index<this->mNodes.size(); local_node_index++)
    {
        forward_star_nodes.clear();

        for (typename Node<SPACE_DIM>::ContainingElementIterator it = this->mNodes[local_node_index]->ContainingElementsBegin();
                it != this->mNodes[local_node_index]->ContainingElementsEnd();
                ++it)
        {
            Element<ELEMENT_DIM, SPACE_DIM>* p_elem = this->GetElement(*it);
            for (unsigned i=0; i<nodes_per_element; i++)
            {
                forward_star_nodes.insert(p_elem->GetNodeGlobalIndex(i));
            }
        }
        if (forward_star_nodes.size() > max_connectivity)
        {
            max_connectivity = forward_star_nodes.size();
        }
    }
    return max_connectivity;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetHaloNodeIndices(std::vector<unsigned>& rHaloIndices) const
{
    // Make sure the output vector is empty
    rHaloIndices.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMesh(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rOtherMesh)
{
    for (unsigned i=0; i<rOtherMesh.GetNumNodes(); i++)
    {
        Node<SPACE_DIM>* p_node =  rOtherMesh.GetNode(i);
        assert(!p_node->IsDeleted());
        const c_vector<double, SPACE_DIM>& location=p_node->rGetLocation();
        bool is_boundary=p_node->IsBoundaryNode();

        Node<SPACE_DIM>* p_node_copy = new Node<SPACE_DIM>(i, location, is_boundary);
        this->mNodes.push_back( p_node_copy );
        if (is_boundary)
        {
            this->mBoundaryNodes.push_back( p_node_copy );
        }
    }

    for (unsigned i=0; i<rOtherMesh.GetNumElements(); i++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_elem = rOtherMesh.GetElement(i);
        assert(!p_elem->IsDeleted());
        std::vector<Node<SPACE_DIM>*> nodes_for_element;
        for (unsigned j=0; j<p_elem->GetNumNodes(); j++)
        {
            nodes_for_element.push_back(this->mNodes[ p_elem->GetNodeGlobalIndex(j) ]);
        }
        Element<ELEMENT_DIM, SPACE_DIM>* p_elem_copy = new Element<ELEMENT_DIM, SPACE_DIM>(i, nodes_for_element);
        p_elem_copy->RegisterWithNodes();
        this->mElements.push_back(p_elem_copy);
    }

    for (unsigned i=0; i<rOtherMesh.GetNumBoundaryElements(); i++)
    {
        BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_b_elem =  rOtherMesh.GetBoundaryElement(i);
        assert(!p_b_elem->IsDeleted());
        std::vector<Node<SPACE_DIM>*> nodes_for_element;
        for (unsigned j=0; j<p_b_elem->GetNumNodes(); j++)
        {
            nodes_for_element.push_back(this->mNodes[ p_b_elem->GetNodeGlobalIndex(j) ]);
        }
        BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_b_elem_copy = new BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>(i, nodes_for_element);
        p_b_elem_copy->RegisterWithNodes();
        this->mBoundaryElements.push_back(p_b_elem_copy);
    }
    this->RefreshMesh();

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateNodeExchange(
                                     std::vector<std::vector<unsigned> >& rNodesToSendPerProcess,
                                     std::vector<std::vector<unsigned> >& rNodesToReceivePerProcess)
{
    assert( rNodesToSendPerProcess.empty() );
    assert( rNodesToReceivePerProcess.empty() );

    //Initialise vectors of sets for the exchange data
    std::vector<std::set<unsigned> > node_sets_to_send_per_process;
    std::vector<std::set<unsigned> > node_sets_to_receive_per_process;

    node_sets_to_send_per_process.resize(PetscTools::GetNumProcs());
    node_sets_to_receive_per_process.resize(PetscTools::GetNumProcs());
    std::vector<unsigned> global_lows = this->GetDistributedVectorFactory()->rGetGlobalLows();

    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        std::vector <unsigned> nodes_on_this_process;
        std::vector <unsigned> nodes_not_on_this_process;
        //Calculate local and non-local node indices
        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            unsigned node_index=iter->GetNodeGlobalIndex(i);
            if (this->GetDistributedVectorFactory()->IsGlobalIndexLocal(node_index))
            {
                nodes_on_this_process.push_back(node_index);
            }
            else
            {
                nodes_not_on_this_process.push_back(node_index);
            }
        }
        /*
         * If this is a TetrahedralMesh (not distributed) then it's possible that we own none
         * of the nodes in this element.  In that case we must skip the element.
         */
        if (nodes_on_this_process.empty())
        {
            continue; //Move on to the next element.
            ///\todo #2426 This is where we should check that we DO HAVE a Tetrahedral (not Distributed mesh)
        }
        // If there are any non-local nodes on this element then we need to add to the data exchange
        if (!nodes_not_on_this_process.empty())
        {
            for (unsigned i=0; i<nodes_not_on_this_process.size(); i++)
            {
                // Calculate who owns this remote node
                unsigned remote_process=global_lows.size()-1;
                for (; global_lows[remote_process] > nodes_not_on_this_process[i]; remote_process--)
                {
                }

                // Add this node to the correct receive set
                node_sets_to_receive_per_process[remote_process].insert(nodes_not_on_this_process[i]);

                // Add all local nodes to the send set
                for (unsigned j=0; j<nodes_on_this_process.size(); j++)
                {
                    node_sets_to_send_per_process[remote_process].insert(nodes_on_this_process[j]);
                }
            }
        }
    }

    for (unsigned process_number = 0; process_number < PetscTools::GetNumProcs(); process_number++)
    {
        std::vector<unsigned> process_send_vector( node_sets_to_send_per_process[process_number].begin(),
                                                   node_sets_to_send_per_process[process_number].end()    );
        std::sort(process_send_vector.begin(), process_send_vector.end());

        rNodesToSendPerProcess.push_back(process_send_vector);

        std::vector<unsigned> process_receive_vector( node_sets_to_receive_per_process[process_number].begin(),
                                                      node_sets_to_receive_per_process[process_number].end()    );
        std::sort(process_receive_vector.begin(), process_receive_vector.end());

        rNodesToReceivePerProcess.push_back(process_receive_vector);
    }

}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 2> AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMinMaxEdgeLengths()
{
    c_vector<double, 2> min_max;
    min_max[0] = DBL_MAX;
    min_max[1] = 0.0;
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator ele_iter = GetElementIteratorBegin();
            ele_iter != GetElementIteratorEnd();
            ++ele_iter)
    {
        c_vector<double, 2> ele_min_max = ele_iter->CalculateMinMaxEdgeLengths();
        if (ele_min_max[0] < min_max[0])
        {
            min_max[0] = ele_min_max[0];
        }
        if (ele_min_max[1] > min_max[1])
        {
            min_max[1] = ele_min_max[1];
        }
    }
    return min_max;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndex(const ChastePoint<SPACE_DIM>& rTestPoint,
                                                                            bool strict,
                                                                            std::set<unsigned> testElements,
                                                                            bool onlyTryWithTestElements)
{
    for (std::set<unsigned>::iterator iter=testElements.begin(); iter!=testElements.end(); iter++)
    {
        assert(*iter<this->GetNumElements());
        if (this->mElements[*iter]->IncludesPoint(rTestPoint, strict))
        {
            assert(!this->mElements[*iter]->IsDeleted());
            return *iter;
        }
    }

    if (!onlyTryWithTestElements)
    {
        for (unsigned i=0; i<this->mElements.size(); i++)
        {
            if (this->mElements[i]->IncludesPoint(rTestPoint, strict))
            {
                assert(!this->mElements[i]->IsDeleted());
                return i;
            }
        }
    }

    // If it's in none of the elements, then throw
    std::stringstream ss;
    ss << "Point [";
    for (unsigned j=0; (int)j<(int)SPACE_DIM-1; j++)
    {
        ss << rTestPoint[j] << ",";
    }
    ss << rTestPoint[SPACE_DIM-1] << "] is not in ";
    if (!onlyTryWithTestElements)
    {
        ss << "mesh - all elements tested";
    }
    else
    {
        ss << "set of elements given";
    }
    EXCEPTION(ss.str());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNearestElementIndexFromTestElements(const ChastePoint<SPACE_DIM>& rTestPoint,
                                                                                                 std::set<unsigned> testElements)
{
    assert(testElements.size() > 0);
    EXCEPT_IF_NOT(ELEMENT_DIM == SPACE_DIM); // LCOV_EXCL_LINE // CalculateInterpolationWeights hits an assertion otherwise

    double max_min_weight = -std::numeric_limits<double>::infinity();
    unsigned closest_index = 0;
    for (std::set<unsigned>::iterator iter = testElements.begin();
        iter != testElements.end();
        iter++)
    {
        c_vector<double, ELEMENT_DIM+1> weight = this->mElements[*iter]->CalculateInterpolationWeights(rTestPoint);
        double neg_weight_sum = 0.0;
        for (unsigned j=0; j<=ELEMENT_DIM; j++)
        {
            if (weight[j] < 0.0)
            {
                neg_weight_sum += weight[j];
            }
        }
        if (neg_weight_sum > max_min_weight)
        {
            max_min_weight = neg_weight_sum;
            closest_index = *iter;
        }
    }
    assert(!this->mElements[closest_index]->IsDeleted());
    return closest_index;
}

// Explicit instantiation
template class AbstractTetrahedralMesh<1,1>;
template class AbstractTetrahedralMesh<1,2>;
template class AbstractTetrahedralMesh<1,3>;
template class AbstractTetrahedralMesh<2,2>;
template class AbstractTetrahedralMesh<2,3>;
template class AbstractTetrahedralMesh<3,3>;
