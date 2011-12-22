/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include <map>
#include <cstring>

#include "MutableMesh.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"

//Jonathan Shewchuk's triangle and Hang Si's tetgen
#define REAL double
#define VOID void
#include "triangle.h"
#include "tetgen.h"
#undef REAL
#undef VOID


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMesh<ELEMENT_DIM, SPACE_DIM>::MutableMesh()
    : mAddedNodes(false)
{
    this->mMeshChangesDuringSimulation = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMesh<ELEMENT_DIM, SPACE_DIM>::MutableMesh(std::vector<Node<SPACE_DIM> *> nodes)
{
    this->mMeshChangesDuringSimulation = true;
    Clear();
    for (unsigned index=0; index<nodes.size(); index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[index];
        this->mNodes.push_back(p_temp_node);
    }
    mAddedNodes = true;
    NodeMap node_map(nodes.size());
    ReMesh(node_map);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMesh<ELEMENT_DIM, SPACE_DIM>::~MutableMesh()
{
    Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
{
    if (mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(this->mNodes.size());
        this->mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index = mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        mDeletedNodeIndices.pop_back();
        delete this->mNodes[index];
        this->mNodes[index] = pNewNode;
    }
    mAddedNodes = true;
    return pNewNode->GetIndex();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    mDeletedElementIndices.clear();
    mDeletedBoundaryElementIndices.clear();
    mDeletedNodeIndices.clear();
    mAddedNodes = false;

    TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
{
    return this->mBoundaryElements.size() - mDeletedBoundaryElementIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return this->mElements.size() - mDeletedElementIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size() - mDeletedNodeIndices.size();
}

/**
 * The RescaleMeshFromBoundaryNode method is only defined for 1D meshes.
 *
 * @param updatedPoint point determining the scale factor
 * @param boundaryNodeIndex index of the boundary node
 */
template<>
void MutableMesh<1, 1>::RescaleMeshFromBoundaryNode(ChastePoint<1> updatedPoint, unsigned boundaryNodeIndex)
{
    assert(this->GetNode(boundaryNodeIndex)->IsBoundaryNode());
    double scaleFactor = updatedPoint[0] / this->GetNode(boundaryNodeIndex)->GetPoint()[0];
    double temp;
    for (unsigned i=0; i < boundaryNodeIndex+1; i++)
    {
        temp = scaleFactor * this->mNodes[i]->GetPoint()[0];
        ChastePoint<1> newPoint(temp);
        this->mNodes[i]->SetPoint(newPoint);
    }
    this->RefreshMesh();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned index,
        ChastePoint<SPACE_DIM> point,
        bool concreteMove)
{
    this->mNodes[index]->SetPoint(point);

    if (concreteMove)
    {
        for (typename Node<SPACE_DIM>::ContainingElementIterator it = this->mNodes[index]->ContainingElementsBegin();
             it != this->mNodes[index]->ContainingElementsEnd();
             ++it)
        {
            if (ELEMENT_DIM == SPACE_DIM)
            {
                try
                {
                    GetElement(*it)->CalculateInverseJacobian(this->mElementJacobians[ (*it) ],
                                                              this->mElementJacobianDeterminants[ (*it) ],
                                                              this->mElementInverseJacobians[ (*it) ]);
                }
                catch (Exception e)
                {
                        EXCEPTION("Moving node caused an element to have a non-positive Jacobian determinant");
                }
            }
            else
            {
                c_vector<double,SPACE_DIM> previous_direction = this->mElementWeightedDirections[ (*it) ];

                GetElement(*it)->CalculateWeightedDirection(this->mElementWeightedDirections[ (*it) ],
                                                            this->mElementJacobianDeterminants[ (*it) ]);

                if ( inner_prod(previous_direction, this->mElementWeightedDirections[ (*it) ]) < 0)
                {
                    EXCEPTION("Moving node caused an subspace element to change direction");
                }

            }
        }
        for (typename Node<SPACE_DIM>::ContainingBoundaryElementIterator it = this->mNodes[index]->ContainingBoundaryElementsBegin();
             it != this->mNodes[index]->ContainingBoundaryElementsEnd();
             ++it)
        {
            try
            {
                GetBoundaryElement(*it)->CalculateWeightedDirection(this->mBoundaryElementWeightedDirections[ (*it) ],
                                                                    this->mBoundaryElementJacobianDeterminants[ (*it) ]);
            }
            catch (Exception e)
            {
                EXCEPTION("Moving node caused a boundary element to have a non-positive Jacobian determinant");
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::DeleteNode(unsigned index)
{
    if (this->mNodes[index]->IsDeleted())
    {
        EXCEPTION("Trying to delete a deleted node");
    }
    unsigned target_index = (unsigned)(-1);
    bool found_target=false;
    for (typename Node<SPACE_DIM>::ContainingElementIterator it = this->mNodes[index]->ContainingElementsBegin();
         !found_target && it != this->mNodes[index]->ContainingElementsEnd();
         ++it)
    {
        Element <ELEMENT_DIM,SPACE_DIM>* p_element = GetElement(*it);
        for (unsigned i=0; i<=ELEMENT_DIM && !found_target; i++)
        {
            target_index = p_element->GetNodeGlobalIndex(i);
            try
            {
                MoveMergeNode(index, target_index, false);
                found_target = true;
            }
            catch (Exception e)
            {
                // Just try the next node
            }
        }
    }
    if (!found_target)
    {
        EXCEPTION("Failure to delete node");
    }

    MoveMergeNode(index, target_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::DeleteNodePriorToReMesh(unsigned index)
{
    this->mNodes[index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::MoveMergeNode(unsigned index,
        unsigned targetIndex,
        bool concreteMove)
{
    if (this->mNodes[index]->IsDeleted())
    {
        EXCEPTION("Trying to move a deleted node");
    }

    if (index == targetIndex)
    {
        EXCEPTION("Trying to merge a node with itself");
    }
    if (this->mNodes[index]->IsBoundaryNode())
    {
        if (!this->mNodes[targetIndex]->IsBoundaryNode())
        {
            EXCEPTION("A boundary node can only be moved on to another boundary node");
        }
    }
    std::set<unsigned> unshared_element_indices;
    std::set_difference(this->mNodes[index]->rGetContainingElementIndices().begin(),
                        this->mNodes[index]->rGetContainingElementIndices().end(),
                        this->mNodes[targetIndex]->rGetContainingElementIndices().begin(),
                        this->mNodes[targetIndex]->rGetContainingElementIndices().end(),
                        std::inserter(unshared_element_indices, unshared_element_indices.begin()));


    if (unshared_element_indices.size() == this->mNodes[index]->rGetContainingElementIndices().size())
    {
        EXCEPTION("These nodes cannot be merged since they are not neighbours");
    }

    std::set<unsigned> unshared_boundary_element_indices;
    std::set_difference(this->mNodes[index]->rGetContainingBoundaryElementIndices().begin(),
                        this->mNodes[index]->rGetContainingBoundaryElementIndices().end(),
                        this->mNodes[targetIndex]->rGetContainingBoundaryElementIndices().begin(),
                        this->mNodes[targetIndex]->rGetContainingBoundaryElementIndices().end(),
                        std::inserter(unshared_boundary_element_indices, unshared_boundary_element_indices.begin()));


    if (this->mNodes[index]->IsBoundaryNode())
    {
        if (unshared_boundary_element_indices.size()
            == this->mNodes[index]->rGetContainingBoundaryElementIndices().size())
        {
            //May be redundant (only thrown in 1D tests)
            EXCEPTION("These nodes cannot be merged since they are not neighbours on the boundary");
        }
    }

    this->mNodes[index]->rGetModifiableLocation() = this->mNodes[targetIndex]->rGetLocation();

    for (std::set<unsigned>::const_iterator element_iter=unshared_element_indices.begin();
             element_iter != unshared_element_indices.end();
             element_iter++)
    {
        try
        {
            if (SPACE_DIM == ELEMENT_DIM)
            {
                this->GetElement(*element_iter)->CalculateInverseJacobian(this->mElementJacobians[(*element_iter)],
                                                                          this->mElementJacobianDeterminants[(*element_iter)],
                                                                          this->mElementInverseJacobians[ (*element_iter) ]);
            }
            else
            {
                this->GetElement(*element_iter)->CalculateWeightedDirection(this->mElementWeightedDirections[(*element_iter)],
                                                                            this->mElementJacobianDeterminants[(*element_iter)]);
            }

            if (concreteMove)
            {
                this->GetElement(*element_iter)->ReplaceNode(this->mNodes[index], this->mNodes[targetIndex]);
            }

        }
        catch (Exception e)
        {
            EXCEPTION("Moving node caused an element to have a non-positive Jacobian determinant");
        }
    }

    for (std::set<unsigned>::const_iterator boundary_element_iter=
                 unshared_boundary_element_indices.begin();
             boundary_element_iter != unshared_boundary_element_indices.end();
             boundary_element_iter++)
    {

        this->GetBoundaryElement(*boundary_element_iter)->CalculateWeightedDirection(this->mBoundaryElementWeightedDirections[(*boundary_element_iter)],
                                                                                     this->mBoundaryElementJacobianDeterminants[(*boundary_element_iter)]);

        if (concreteMove)
        {
            this->GetBoundaryElement(*boundary_element_iter)->ReplaceNode(this->mNodes[index], this->mNodes[targetIndex]);
        }
    }

    std::set<unsigned> shared_element_indices;
    std::set_intersection(this->mNodes[index]->rGetContainingElementIndices().begin(),
                          this->mNodes[index]->rGetContainingElementIndices().end(),
                          this->mNodes[targetIndex]->rGetContainingElementIndices().begin(),
                          this->mNodes[targetIndex]->rGetContainingElementIndices().end(),
                          std::inserter(shared_element_indices, shared_element_indices.begin()));

    for (std::set<unsigned>::const_iterator element_iter=shared_element_indices.begin();
             element_iter != shared_element_indices.end();
             element_iter++)
    {
        if (concreteMove)
        {
            this->GetElement(*element_iter)->MarkAsDeleted();
            this->mElementJacobianDeterminants[ (*element_iter) ] = 0.0; //This used to be done in MarkAsDeleted
            mDeletedElementIndices.push_back(*element_iter);
        }
        else
        {
            this->mElementJacobianDeterminants[ (*element_iter) ] = 0.0;
        }
    }


    std::set<unsigned> shared_boundary_element_indices;
    std::set_intersection(this->mNodes[index]->rGetContainingBoundaryElementIndices().begin(),
                          this->mNodes[index]->rGetContainingBoundaryElementIndices().end(),
                          this->mNodes[targetIndex]->rGetContainingBoundaryElementIndices().begin(),
                          this->mNodes[targetIndex]->rGetContainingBoundaryElementIndices().end(),
                          std::inserter(shared_boundary_element_indices, shared_boundary_element_indices.begin()));

    for (std::set<unsigned>::const_iterator boundary_element_iter=shared_boundary_element_indices.begin();
             boundary_element_iter != shared_boundary_element_indices.end();
             boundary_element_iter++)
    {
        if (concreteMove)
        {
            this->GetBoundaryElement(*boundary_element_iter)->MarkAsDeleted();
            this->mBoundaryElementJacobianDeterminants[ (*boundary_element_iter) ] = 0.0; //This used to be done in MarkAsDeleted
            mDeletedBoundaryElementIndices.push_back(*boundary_element_iter);
        }
        else
        {
            this->mBoundaryElementJacobianDeterminants[ (*boundary_element_iter) ] = 0.0;
            this->mBoundaryElementWeightedDirections[ (*boundary_element_iter) ] = zero_vector<double>(SPACE_DIM);
        }
    }

    if (concreteMove)
    {
        this->mNodes[index]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(index);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMesh<ELEMENT_DIM, SPACE_DIM>::RefineElement(
    Element<ELEMENT_DIM,SPACE_DIM>* pElement,
    ChastePoint<SPACE_DIM> point)
{
    //Check that the point is in the element
    if (pElement->IncludesPoint(point, true) == false)
    {
        EXCEPTION("RefineElement could not be started (point is not in element)");
    }

    // Add a new node from the point that is passed to RefineElement
    unsigned new_node_index = AddNode(new Node<SPACE_DIM>(0, point.rGetLocation()));
    // Note: the first argument is the index of the node, which is going to be
    //       overriden by AddNode, so it can safely be ignored

    // This loop constructs the extra elements which are going to fill the space
    for (unsigned i = 0; i < ELEMENT_DIM; i++)
    {

        // First, make a copy of the current element making sure we update its index
        unsigned new_elt_index;
        if (mDeletedElementIndices.empty())
        {
            new_elt_index = this->mElements.size();
        }
        else
        {
            new_elt_index = mDeletedElementIndices.back();
            mDeletedElementIndices.pop_back();
        }

        Element<ELEMENT_DIM,SPACE_DIM>* p_new_element=
            new Element<ELEMENT_DIM,SPACE_DIM>(*pElement, new_elt_index);

        // Second, update the node in the element with the new one
        p_new_element->UpdateNode(ELEMENT_DIM-1-i, this->mNodes[new_node_index]);

        // Third, add the new element to the set
        if ((unsigned) new_elt_index == this->mElements.size())
        {
            this->mElements.push_back(p_new_element);
        }
        else
        {
            delete this->mElements[new_elt_index];
            this->mElements[new_elt_index] = p_new_element;
        }
    }

    // Lastly, update the last node in the element to be refined
    pElement->UpdateNode(ELEMENT_DIM, this->mNodes[new_node_index]);

    return new_node_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::DeleteBoundaryNodeAt(unsigned index)
{
    if (!this->mNodes[index]->IsBoundaryNode() )
    {
        EXCEPTION(" You may only delete a boundary node ");
    }

    this->mNodes[index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(index);

    // Update the boundary node vector
    typename std::vector<Node<SPACE_DIM>*>::iterator b_node_iter
    = std::find(this->mBoundaryNodes.begin(), this->mBoundaryNodes.end(), this->mNodes[index]);
    this->mBoundaryNodes.erase(b_node_iter);

    // Remove boundary elements containing this node
    std::set<unsigned> boundary_element_indices = this->mNodes[index]->rGetContainingBoundaryElementIndices();
    std::set<unsigned>::const_iterator boundary_element_indices_iterator = boundary_element_indices.begin();
    while (boundary_element_indices_iterator != boundary_element_indices.end())
    {
        BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_boundary_element = this->GetBoundaryElement(*boundary_element_indices_iterator);
        p_boundary_element->MarkAsDeleted();
        mDeletedBoundaryElementIndices.push_back(*boundary_element_indices_iterator);
        boundary_element_indices_iterator++;
    }

    // Remove elements containing this node
    std::set<unsigned> element_indices = this->mNodes[index]->rGetContainingElementIndices();
    std::set<unsigned>::const_iterator element_indices_iterator = element_indices.begin();
    while (element_indices_iterator != element_indices.end())
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*element_indices_iterator);
        for (unsigned i=0; i<p_element->GetNumNodes(); i++)
        {
            Node<SPACE_DIM>* p_node = p_element->GetNode(i);
            if (!p_node->IsDeleted())
            {
                p_node->SetAsBoundaryNode();
                // Update the boundary node vector
                this->mBoundaryNodes.push_back(p_node);
            }
        }
        p_element->MarkAsDeleted();
        mDeletedElementIndices.push_back(p_element->GetIndex());
        element_indices_iterator++;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::ReIndex(NodeMap& map)
{
    assert(!mAddedNodes);
    map.Resize(this->GetNumAllNodes());

    std::vector<Element<ELEMENT_DIM, SPACE_DIM> *> live_elements;

    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        assert(i==this->mElements[i]->GetIndex()); // We need this to be true to be able to reindex the Jacobian cache
        if (this->mElements[i]->IsDeleted())
        {
            delete this->mElements[i];
        }
        else
        {
            live_elements.push_back(this->mElements[i]);

            unsigned this_element_index = this->mElements[i]->GetIndex();
            if (SPACE_DIM == ELEMENT_DIM)
            {
                this->mElementJacobians[live_elements.size()-1] = this->mElementJacobians[this_element_index];
                this->mElementInverseJacobians[live_elements.size()-1] = this->mElementInverseJacobians[this_element_index];
            }
            else
            {
                this->mElementWeightedDirections[live_elements.size()-1] = this->mElementWeightedDirections[this_element_index];
            }
            this->mElementJacobianDeterminants[live_elements.size()-1] = this->mElementJacobianDeterminants[this_element_index];
        }
    }

    assert(mDeletedElementIndices.size() == this->mElements.size()-live_elements.size());
    mDeletedElementIndices.clear();
    this->mElements = live_elements;
    unsigned num_elements = this->mElements.size();

    if (SPACE_DIM == ELEMENT_DIM)
    {
        this->mElementJacobians.resize(num_elements);
        this->mElementInverseJacobians.resize(num_elements);
    }
    else
    {
        this->mElementWeightedDirections.resize(num_elements);
    }
    this->mElementJacobianDeterminants.resize(num_elements);

    std::vector<Node<SPACE_DIM> *> live_nodes;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->mNodes[i]->IsDeleted())
        {
            delete this->mNodes[i];
            map.SetDeleted(i);
        }
        else
        {
            live_nodes.push_back(this->mNodes[i]);
            // the nodes will have their index set to be the index into the live_nodes
            // vector further down
            map.SetNewIndex(i, (unsigned)(live_nodes.size()-1));
        }
    }

    assert(mDeletedNodeIndices.size() == this->mNodes.size()-live_nodes.size());
    this->mNodes = live_nodes;
    mDeletedNodeIndices.clear();

    std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *> live_boundary_elements;
    for (unsigned i=0; i<this->mBoundaryElements.size(); i++)
    {
        if (this->mBoundaryElements[i]->IsDeleted())
        {
            delete this->mBoundaryElements[i];
        }
        else
        {
            live_boundary_elements.push_back(this->mBoundaryElements[i]);

            this->mBoundaryElementWeightedDirections[live_boundary_elements.size()-1] = this->mBoundaryElementWeightedDirections[this->mBoundaryElements[i]->GetIndex()];
            this->mBoundaryElementJacobianDeterminants[live_boundary_elements.size()-1] = this->mBoundaryElementJacobianDeterminants[this->mBoundaryElements[i]->GetIndex()];
        }
    }

    assert(mDeletedBoundaryElementIndices.size() == this->mBoundaryElements.size()-live_boundary_elements.size());
    this->mBoundaryElements = live_boundary_elements;
    mDeletedBoundaryElementIndices.clear();

    unsigned num_boundary_elements = this->mBoundaryElements.size();

    this->mBoundaryElementWeightedDirections.resize(num_boundary_elements);
    this->mBoundaryElementJacobianDeterminants.resize(num_boundary_elements);

    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->SetIndex(i);
    }
    for (unsigned i=0; i<this->mElements.size(); i++)
    {

        this->mElements[i]->ResetIndex(i);
    }
    for (unsigned i=0; i<this->mBoundaryElements.size(); i++)
    {
        this->mBoundaryElements[i]->ResetIndex(i);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(NodeMap& map)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE

    // Avoid some triangle/tetgen errors: need at least four
    // nodes for tetgen, and at least three for triangle
    if (GetNumNodes() <= SPACE_DIM)
    {
        EXCEPTION("The number of nodes must exceed the spatial dimension.");
    }

    // Make sure the map is big enough
    map.Resize(this->GetNumAllNodes());
    if (mAddedNodes || !mDeletedNodeIndices.empty())
    {
        // Size of mesh is about to change
        if (this->mpDistributedVectorFactory)
        {
            delete this->mpDistributedVectorFactory;
            this->mpDistributedVectorFactory = new DistributedVectorFactory(this->GetNumNodes());
        }
    }
    if (SPACE_DIM==1)
    {
        // Store the node locations
        std::vector<c_vector<double, SPACE_DIM> > old_node_locations;
        unsigned new_index = 0;
        for (unsigned i=0; i<this->GetNumAllNodes(); i++)
        {
            if (this->mNodes[i]->IsDeleted())
            {
                map.SetDeleted(i);
            }
            else
            {
                map.SetNewIndex(i, new_index);
                old_node_locations.push_back(this->mNodes[i]->rGetLocation());
                new_index++;
            }
        }

        // Remove current data
        Clear();

        // Construct the nodes and boundary nodes
        for (unsigned node_index=0; node_index<old_node_locations.size(); node_index++)
        {
            Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, old_node_locations[node_index], false);
            this->mNodes.push_back(p_node);

            // As we're in 1D, the boundary nodes are simply at either end of the mesh
            if ( node_index==0 || node_index==old_node_locations.size()-1 )
            {
                this->mBoundaryNodes.push_back(p_node);
            }
        }

        // Create a map between node indices and node locations
        std::map<double, unsigned> location_index_map;
        for (unsigned i=0; i<this->mNodes.size(); i++)
        {
            location_index_map[this->mNodes[i]->rGetLocation()[0]] = this->mNodes[i]->GetIndex();
        }

        // Use this map to generate a vector of node indices that are ordered spatially
        std::vector<unsigned> node_indices_ordered_spatially;
        for (std::map<double, unsigned>::iterator iter = location_index_map.begin();
             iter != location_index_map.end();
             ++iter)
        {
            node_indices_ordered_spatially.push_back(iter->second);
        }

        // Construct the elements
        this->mElements.reserve(old_node_locations.size()-1);
        for (unsigned element_index=0; element_index<old_node_locations.size()-1; element_index++)
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            for (unsigned j=0; j<2; j++)
            {
                unsigned global_node_index = node_indices_ordered_spatially[element_index + j];
                assert(global_node_index < this->mNodes.size());
                nodes.push_back(this->mNodes[global_node_index]);
            }
            this->mElements.push_back(new Element<ELEMENT_DIM, SPACE_DIM>(element_index, nodes));
        }

        // Construct the two boundary elements - as we're in 1D, these are simply at either end of the mesh
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[0]);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>(0, nodes));

        nodes.clear();
        nodes.push_back(this->mNodes[old_node_locations.size()-1]);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>(1, nodes));

        this->RefreshJacobianCachedData();
    }
    else if (SPACE_DIM==2)  // In 2D, remesh using triangle via library calls
    {
        struct triangulateio mesher_input, mesher_output;
        this->InitialiseTriangulateIo(mesher_input);
        this->InitialiseTriangulateIo(mesher_output);

        this->ExportToMesher(map, mesher_input);

        // Library call
        triangulate((char*)"Qze", &mesher_input, &mesher_output, NULL);

        try
        {
            this->ImportFromMesher(mesher_output, mesher_output.numberoftriangles, mesher_output.trianglelist, mesher_output.numberofedges, mesher_output.edgelist, mesher_output.edgemarkerlist);
        }
        catch (Exception& e)
        {
            //A zero-area element was detected so we couldn't carry on
            //Tidy up triangle
            this->FreeTriangulateIo(mesher_input);
            this->FreeTriangulateIo(mesher_output);
            throw e;
        }
        //Tidy up triangle
        this->FreeTriangulateIo(mesher_input);
        this->FreeTriangulateIo(mesher_output);
    }
    else // in 3D, remesh using tetgen
    {

        struct tetgen::tetgenio mesher_input, mesher_output;

        this->ExportToMesher(map, mesher_input);

        // Library call
        tetgen::tetrahedralize((char*)"Qz", &mesher_input, &mesher_output);

        this->ImportFromMesher(mesher_output, mesher_output.numberoftetrahedra, mesher_output.tetrahedronlist, mesher_output.numberoftrifaces, mesher_output.trifacelist, NULL);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh()
{
    NodeMap map(GetNumNodes());
    ReMesh(map);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableMesh<ELEMENT_DIM, SPACE_DIM>::CheckIsVoronoi(Element<ELEMENT_DIM, SPACE_DIM>* pElement, double maxPenetration)
{
    assert(ELEMENT_DIM == SPACE_DIM);
    unsigned num_nodes = pElement->GetNumNodes();
    std::set<unsigned> neighbouring_elements_indices;
    std::set< Element<ELEMENT_DIM,SPACE_DIM> *> neighbouring_elements;
    std::set<unsigned> neighbouring_nodes_indices;

    // Form a set of neighbouring elements via the nodes
    for (unsigned i=0; i<num_nodes; i++)
    {
        Node<SPACE_DIM>* p_node = pElement->GetNode(i);
        neighbouring_elements_indices = p_node->rGetContainingElementIndices();
        for (std::set<unsigned>::const_iterator it = neighbouring_elements_indices.begin();
             it != neighbouring_elements_indices.end();
             ++it)
        {
            neighbouring_elements.insert(this->GetElement(*it));
        }
    }
    neighbouring_elements.erase(pElement);

    // For each neighbouring element find the supporting nodes
    typedef typename std::set<Element<ELEMENT_DIM,SPACE_DIM> *>::const_iterator ElementIterator;

    for (ElementIterator it = neighbouring_elements.begin();
         it != neighbouring_elements.end();
         ++it)
    {
        for (unsigned i=0; i<num_nodes; i++)
        {
            neighbouring_nodes_indices.insert((*it)->GetNodeGlobalIndex(i));
        }
    }

    // Remove the nodes that support this element
    for (unsigned i = 0; i < num_nodes; i++)
    {
        neighbouring_nodes_indices.erase(pElement->GetNodeGlobalIndex(i));
    }

    // Get the circumsphere information
    c_vector<double, SPACE_DIM+1> this_circum_centre;

    this_circum_centre = pElement->CalculateCircumsphere(this->mElementJacobians[pElement->GetIndex()], this->mElementInverseJacobians[pElement->GetIndex()]);

    // Copy the actualy circumcentre into a smaller vector
    c_vector<double, ELEMENT_DIM> circum_centre;
    for (unsigned i=0; i<ELEMENT_DIM; i++)
    {
        circum_centre[i] = this_circum_centre[i];
    }

    for (std::set<unsigned>::const_iterator it = neighbouring_nodes_indices.begin();
         it != neighbouring_nodes_indices.end();
         ++it)
    {
        c_vector<double, ELEMENT_DIM> node_location = this->GetNode(*it)->rGetLocation();

        // Calculate vector from circumcenter to node
        node_location -= circum_centre;

        // This is to calculate the squared distance betweeen them
        double squared_distance = inner_prod(node_location, node_location);

        // If the squared idstance is less than the elements circum-radius(squared),
        // then the Voronoi property is violated.
        if (squared_distance < this_circum_centre[ELEMENT_DIM])
        {
            // We know the node is inside the circumsphere, but we don't know how far
            double radius = sqrt(this_circum_centre[ELEMENT_DIM]);
            double distance = radius - sqrt(squared_distance);

            // If the node penetration is greater than supplied maximum penetration factor
            if (distance/radius > maxPenetration)
            {
                return false;
            }
        }
    }
    return true;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableMesh<ELEMENT_DIM, SPACE_DIM>::CheckIsVoronoi(double maxPenetration)
{
    // Looping through all the elements in the mesh
    /// \todo use ElementIterator here?
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        // Check if the element is not deleted
        if (!this->mElements[i]->IsDeleted())
        {
            // Checking the Voronoi of the Element
            if (CheckIsVoronoi(this->mElements[i], maxPenetration) == false)
            {
                return false;
            }
        }
    }
    return true;
}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class MutableMesh<1,1>;
template class MutableMesh<1,2>;
template class MutableMesh<1,3>;
template class MutableMesh<2,2>;
template class MutableMesh<2,3>;
template class MutableMesh<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MutableMesh)
