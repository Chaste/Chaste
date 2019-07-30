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

#include "TetrahedralMesh.hpp"

#include <cassert>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>

#include "BoundaryElement.hpp"
#include "Element.hpp"
#include "Exception.hpp"
#include "Node.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

// Jonathan Shewchuk's triangle and Hang Si's tetgen
#define REAL double
#define VOID void
#include "tetgen.h"
#include "triangle.h"
#undef REAL
#undef VOID

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::TetrahedralMesh()
{
    Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)
{
    assert(rMeshReader.HasNodePermutation() == false);
    this->mMeshFileBaseName = rMeshReader.GetMeshFileBaseName();

    // Record number of corner nodes
    unsigned num_nodes = rMeshReader.GetNumNodes();

    /*
     * Reserve memory for nodes, so we don't have problems with
     * pointers stored in elements becoming invalid.
     */
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    //typename std::map<std::pair<unsigned,unsigned>,unsigned>::const_iterator iterator;
    //std::map<std::pair<unsigned,unsigned>,unsigned> internal_nodes_map;

    // Add nodes
    std::vector<double> coords;
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        coords = rMeshReader.GetNextNode();
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, coords, false);

        for (unsigned i = 0; i < rMeshReader.GetNodeAttributes().size(); i++)
        {
            double attribute = rMeshReader.GetNodeAttributes()[i];
            p_node->AddNodeAttribute(attribute);
        }
        this->mNodes.push_back(p_node);
    }

    //unsigned new_node_index = mNumCornerNodes;

    rMeshReader.Reset();
    // Add elements
    //new_node_index = mNumCornerNodes;
    this->mElements.reserve(rMeshReader.GetNumElements());

    for (unsigned element_index = 0; element_index < (unsigned)rMeshReader.GetNumElements(); element_index++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();
        std::vector<Node<SPACE_DIM>*> nodes;

        /*
         * NOTE: currently just reading element vertices from mesh reader - even if it
         * does contain information about internal nodes (ie for quadratics) this is
         * ignored here and used elsewhere: ie don't do this:
         *   unsigned nodes_size = node_indices.size();
         */
        for (unsigned j = 0; j < ELEMENT_DIM + 1; j++) // num vertices=ELEMENT_DIM+1, may not be equal to nodes_size.
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        Element<ELEMENT_DIM, SPACE_DIM>* p_element = new Element<ELEMENT_DIM, SPACE_DIM>(element_index, nodes);

        this->mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            double attribute_value = element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }

    // Add boundary elements and nodes
    for (unsigned face_index = 0; face_index < (unsigned)rMeshReader.GetNumFaces(); face_index++)
    {
        ElementData face_data = rMeshReader.GetNextFaceData();
        std::vector<unsigned> node_indices = face_data.NodeIndices;

        /*
         * NOTE: unlike the above where we just read element *vertices* from mesh reader, here we are
         * going to read a quadratic mesh with internal elements.
         * (There are only a few meshes with internals in the face file that we might as well use them.)
         *
         */
        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned node_index = 0; node_index < node_indices.size(); node_index++)
        {
            assert(node_indices[node_index] < this->mNodes.size());
            // Add Node pointer to list for creating an element
            nodes.push_back(this->mNodes[node_indices[node_index]]);
        }

        // This is a boundary face, so ensure all its nodes are marked as boundary nodes
        for (unsigned j = 0; j < nodes.size(); j++)
        {
            if (!nodes[j]->IsBoundaryNode())
            {
                nodes[j]->SetAsBoundaryNode();
                this->mBoundaryNodes.push_back(nodes[j]);
            }

            // Register the index that this bounday element will have with the node
            nodes[j]->AddBoundaryElement(face_index);
        }

        // The added elements will be deleted in our destructor
        BoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* p_boundary_element = new BoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>(face_index, nodes);
        this->mBoundaryElements.push_back(p_boundary_element);

        if (rMeshReader.GetNumFaceAttributes() > 0)
        {
            assert(rMeshReader.GetNumFaceAttributes() == 1);
            double attribute_value = face_data.AttributeValue;
            p_boundary_element->SetAttribute(attribute_value);
        }
    }

    RefreshJacobianCachedData();

    rMeshReader.Reset();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReadNodesPerProcessorFile(const std::string& rNodesPerProcessorFile)
{
    std::vector<unsigned> nodes_per_processor_vec;

    std::ifstream file_stream(rNodesPerProcessorFile.c_str());
    if (file_stream.is_open())
    {
        while (file_stream)
        {
            unsigned nodes_per_processor;
            file_stream >> nodes_per_processor;

            if (file_stream)
            {
                nodes_per_processor_vec.push_back(nodes_per_processor);
            }
        }
    }
    else
    {
        EXCEPTION("Unable to read nodes per processor file " + rNodesPerProcessorFile);
    }

    unsigned sum = 0;
    for (unsigned i = 0; i < nodes_per_processor_vec.size(); i++)
    {
        sum += nodes_per_processor_vec[i];
    }

    if (sum != this->GetNumNodes())
    {
        EXCEPTION("Sum of nodes per processor, " << sum
                                                 << ", not equal to number of nodes in mesh, " << this->GetNumNodes());
    }

    unsigned num_owned = nodes_per_processor_vec[PetscTools::GetMyRank()];

    if (nodes_per_processor_vec.size() != PetscTools::GetNumProcs())
    {
        EXCEPTION("Number of processes doesn't match the size of the nodes-per-processor file");
    }
    delete this->mpDistributedVectorFactory;
    this->mpDistributedVectorFactory = new DistributedVectorFactory(this->GetNumNodes(), num_owned);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CheckIsConforming()
{
    /*
     * Each face of each element is a set of node indices.
     * We form a set of these in order to get their parity:
     *   all faces which appear once are inserted into the set;
     *   all faces which appear twice are inserted and then removed from the set;
     *   we're assuming that faces never appear more than twice.
     */
    std::set<std::set<unsigned> > odd_parity_faces;

    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        for (unsigned face_index = 0; face_index <= ELEMENT_DIM; face_index++)
        {
            std::set<unsigned> face_info;
            for (unsigned node_index = 0; node_index <= ELEMENT_DIM; node_index++)
            {
                // Leave one index out each time
                if (node_index != face_index)
                {
                    face_info.insert(iter->GetNodeGlobalIndex(node_index));
                }
            }
            // Face is now formed - attempt to find it
            std::set<std::set<unsigned> >::iterator find_face = odd_parity_faces.find(face_info);
            if (find_face != odd_parity_faces.end())
            {
                // Face was in set, so it now has even parity.
                // Remove it via the iterator
                odd_parity_faces.erase(find_face);
            }
            else
            {
                // Face is not in set so it now has odd parity. Insert it
                odd_parity_faces.insert(face_info);
            }
        }
    }

    /*
     * At this point the odd parity faces should be the same as the
     * boundary elements. We could check this explicitly or we could
     * just count them.
     */
    return (odd_parity_faces.size() == this->GetNumBoundaryElements());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetVolume()
{
    double mesh_volume = 0.0;

    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        mesh_volume += iter->GetVolume(mElementJacobianDeterminants[iter->GetIndex()]);
    }

    return mesh_volume;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetSurfaceArea()
{
    // ELEMENT_DIM-1 is the dimension of the boundary element
    assert(ELEMENT_DIM >= 1);
    const unsigned bound_element_dim = ELEMENT_DIM - 1;
    assert(bound_element_dim < 3);
    if (bound_element_dim == 0)
    {
        return 0.0;
    }

    double mesh_surface = 0.0;
    typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator it = this->GetBoundaryElementIteratorBegin();

    while (it != this->GetBoundaryElementIteratorEnd())
    {
        mesh_surface += mBoundaryElementJacobianDeterminants[(*it)->GetIndex()];
        it++;
    }

    if (bound_element_dim == 2)
    {
        mesh_surface /= 2.0;
    }

    return mesh_surface;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
{
    // Make a permutation vector of the identity
    RandomNumberGenerator* p_rng = RandomNumberGenerator::Instance();
    std::vector<unsigned> perm;
    p_rng->Shuffle(this->mNodes.size(), perm);

    // Call the non-random version
    PermuteNodes(perm);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes(const std::vector<unsigned>& perm)
{
    // Let's not do this if there are any deleted nodes
    assert(this->GetNumAllNodes() == this->GetNumNodes());

    assert(perm.size() == this->mNodes.size());

    // Copy the node pointers
    std::vector<Node<SPACE_DIM>*> copy_m_nodes;
    copy_m_nodes.assign(this->mNodes.begin(), this->mNodes.end());

    for (unsigned original_index = 0; original_index < this->mNodes.size(); original_index++)
    {
        assert(perm[original_index] < this->mNodes.size());
        //perm[original_index] holds the new assigned index of that node
        this->mNodes[perm[original_index]] = copy_m_nodes[original_index];
    }

    // Update indices
    for (unsigned index = 0; index < this->mNodes.size(); index++)
    {
        this->mNodes[index]->SetIndex(index);
    }

    // Copy the permutation vector into the mesh
    this->mNodePermutation = perm;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndexWithInitialGuess(const ChastePoint<SPACE_DIM>& rTestPoint, unsigned startingElementGuess, bool strict)
{
    assert(startingElementGuess < this->GetNumElements());

    /*
     * Let m=startingElementGuess, N=num_elem-1.
     * We search from in this order: m, m+1, m+2, .. , N, 0, 1, .., m-1.
     */
    unsigned i = startingElementGuess;
    bool reached_end = false;

    while (!reached_end)
    {
        if (this->mElements[i]->IncludesPoint(rTestPoint, strict))
        {
            assert(!this->mElements[i]->IsDeleted());
            return i;
        }

        // Increment
        i++;
        if (i == this->GetNumElements())
        {
            i = 0;
        }

        // Back to the beginning yet?
        if (i == startingElementGuess)
        {
            reached_end = true;
        }
    }

    // If it's in none of the elements, then throw
    std::stringstream ss;
    ss << "Point [";
    for (unsigned j = 0; (int)j < (int)SPACE_DIM - 1; j++)
    {
        ss << rTestPoint[j] << ",";
    }
    ss << rTestPoint[SPACE_DIM - 1] << "] is not in mesh - all elements tested";
    EXCEPTION(ss.str());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNearestElementIndex(const ChastePoint<SPACE_DIM>& rTestPoint)
{
    EXCEPT_IF_NOT(ELEMENT_DIM == SPACE_DIM); // LCOV_EXCL_LINE // CalculateInterpolationWeights hits an assertion otherwise
    double max_min_weight = -std::numeric_limits<double>::infinity();
    unsigned closest_index = 0;
    for (unsigned i = 0; i < this->mElements.size(); i++)
    {
        c_vector<double, ELEMENT_DIM + 1> weight = this->mElements[i]->CalculateInterpolationWeights(rTestPoint);
        double neg_weight_sum = 0.0;
        for (unsigned j = 0; j <= ELEMENT_DIM; j++)
        {
            if (weight[j] < 0.0)
            {
                neg_weight_sum += weight[j];
            }
        }
        if (neg_weight_sum > max_min_weight)
        {
            max_min_weight = neg_weight_sum;
            closest_index = i;
        }
    }
    assert(!this->mElements[closest_index]->IsDeleted());
    return closest_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndices(const ChastePoint<SPACE_DIM>& rTestPoint)
{
    std::vector<unsigned> element_indices;
    for (unsigned i = 0; i < this->mElements.size(); i++)
    {
        if (this->mElements[i]->IncludesPoint(rTestPoint))
        {
            assert(!this->mElements[i]->IsDeleted());
            element_indices.push_back(i);
        }
    }
    return element_indices;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    // Three loops, just like the destructor. note we don't delete boundary nodes.
    for (unsigned i = 0; i < this->mBoundaryElements.size(); i++)
    {
        delete this->mBoundaryElements[i];
    }
    for (unsigned i = 0; i < this->mElements.size(); i++)
    {
        delete this->mElements[i];
    }
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }

    this->mNodes.clear();
    this->mElements.clear();
    this->mBoundaryElements.clear();
    this->mBoundaryNodes.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetAngleBetweenNodes(unsigned indexA, unsigned indexB)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE
    assert(SPACE_DIM == ELEMENT_DIM); // LCOV_EXCL_LINE

    double x_difference = this->mNodes[indexB]->rGetLocation()[0] - this->mNodes[indexA]->rGetLocation()[0];
    double y_difference = this->mNodes[indexB]->rGetLocation()[1] - this->mNodes[indexA]->rGetLocation()[1];

    if (x_difference == 0)
    {
        if (y_difference > 0)
        {
            return M_PI / 2.0;
        }
        else if (y_difference < 0)
        {
            return -M_PI / 2.0;
        }
        else
        {
            EXCEPTION("Tried to compute polar angle of (0,0)");
        }
    }

    double angle = atan2(y_difference, x_difference);
    return angle;
}

//////////////////////////////////////////////////////////////////////////////
//                          Edge iterator class                             //
//////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::GetNodeA()
{
    assert((*this) != mrMesh.EdgesEnd());
    Element<ELEMENT_DIM, SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
    return p_element->GetNode(mNodeALocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::GetNodeB()
{
    assert((*this) != mrMesh.EdgesEnd());
    Element<ELEMENT_DIM, SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
    return p_element->GetNode(mNodeBLocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::operator!=(const typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& rOther)
{
    return (mElemIndex != rOther.mElemIndex || mNodeALocalIndex != rOther.mNodeALocalIndex || mNodeBLocalIndex != rOther.mNodeBLocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::operator++()
{
    bool already_seen_this_edge;

    unsigned num_elements = mrMesh.GetNumAllElements();
    std::pair<unsigned, unsigned> current_node_pair;
    do
    {
        /*
         * Advance to the next edge in the mesh.
         * Node indices are incremented modulo #nodes_per_elem.
         */
        mNodeBLocalIndex = (mNodeBLocalIndex + 1) % (ELEMENT_DIM + 1);
        if (mNodeBLocalIndex == mNodeALocalIndex)
        {
            mNodeALocalIndex = (mNodeALocalIndex + 1) % (ELEMENT_DIM + 1);
            mNodeBLocalIndex = (mNodeALocalIndex + 1) % (ELEMENT_DIM + 1);
        }

        if (mNodeALocalIndex == 0 && mNodeBLocalIndex == 1) // advance to next element...
        {
            // ...skipping deleted ones
            do
            {
                mElemIndex++;
            } while (mElemIndex != num_elements && mrMesh.GetElement(mElemIndex)->IsDeleted());
        }

        if (mElemIndex != num_elements)
        {
            Element<ELEMENT_DIM, SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
            unsigned node_a_global_index = p_element->GetNodeGlobalIndex(mNodeALocalIndex);
            unsigned node_b_global_index = p_element->GetNodeGlobalIndex(mNodeBLocalIndex);
            if (node_b_global_index < node_a_global_index)
            {
                // Swap them over
                unsigned temp = node_a_global_index;
                node_a_global_index = node_b_global_index;
                node_b_global_index = temp;
            }

            // Check we haven't seen it before
            current_node_pair = std::pair<unsigned, unsigned>(node_a_global_index, node_b_global_index);
            already_seen_this_edge = (mEdgesVisited.count(current_node_pair) != 0);
        }
        else
        {
            already_seen_this_edge = false;
        }
    }

    while (already_seen_this_edge);
    mEdgesVisited.insert(current_node_pair);

    return (*this);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::EdgeIterator(TetrahedralMesh& rMesh, unsigned elemIndex)
        : mrMesh(rMesh),
          mElemIndex(elemIndex),
          mNodeALocalIndex(0),
          mNodeBLocalIndex(1)
{
    if (elemIndex == mrMesh.GetNumAllElements())
    {
        return;
    }

    mEdgesVisited.clear();

    // Add the current node pair to the store
    unsigned node_a_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeALocalIndex);
    unsigned node_b_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeBLocalIndex);
    if (node_b_global_index < node_a_global_index)
    {
        // Swap them over
        unsigned temp = node_a_global_index;
        node_a_global_index = node_b_global_index;
        node_b_global_index = temp;
    }

    // Check we haven't seen it before
    std::pair<unsigned, unsigned> current_node_pair = std::pair<unsigned, unsigned>(node_a_global_index, node_b_global_index);
    mEdgesVisited.insert(current_node_pair);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgesBegin()
{
    unsigned first_element_index = 0;
    while (first_element_index != this->GetNumAllElements() && this->GetElement(first_element_index)->IsDeleted())
    {
        first_element_index++;
    }
    return EdgeIterator(*this, first_element_index);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgesEnd()
{
    return EdgeIterator(*this, this->GetNumAllElements());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshMesh()
{
    RefreshJacobianCachedData();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    assert(index < this->mBoundaryElements.size());
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshJacobianCachedData()
{
    unsigned num_elements = this->GetNumAllElements();
    unsigned num_boundary_elements = this->GetNumAllBoundaryElements();

    // Make sure we have enough space
    this->mElementJacobians.resize(num_elements);
    this->mElementInverseJacobians.resize(num_elements);

    if (ELEMENT_DIM < SPACE_DIM)
    {
        this->mElementWeightedDirections.resize(num_elements);
    }

    this->mBoundaryElementWeightedDirections.resize(num_boundary_elements);

    this->mElementJacobianDeterminants.resize(num_elements);
    this->mBoundaryElementJacobianDeterminants.resize(num_boundary_elements);

    // Update caches
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        unsigned index = iter->GetIndex();
        iter->CalculateInverseJacobian(this->mElementJacobians[index], this->mElementJacobianDeterminants[index], this->mElementInverseJacobians[index]);
    }

    if (ELEMENT_DIM < SPACE_DIM)
    {
        for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
             iter != this->GetElementIteratorEnd();
             ++iter)
        {
            unsigned index = iter->GetIndex();
            iter->CalculateWeightedDirection(this->mElementWeightedDirections[index], this->mElementJacobianDeterminants[index]);
        }
    }

    for (typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator itb = this->GetBoundaryElementIteratorBegin();
         itb != this->GetBoundaryElementIteratorEnd();
         itb++)
    {
        unsigned index = (*itb)->GetIndex();
        (*itb)->CalculateWeightedDirection(this->mBoundaryElementWeightedDirections[index], this->mBoundaryElementJacobianDeterminants[index]);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double& rJacobianDeterminant) const
{
    assert(ELEMENT_DIM <= SPACE_DIM);
    assert(elementIndex < this->mElementJacobians.size());
    rJacobian = this->mElementJacobians[elementIndex];
    rJacobianDeterminant = this->mElementJacobianDeterminants[elementIndex];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double& rJacobianDeterminant, c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const
{
    assert(ELEMENT_DIM <= SPACE_DIM); // LCOV_EXCL_LINE
    assert(elementIndex < this->mElementInverseJacobians.size());
    rInverseJacobian = this->mElementInverseJacobians[elementIndex];
    rJacobian = this->mElementJacobians[elementIndex];
    rJacobianDeterminant = this->mElementJacobianDeterminants[elementIndex];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const
{
    assert(ELEMENT_DIM < SPACE_DIM); // LCOV_EXCL_LINE
    assert(elementIndex < this->mElementWeightedDirections.size());
    rWeightedDirection = this->mElementWeightedDirections[elementIndex];
    rJacobianDeterminant = this->mElementJacobianDeterminants[elementIndex];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForBoundaryElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const
{
    assert(elementIndex < this->mBoundaryElementWeightedDirections.size());
    rWeightedDirection = this->mBoundaryElementWeightedDirections[elementIndex];
    rJacobianDeterminant = this->mBoundaryElementJacobianDeterminants[elementIndex];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::InitialiseTriangulateIo(triangulateio& mesherIo)
{
    mesherIo.numberofpoints = 0;
    mesherIo.pointlist = nullptr;
    mesherIo.numberofpointattributes = 0;
    mesherIo.pointattributelist = (double*)nullptr;
    mesherIo.pointmarkerlist = (int*)nullptr;
    mesherIo.numberofsegments = 0;
    mesherIo.numberofholes = 0;
    mesherIo.numberofregions = 0;
    mesherIo.trianglelist = (int*)nullptr;
    mesherIo.triangleattributelist = (double*)nullptr;
    mesherIo.numberoftriangleattributes = 0;
    mesherIo.edgelist = (int*)nullptr;
    mesherIo.edgemarkerlist = (int*)nullptr;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::FreeTriangulateIo(triangulateio& mesherIo)
{
    if (mesherIo.numberofpoints != 0)
    {
        mesherIo.numberofpoints = 0;
        free(mesherIo.pointlist);
    }

    // These (and the above) should actually be safe since we explicity set to NULL above
    free(mesherIo.pointattributelist);
    free(mesherIo.pointmarkerlist);
    free(mesherIo.trianglelist);
    free(mesherIo.triangleattributelist);
    free(mesherIo.edgelist);
    free(mesherIo.edgemarkerlist);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
template <class MESHER_IO>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ExportToMesher(NodeMap& map, MESHER_IO& mesherInput, int* elementList)
{
    if (SPACE_DIM == 2)
    {
        mesherInput.pointlist = (double*)malloc(this->GetNumNodes() * SPACE_DIM * sizeof(double));
    }
    else
    {
        mesherInput.pointlist = new double[this->GetNumNodes() * SPACE_DIM];
    }

    mesherInput.numberofpoints = this->GetNumNodes();
    unsigned new_index = 0;
    for (unsigned i = 0; i < this->GetNumAllNodes(); i++)
    {
        if (this->mNodes[i]->IsDeleted())
        {
            map.SetDeleted(i);
        }
        else
        {
            map.SetNewIndex(i, new_index);
            for (unsigned j = 0; j < SPACE_DIM; j++)
            {
                mesherInput.pointlist[SPACE_DIM * new_index + j] = this->mNodes[i]->rGetLocation()[j];
            }
            new_index++;
        }
    }
    if (elementList != nullptr)
    {
        unsigned element_index = 0;

        // Assume there is enough space for this
        mesherInput.numberofcorners = ELEMENT_DIM + 1;
        for (typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->GetElementIteratorBegin();
             elem_iter != this->GetElementIteratorEnd();
             ++elem_iter)
        {

            for (unsigned j = 0; j <= ELEMENT_DIM; j++)
            {
                elementList[element_index * (ELEMENT_DIM + 1) + j] = (*elem_iter).GetNodeGlobalIndex(j);
            }
            element_index++;
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
template <class MESHER_IO>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ImportFromMesher(MESHER_IO& mesherOutput, unsigned numberOfElements, int* elementList, unsigned numberOfFaces, int* faceList, int* edgeMarkerList)
{
    unsigned nodes_per_element = mesherOutput.numberofcorners;

    assert(nodes_per_element == ELEMENT_DIM + 1 || nodes_per_element == (ELEMENT_DIM + 1) * (ELEMENT_DIM + 2) / 2);

    Clear();

    // Construct the nodes
    for (unsigned node_index = 0; node_index < (unsigned)mesherOutput.numberofpoints; node_index++)
    {
        this->mNodes.push_back(new Node<SPACE_DIM>(node_index, &mesherOutput.pointlist[node_index * SPACE_DIM], false));
    }

    // Construct the elements
    this->mElements.reserve(numberOfElements);

    unsigned real_element_index = 0;
    for (unsigned element_index = 0; element_index < numberOfElements; element_index++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned j = 0; j < ELEMENT_DIM + 1; j++)
        {
            unsigned global_node_index = elementList[element_index * (nodes_per_element) + j];
            assert(global_node_index < this->mNodes.size());
            nodes.push_back(this->mNodes[global_node_index]);
        }

        /*
         * For some reason, tetgen in library mode makes its initial Delaunay mesh
         * with very thin slivers. Hence we expect to ignore some of the elements!
         */
        Element<ELEMENT_DIM, SPACE_DIM>* p_element;
        try
        {
            p_element = new Element<ELEMENT_DIM, SPACE_DIM>(real_element_index, nodes);

            // Shouldn't throw after this point
            this->mElements.push_back(p_element);

            // Add the internals to quadratics
            for (unsigned j = ELEMENT_DIM + 1; j < nodes_per_element; j++)
            {
                unsigned global_node_index = elementList[element_index * nodes_per_element + j];
                assert(global_node_index < this->mNodes.size());
                this->mElements[real_element_index]->AddNode(this->mNodes[global_node_index]);
                this->mNodes[global_node_index]->AddElement(real_element_index);
                this->mNodes[global_node_index]->MarkAsInternal();
            }
            real_element_index++;
        }
        catch (Exception&)
        {
            if (SPACE_DIM == 2)
            {
                WARNING("Triangle has produced a zero area (collinear) element");
            }
            else
            {
                WARNING("Tetgen has produced a zero volume (coplanar) element");
            }
        }
    }

    // Construct the BoundaryElements (and mark boundary nodes)
    unsigned next_boundary_element_index = 0;
    for (unsigned boundary_element_index = 0; boundary_element_index < numberOfFaces; boundary_element_index++)
    {
        /*
         * Tetgen produces only boundary faces (set edgeMarkerList to NULL).
         * Triangle marks which edges are on the boundary.
         */
        if (edgeMarkerList == nullptr || edgeMarkerList[boundary_element_index] == 1)
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            for (unsigned j = 0; j < ELEMENT_DIM; j++)
            {
                unsigned global_node_index = faceList[boundary_element_index * ELEMENT_DIM + j];
                assert(global_node_index < this->mNodes.size());
                nodes.push_back(this->mNodes[global_node_index]);
                if (!nodes[j]->IsBoundaryNode())
                {
                    nodes[j]->SetAsBoundaryNode();
                    this->mBoundaryNodes.push_back(nodes[j]);
                }
            }

            /*
             * For some reason, tetgen in library mode makes its initial Delaunay mesh
             * with very thin slivers. Hence we expect to ignore some of the elements!
             */
            BoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>* p_b_element;
            try
            {
                p_b_element = new BoundaryElement<ELEMENT_DIM - 1, SPACE_DIM>(next_boundary_element_index, nodes);
                this->mBoundaryElements.push_back(p_b_element);
                next_boundary_element_index++;
            }
            // LCOV_EXCL_START
            catch (Exception&)
            {
                // Tetgen is feeding us lies
                /*
                 *  Note: this code is covered in profiling (Test3dOffLatticeRepresentativeSimulation).
                 *  It's hard to replicate Tetgen's behaviour with a unit test.
                 */
                assert(SPACE_DIM == 3);
            }
            // LCOV_EXCL_STOP
        }
    }

    this->RefreshJacobianCachedData();
}

// Explicit instantiation
template class TetrahedralMesh<1, 1>;
template class TetrahedralMesh<1, 2>;
template class TetrahedralMesh<1, 3>;
template class TetrahedralMesh<2, 2>;
template class TetrahedralMesh<2, 3>;
template class TetrahedralMesh<3, 3>;

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template void TetrahedralMesh<2, 2>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
template void TetrahedralMesh<2, 2>::ImportFromMesher<triangulateio>(triangulateio&, unsigned, int*, unsigned, int*, int*);

template void TetrahedralMesh<3, 3>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
template void TetrahedralMesh<3, 3>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned, int*, unsigned, int*, int*);

//The following don't ever need to be instantiated, but are needed to keep some compilers happy
template void TetrahedralMesh<1, 2>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
template void TetrahedralMesh<1, 2>::ImportFromMesher<triangulateio>(triangulateio&, unsigned, int*, unsigned, int*, int*);

template void TetrahedralMesh<1, 3>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
template void TetrahedralMesh<1, 3>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned, int*, unsigned, int*, int*);
template void TetrahedralMesh<2, 3>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
template void TetrahedralMesh<2, 3>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned, int*, unsigned, int*, int*);

// Intel compilation with IPO thinks that it's missing some bizarre instantiations
template void TetrahedralMesh<3u, 3u>::ImportFromMesher<triangulateio>(triangulateio&, unsigned int, int*, unsigned int, int*, int*);
template void TetrahedralMesh<1u, 1u>::ImportFromMesher<triangulateio>(triangulateio&, unsigned int, int*, unsigned int, int*, int*);
template void TetrahedralMesh<1u, 1u>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned int, int*, unsigned int, int*, int*);
template void TetrahedralMesh<2u, 2u>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned int, int*, unsigned int, int*, int*);
template void TetrahedralMesh<1u, 1u>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);

// Compiler on ARC cluster HAL requires the following
template void TetrahedralMesh<3u, 3u>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
template void TetrahedralMesh<1u, 1u>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
template void TetrahedralMesh<2u, 2u>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);

// Intel v11 compilation thinks that it's missing even more bizarre instantiations
//template void TetrahedralMesh<2,2>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
//template void TetrahedralMesh<3,3>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
//template void TetrahedralMesh<1,3>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
//template void TetrahedralMesh<1,1>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
//template void TetrahedralMesh<1,2>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
//template void TetrahedralMesh<2,3>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
//template void TetrahedralMesh<1,3>::ImportFromMesher<triangulateio>(triangulateio&, unsigned, int *, unsigned, int *, int *);
//template void TetrahedralMesh<2,3>::ImportFromMesher<triangulateio>(triangulateio&, unsigned, int *, unsigned, int *, int *);
//template void TetrahedralMesh<1,2>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned, int *, unsigned, int *, int *);
#ifdef _MSC_VER
//MSVC requires these
template void TetrahedralMesh<2, 2>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
template void TetrahedralMesh<3, 3>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
template void TetrahedralMesh<1, 3>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
template void TetrahedralMesh<1, 1>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
template void TetrahedralMesh<1, 2>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
template void TetrahedralMesh<2, 3>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
template void TetrahedralMesh<1, 3>::ImportFromMesher<triangulateio>(triangulateio&, unsigned, int*, unsigned, int*, int*);
template void TetrahedralMesh<2, 3>::ImportFromMesher<triangulateio>(triangulateio&, unsigned, int*, unsigned, int*, int*);
template void TetrahedralMesh<1, 2>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned, int*, unsigned, int*, int*);
#endif
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TetrahedralMesh)
