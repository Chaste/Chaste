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

#ifndef EDGE_HPP_
#define EDGE_HPP_

#include <set>
#include <vector>
#include "ChasteSerialization.hpp"
#include <boost/serialization/set.hpp>
#include "Node.hpp"

typedef std::pair<unsigned ,unsigned> UIndexPair;

/**
 * An edge in a finite element mesh.
 */
template<unsigned SPACE_DIM>
class Edge
{
private:

    /** Index of this edge within the mesh **/
    unsigned mIndex;

    /** Indicates whether this edge is deleted from the mesh **/
    bool mIsDeleted;

    /** Nodes that form this edge **/
    std::vector<Node<SPACE_DIM>*> mNodes;

    /** Elements that this edge belongs to **/
    std::set<unsigned> mElementIndices;
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mElementIndices;
    }
public:

    /**
     * Create an Edge with only the local index within the mesh. Nodes must be set using SetNode() after the
     * constructor to make this a valid edge.
     *
     * @param index Index of this edge within the mesh
     */
    Edge(unsigned index);

    /**
     * Create an Edge that has an index and associated nodes.
     *
     * @param index Index of this edge within the mesh
     * @param pNode0
     * @param pNode1
     */
    Edge(unsigned index, Node<SPACE_DIM>* pNode0, Node<SPACE_DIM>* pNode1);

    /**
     * Destructor.
     */
    ~Edge();

    /**
     * Mark the Edge to be deleted
     */
    void MarkAsDeleted();

    /**
     *
     * @return True if Edge has been marked as deleted
     */
    bool IsDeleted();

    /**
     * Sets the index of this edge within the mesh
     * @param index The index of this edge within the mesh
     */
    void SetIndex(unsigned index);

    /**
     * Gets the index of this edge within the mesh
     * @return
     */
    unsigned GetIndex() const;

    /**
     * Obtains a pair of associated nodes' indices
     * @return Obtains a pair of associated nodes' indices
     */
    UIndexPair GetMapIndex();


    /**
     * Clear all associated nodes
     */
    void RemoveNodes();

    /**
     * Set the Edge's associated nodes
     * @param pNode0 A Node that forms one point of the edge
     * @param pNode1 A different Node that forms the other point of the edge
     */
    void SetNodes(Node<SPACE_DIM>* pNode0, Node<SPACE_DIM>* pNode1);

    /**
     * Replace a Node in this Edge with another
     * @param pOldNode The old Node to be replaced
     * @param pNewNode New Node to replace the old Node
     */
    void ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode);

    /**
     * Gets the Node at index
     * @param index
     * @return
     */
    Node<SPACE_DIM>* GetNode(unsigned index) const;

    /**
     * Gets the number of Nodes associated with this edge
     * @return
     */
    unsigned GetNumNodes();

    /**
     * Checks that the Edge contains pNode
     * @param pNode
     * @return true if pNode is containd in Edge, otherwise false
     */
    bool ContainsNode(Node<SPACE_DIM>* pNode) const;

    /**
     *
     * @return The centre location of this edge, which is middle of two Nodes associated with the Edge
     */
    c_vector<double, SPACE_DIM> rGetCentreLocation();

    /**
     *
     * @return The length of the Edge, i.e. the distance between the Edge's two Nodes
     */
    double rGetLength();

    /**
     * Gets other Element indices that the edge is associated to.
     * element.
     * @param elementIndex The Element index to exclude
     * @return A set of Element indices or an empty set if there's no association.
     */
    std::set<unsigned> GetOtherElements(unsigned elementIndex);

    /**
     * Add an Element index that the Edge is associated to
     * @param elementIndex
     */
    void AddElement(unsigned elementIndex);

    /**
     * Remove an Element index association from the Edge
     * @param elementIndex
     */
    void RemoveElement(unsigned elementIndex);

    /**
    * Gets all Element indices that the edge is associated to. Used for testing the accounting of Add and Remove
    * element.
    *
    * @return A set of Element indices or an empty set if there's no association.
    */
    std::set<unsigned> GetNeighbouringElementIndices();

    /**
     *
     * @return The number of Elements associated with this Edge
     */
    unsigned GetNumElements();

    /**
     * Checks whether this edge is valid i.e. must have 2 valid Nodes, can't have associated Elements if 1D
     * and maximum of 2 associated Elements if 2D. A consistency check is also performed for Edge's associated elements
     * and contain Nodes' associated elements.
     * @return
     */
    bool IsEdgeValid();

    /**
     * Checks whether the edge is on the boundary
     * @return true if on boundary
     */
    bool IsBoundaryEdge() const;

    /**
     * Comparison operator.
     * @param edge_1
     * @param edge_2
     * @return true if edges are equal
     */
    bool operator==(const Edge<SPACE_DIM>& edge) const;
};
#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Edge)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Node.
 */
template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const Edge<SPACE_DIM> * t, const unsigned int file_version)
{
    const Node<SPACE_DIM>* const p_Node0 = t->GetNode(0);
    ar & p_Node0;
    const Node<SPACE_DIM>* const p_Node1 = t->GetNode(1);
    ar & p_Node1;

    unsigned index = t->GetIndex();
    ar << index;
}

/**
 * De-serialize constructor parameters and initialize an Edge.
 */
template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, Edge<SPACE_DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance of Edge
    Node<SPACE_DIM>* p_Node0;
    ar & p_Node0;
    Node<SPACE_DIM>* p_Node1;
    ar & p_Node1;

    unsigned index;
    ar >> index;

    // Invoke inplace constructor to initialise instance
    ::new(t)Edge<SPACE_DIM>(index, p_Node0, p_Node1);
}

}
} // namespace ...
#endif //EDGE_HPP_
