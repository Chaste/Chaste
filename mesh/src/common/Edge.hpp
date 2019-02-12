//
// Created by twin on 25/01/19.
//

#ifndef EDGE_HPP_
#define EDGE_HPP_

#include <set>
#include <vector>

#include "Node.hpp"

typedef std::pair<unsigned ,unsigned> UIndexPair;

/**
 * An edge in a finite element mesh
 */
template<unsigned SPACE_DIM>
class Edge {

private:

    /** Index of this edge within the mesh **/
    unsigned mIndex;

    bool mIsDeleted;

    /** Nodes that form this edge **/
    std::vector<Node<SPACE_DIM>*> mNodes;


    /** Elements that this edge belongs to **/
    std::set<unsigned> mElementIndices;

public:

    Edge(unsigned index) : mIndex(index), mIsDeleted(false)
    {
        this->mIndex = index;
    }

    Edge(unsigned index, Node<SPACE_DIM>* node0, Node<SPACE_DIM>* node1) : mIndex(index), mIsDeleted(false)
    {
        this->SetNodes(node0, node1);
    }

    ~Edge()
    {
        //Remove all previous nodes and references
        for(auto node: mNodes)
            node->RemoveEdge(this->GetIndex());
        mNodes.clear();

    }

    void MarkDeleted()
    {
        mIsDeleted = true;
    }

    bool IsDeleted()
    {
        return mIsDeleted;
    }

    void SetIndex(unsigned index)
    {
        mIndex = index;
    }

    unsigned GetIndex()
    {
        return mIndex;
    }

    UIndexPair GetMapIndex()
    {
        assert(mNodes.size() == 2);
        auto index0 = mNodes[0]->GetIndex();
        auto index1 = mNodes[1]->GetIndex();
        if(index0 > index1)
        {
            auto indexSwap = index0;
            index0 = index1;
            index1 = indexSwap;
        }

        return UIndexPair(index0, index1);
    }

    void RemoveNodes(){

        //Remove all previous nodes and references
        for(auto node: mNodes)
            node->RemoveEdge(this->GetIndex());
        mNodes.clear();
    }

    void SetNodes(Node<SPACE_DIM>* node0, Node<SPACE_DIM>* node1)
    {
        //Clear the nodes first
        this->RemoveNodes();

        //Add nodes
        mNodes.push_back(node0);
        mNodes.push_back(node1);

        for(auto node: mNodes)
            node->AddEdge(this->GetIndex());

    }

    void ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode)
    {
        for(unsigned i = 0 ; i < mNodes.size(); i++)
        {
            if(mNodes[i] == pOldNode)
            {
                pOldNode->RemoveEdge(this->mIndex);
                pNewNode->AddEdge(this->mIndex);
                mNodes[i] = pNewNode;
            }
        }
    }

    Node<SPACE_DIM>* GetNode(unsigned index)
    {
        return mNodes[index];
    }

    unsigned GetNumNodes()
    {
        return mNodes.size();
    }

    bool ContainsNode(Node<SPACE_DIM>* pNode)
    {
        for(auto node: mNodes)
        {
            if(node->GetIndex() == pNode->GetIndex())
                return true;
        }

        return false;
    }

    void AddElement(unsigned elementIndex)
    {
        mElementIndices.insert(elementIndex);
    }

    void RemoveElement(unsigned elementIndex)
    {
        mElementIndices.erase(elementIndex);
    }

    void GetNeighbouringElementIndices(std::set<unsigned>& neighbouring_element_indices)
    {

        neighbouring_element_indices.clear();

        // Loop over nodes owned by this element
        for (auto p_node: mNodes)
        {
            // Find the indices of the elements owned by this node
            std::set<unsigned> containing_elem_indices = p_node->rGetContainingElementIndices();

            // Form the union of this set with the current set of neighbouring element indices
            std::set<unsigned> all_elements;
            std::set_union(neighbouring_element_indices.begin(), neighbouring_element_indices.end(),
                           containing_elem_indices.begin(), containing_elem_indices.end(),
                           std::inserter(all_elements, all_elements.begin()));

            // Update the set of neighbouring element indices
            neighbouring_element_indices = all_elements;
        }
    }

    unsigned GetNumElements()
    {
        return mElementIndices.size();
    }

    bool IsEdgeValid()
    {
        //MUST have 2 existing nodes to form an edge
        if(mNodes.size() != 2)
            return false;


        //Nodes should not be nullptr
        for(auto node: mNodes)
        {
            if(node == nullptr)
                return false;
        }

        //Can't have associated elements if we're less than 2D
        if(SPACE_DIM <= 1 && mElementIndices.size() > 0)
        {
            return false;
        }

        //An ege can only have a maximum of two elements in 2D
        if(SPACE_DIM == 2 && mElementIndices.size() > 2)
        {
            return false;
        }


        return true;
    }



};


#endif //EDGE_HPP_
