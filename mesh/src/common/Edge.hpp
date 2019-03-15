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

    void RemoveNodes()
    {
        mNodes.clear();
    }

    void SetNodes(Node<SPACE_DIM>* node0, Node<SPACE_DIM>* node1)
    {
        //Clear the nodes first
        this->RemoveNodes();

        //Add nodes
        this->mNodes.push_back(node0);
        this->mNodes.push_back(node1);

    }

    void ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode)
    {
        for(unsigned i = 0 ; i < mNodes.size(); i++)
        {
            if(this->mNodes[i] == pOldNode)
            {
                this->mNodes[i] = pNewNode;
                break;
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

    c_vector<double, SPACE_DIM> rGetCentreLocation()
    {
        assert(mNodes.size() == 2);

        return (mNodes[0]->rGetLocation() + mNodes[1]->rGetLocation()) /2.0;
    }



    std::set<unsigned> GetOtherElements(unsigned elementIndex)
    {
        std::set<unsigned> otherElements;
        std::set<unsigned> currentElem;
        currentElem.insert(elementIndex);
        std::set_difference(currentElem.begin(), currentElem.end(),
                mElementIndices.begin(), mElementIndices.end(),
                std::inserter(otherElements, otherElements.begin()));

        return otherElements;
    }

    void AddElement(unsigned elementIndex)
    {
        mElementIndices.insert(elementIndex);
    }

    void RemoveElement(unsigned elementIndex)
    {
        mElementIndices.erase(elementIndex);
    }

    std::set<unsigned> GetNeighbouringElementIndices()
    {

        std::set<unsigned> neighbouring_element_indices;

        auto elem_indices0 = mNodes[0]->rGetContainingElementIndices();
        auto elem_indices1 = mNodes[1]->rGetContainingElementIndices();

        std::set_intersection(elem_indices0.begin(), elem_indices0.end(),
                elem_indices1.begin(), elem_indices1.end(),
                std::inserter(neighbouring_element_indices, neighbouring_element_indices.begin()));

        return neighbouring_element_indices;
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

        auto neighbour_indices = GetNeighbouringElementIndices();
        if(neighbour_indices != mElementIndices)
        {
            return false;
        }

        return true;
    }



};


#endif //EDGE_HPP_
