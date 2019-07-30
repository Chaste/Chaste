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


#ifndef AIRWAY_TREE_WALKER_HPP_
#define AIRWAY_TREE_WALKER_HPP_

#include <map>
#include "TetrahedralMesh.hpp"

/**
 * Utility class to facilitate traversal of airway trees
 *
 * Chaste represents airway trees using TetrahedralMesh<1,3>. This works well for a number of cases, including
 * finite element simulations, however the lack of meta data to allow easy access to parent/children of each
 * branch complicates relatively simple recursive calculations. Note that this class is only valid
 * on a proper tree and doesn't do much sanity checking - results will be undefined on trees with cycles in.
 */
class AirwayTreeWalker
{
public:
    /**
     * Constructor
     *
     * @param rAirwaysMesh A mesh containing the airway tree to analyze
     * @param rootIndex The index of the inlet node (trachea entrance)
     */
    AirwayTreeWalker(AbstractTetrahedralMesh<1,3>& rAirwaysMesh,
                     unsigned rootIndex);

    /*
     * @return the outlet element index (i.e. the element associated with the root of the tree)
     */
    unsigned GetOutletElementIndex()
    {
        return mOutletElementIndex;
    }

    /*
     * @return true if nodes appear in graph/topological order (a parent should have a lower index than its children)
     */
    bool GetNodesAreGraphOrdered()
    {
        return mNodesAreGraphOrdered;
    }

    /* Gets the parent element of a given element.
     *
     * Returns a null pointer if the parent is the root of the tree
     *
     * @param pElement The element to find the parent of
     * @return The parent element
     */
    Element<1,3>* GetParentElement(Element<1,3>* pElement);

    /**
     * Gets the index of the parent element of a given element.
     *
     * Throws an exception if the element has no parent
     *
     * @param pElement The element to find the parent of
     * @return The index of the parent element
     */
    unsigned GetParentElementIndex(Element<1,3>* pElement);

    /**
     * Gets the index of the parent element of a given element.
     *
     * Throws an exception if the element has no parent
     *
     * @param index The index of the element to find the parent of
     * @return The parent element
     */
    Element<1,3>* GetParentElement(unsigned index);

    /**
     * Gets the index of the parent element of a given element.
     *
     * Throws an exception if the element has no parent
     *
     * @param index The index of the element to find the parent of
     * @return The index of the parent element
     */
    unsigned GetParentElementIndex(unsigned index);

    /**
     * Returns the number of child elements of a given element
     *
     * @param pElement The element to find the number of children of
     * @return The number of child elements
     */
    unsigned GetNumberOfChildElements(Element<1,3>* pElement);

    /**
     * Returns the number of child elements of a given element
     *
     * @param index The index of the element to find the number of children of
     * @return The number of child elements
     */
    unsigned GetNumberOfChildElements(unsigned index);


    /**
     * Returns a vector of pointers to child elements of a given element.
     *
     * @param pElement The element to find the children of
     * @return A vector of child elements
     */
    std::vector<Element<1,3>* > GetChildElements(Element<1,3>* pElement);

    /**
     * Returns a vector of pointers to indices child elements of a given element.
     *
     * @param pElement The element to find the children of
     * @return A vector of indices of child elements
     */
    std::vector<unsigned> GetChildElementIndices(Element<1,3>* pElement);

    /**
     * Returns the distal most node of an element
     *
     * @param pElement The element to find the child of
     * @return The distal node of the element
     */
    Node<3>* GetDistalNode(Element<1,3>* pElement);

    /**
     * Returns the distal most node of an element
     *
     * @param index The index of the element to find the child of
     * @return The distal node of the element
     */
    Node<3>* GetDistalNode(unsigned index);

    /**
     * Returns the index of the distal most node of an element
     *
     * @param pElement The element to find the child of
     * @return The index of the distal node of the element
     */
    unsigned GetDistalNodeIndex(Element<1,3>* pElement);

    /**
     * Returns the index of the distal most node of an element
     *
     * @param index The index of the element to find the child of
     * @return The index of the distal node of the element
     */
    unsigned GetDistalNodeIndex(unsigned index);


    /**
     * Returns the generation number of a given element
     *
     * @param pElement The element to get the generation for
     * @return The generation number of the given element
     */
    unsigned GetElementGeneration(Element<1,3>* pElement);

    /**
     * Returns the generation number of a given element
     *
     * @param element_index The element index to get the generation for
     * @return The generation number of the given element
     */
    unsigned GetElementGeneration(unsigned element_index);

    /**
     * Returns the maximum airway generation
     *
     * @return The maximum generation number of all the elements
     */
    unsigned GetMaxElementGeneration();

    /**
     * Get Element Horsfield Order
     *
     * @param pElement The element to get the Horsfield order for
     * @return The element's Horsfield order
     */
    unsigned GetElementHorsfieldOrder(Element<1,3>* pElement);

    /**
     * Get Element Horsfield Order
     *
     * @param element_index The element index to get the Horsfield order for
     * @return The element's Horsfield order
     */
    unsigned GetElementHorsfieldOrder(unsigned element_index);

    /**
     * Returns the maximum Horsfield order
     *
     * @return The maximum Horsfield order of all the elements
     */
    unsigned GetMaxElementHorsfieldOrder();

    /**
     * Get Element Strahler Order
     *
     * @param pElement The element to get the Strahler order for
     * @return The element's Strahler order
     */
    unsigned GetElementStrahlerOrder(Element<1,3>* pElement);

    /**
     * Get Element Strahler Order
     *
     * @param element_index The element index to get the Strahler order for
     * @return The element's Strahler order
     */
    unsigned GetElementStrahlerOrder(unsigned element_index);

    /**
     * Returns the maximum Strahler order
     *
     * @return The maximum Strahler order of all the elements
     */
    unsigned GetMaxElementStrahlerOrder();

private:

    /** A mesh containing the airways tree.  */
    AbstractTetrahedralMesh<1,3>& mMesh;

    unsigned mOutletNodeIndex; /**< The outlet node is the root of the branching tree structure */
    unsigned mOutletElementIndex; /**< The outlet element is associated to the outlet node and is at the root */
    /** A check that nodes appear in graph/topological order (a parent should have a lower index than its children) */
    bool mNodesAreGraphOrdered;

    /** A map to facilitate finding the parent element of an element */
    std::map<unsigned, unsigned> mParentElementMap;

    /** A map to facilitate finding child elements of an element */
    std::map<unsigned, std::vector<unsigned> > mChildElementsMap;

    /** A map to facilitate finding distal nodes of an element */
    std::map<unsigned, unsigned> mDistalNodeMap;

    /** Maps an element ID to that element's generation number */
    std::map<unsigned, unsigned> mElementGenerations;

    /** Maps an element ID to that element's Horsfield order */
    std::map<unsigned, unsigned> mElementHorsfieldOrder;

    /** Maps an element ID to that element's Strahler order */
    std::map<unsigned, unsigned> mElementStrahlerOrder;

    /**
     * Utility method to recursively process the tree
     *
     * Assumes that the parent of the given element has already been processed.
     *
     * @param pElement The element to process.
     * @param pParentNode The parent node of the element being processed (prevents backtracking)
     */
    void ProcessElement(Element<1,3>* pElement, Node<3>* pParentNode);

    /**
     * Recursively process elements to calculate their properties
     *
     * @param pElement The element to process
     */
    void CalculateElementProperties(Element<1,3>* pElement);
};

#endif // AIRWAY_TREE_WALKER
