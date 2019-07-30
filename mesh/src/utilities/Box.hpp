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
#ifndef BOX_HPP_
#define BOX_HPP_

#include <set>

#include "UblasVectorInclude.hpp"
#include "Node.hpp"
#include "Element.hpp"


/**
 * A small class for a nD 'box' defined by its min/max x/y/z values which
 * can contains a list of nodes and elements located in that box
 */
template<unsigned DIM>
class Box
{
private:

    /** Nodes contained in this box. */
    std::set< Node<DIM>* > mNodesContained;

    /** Elements contained in this box. */
    std::set< Element<DIM,DIM>* > mElementsContained;

public:

    /**
     * Add a node to this box.
     * @param pNode address of the node to be added
     */
    void AddNode(Node<DIM>* pNode);

    /**
     * Remove a node from this box.
     * @param pNode address of the node to be removed
     */
    void RemoveNode(Node<DIM>* pNode);

    /**
     * Remove all nodes from the box.
     */
    void ClearNodes();

    /**
     * An element to this box.
     * @param pElement address of the element to be added
     */
    void AddElement(Element<DIM,DIM>* pElement);

    /** @return all the nodes in this box. */
    std::set< Node<DIM>* >& rGetNodesContained();

    /** @return all the elements in this box. */
    std::set< Element<DIM,DIM>* >& rGetElementsContained();
};

#endif /*BOX_HPP_*/
