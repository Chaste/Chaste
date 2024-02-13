/*

Copyright (c) 2005-2023, University of Oxford.
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
#ifndef SEMELEMENT_HPP_
#define SEMELEMENT_HPP_

#include "AbstractElement.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <map>

/**
 * An element class for use in the SemMesh class.
 */
template<unsigned DIM>
class SemElement : public AbstractElement<DIM, DIM>
{
private:

    /**
     * The id of the cell
    */
    unsigned mCellId;
  
    /**
     * A collection of interaction layers. Each layer is identified by a name.
     * The map stores a vector containing the indices of the nodes which make up the layer.
     * SemForces produce interactions between or within named layers.
    */
    std::map<std::string, std::vector<unsigned int>> mInteractionLayers;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractElement<DIM, DIM> >(*this);
    }

public:

    /**
     * Default constructor, which doesn't add any nodes: they must be added later.
     *
     * @param index global index of the element
     */
    SemElement(unsigned index);

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index global index of the element
     * @param rNodes vector of Nodes associated with the element
     */
    SemElement(unsigned index,
               const std::vector<Node<DIM>*>& rNodes);
    

    /**
     * Destructor.
     */
    ~SemElement();
    
    /**
     * Sets the cell id of the element
     * 
     * @param id the new id
    */
    void SetCellId(unsigned int id);
    
    /**
     * Adds a new interaction layer
     * 
     * @param layerName the name of the layer
     * @param nodes indices of the nodes which make up the layer
    */
    void AddInteractionLayer(const std::string layerName, std::vector<unsigned int>& nodeIndices);
    
    /**
     * Gets a const reference to the interaction layers
     * 
     * @return a const reference to the interaction layers
    */
    const std::map<std::string, std::vector<unsigned int>>& rGetInteractionLayers();
    
    std::vector<Node<DIM>*>& rGetNodes();
    
    void UpdateNode(const unsigned& rIndex, Node<DIM>* pNode) override;
    void MarkAsDeleted() override;
    void RegisterWithNodes() override;
};

#endif /*SEMELEMENT_HPP_*/