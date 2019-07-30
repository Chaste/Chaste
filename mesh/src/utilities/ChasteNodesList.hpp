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


#ifndef CHASTENODESLIST_HPP_
#define CHASTENODESLIST_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractChasteRegion.hpp"
#include "Node.hpp"
#include "ChastePoint.hpp"

#include <vector>
using namespace std;
/**
 * This class defines a list of nodes and provides a method to check
 * whether a point is contained in the list.
 */
template <unsigned SPACE_DIM>
class ChasteNodesList : public AbstractChasteRegion<SPACE_DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractChasteRegion<SPACE_DIM> >(*this);
    }

private:

    /** A vector to store the list of nodes*/
    std::vector< Node<SPACE_DIM>*> mListOfNodes;

    /** Whether we own the Node objects and should free the memory on destruction */
    bool mOwnNodes;

public:

    /**
     * Constructor
     *
     * @param rNodesList  a standard vector of (pointer to) nodes
     * @param ownNodes  whether we own the Node objects and should free the memory on destruction
     */
    ChasteNodesList(const std::vector<Node<SPACE_DIM>*> rNodesList,
                    bool ownNodes=false);

    /**
     * Clean the memory used by the nodes in this node list
     */
    ~ChasteNodesList();

    /** @return the list of nodes in this nodes list */
    const std::vector< Node<SPACE_DIM>*>& rGetNodesList() const;

    /**
     * @return true if a given point is contained in the node list.
     *
     * @param rPointToCheck Point to be checked whether it is a node in the list.
     */
    bool DoesContain(const ChastePoint<SPACE_DIM>& rPointToCheck) const;

    /**
     * @return the size of the nodes list
     */
    unsigned GetSize() const;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChasteNodesList)

namespace boost
{
namespace serialization
{

template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const ChasteNodesList<SPACE_DIM> * t, const unsigned int file_version)
{
    const std::vector<Node<SPACE_DIM>*>& node_list = t->rGetNodesList();

    // Archive the size first
    unsigned size = t->GetSize();
    ar & size;

    // The ChastePoints have to be at unique addresses, otherwise Boost thinks that the first archive is sufficient
    std::vector <ChastePoint<SPACE_DIM>* > point_list;

    for (unsigned i = 0; i < node_list.size(); i++)
    {
        c_vector<double, SPACE_DIM> loc =  node_list[i]->rGetLocation();
        point_list.push_back(new  ChastePoint<SPACE_DIM>(loc));
        ar & point_list[i];
        unsigned index = node_list[i]->GetIndex();
        ar & index;
    }

    // Clean memory
    for (unsigned i = 0; i < point_list.size(); i++)
    {
        delete point_list[i];
    }
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor)
 */
template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, ChasteNodesList<SPACE_DIM> * t, const unsigned int file_version)
{
    // Unarchive the size first
    unsigned size;
    ar & size;

    // Rebuild the node list based on the unarchived points and indices
    std::vector<Node<SPACE_DIM>* > node_list;
    for (unsigned i = 0; i < size; i++)
    {
        ChastePoint<SPACE_DIM>* p_point;
        unsigned index;
        ar & p_point;
        ar & index;

        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>( index, *(p_point));
        node_list.push_back(p_node);

        delete p_point; // not needed any longer
    }

    ::new(t)ChasteNodesList<SPACE_DIM>(node_list, true);
}
}
} // namespace ...

#endif /*CHASTENODESLIST_HPP_*/
