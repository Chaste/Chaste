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


#ifndef MONODOMAINTISSUE_HPP_
#define MONODOMAINTISSUE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

#include <vector>
#include "AbstractCardiacTissue.hpp"


/**
 *  MonodomainTissue class.
 *
 *  Essentially identical to AbstractCardiacTissue - see documentation for
 *  AbstractCardiacTissue.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class MonodomainTissue : public virtual AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>
{
private:
    friend class TestMonodomainTissue;

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
        archive & boost::serialization::base_object<AbstractCardiacTissue<ELEMENT_DIM, SPACE_DIM> >(*this);
    }


public:
    /**
     *  Constructor
     *
     * @param pCellFactory  Provides the mesh and cells
     * @param exchangeHalos used in state-variable interpolation.  Defaults to false.
     */
    MonodomainTissue(AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory, bool exchangeHalos=false);

    /**
     * Another constructor (for archiving)
     *
     * @param pMesh the mesh (recovered from archive)
     */
    MonodomainTissue(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);
};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp" // Must be last
EXPORT_TEMPLATE_CLASS2(MonodomainTissue, 1, 1)
EXPORT_TEMPLATE_CLASS2(MonodomainTissue, 1, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainTissue, 1, 3)
EXPORT_TEMPLATE_CLASS2(MonodomainTissue, 2, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainTissue, 3, 3)

namespace boost
{
namespace serialization
{

template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const MonodomainTissue<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    const AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* p_mesh = t->pGetMesh();
    ar & p_mesh;

    // CreateIntracellularConductivityTensor() is called by constructor and uses HeartConfig. So make sure that it is
    // archived too (needs doing before construction so appears here instead of usual archive location).
    HeartConfig* p_config = HeartConfig::Instance();
    ar & *p_config;
    ar & p_config;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor)
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, MonodomainTissue<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* p_mesh;
    ar & p_mesh;

    // CreateIntracellularConductivityTensor() is called by AbstractCardiacTissue constructor and uses HeartConfig.
    // So make sure that it is archived too (needs doing before construction so appears here instead of usual archive location).
    HeartConfig* p_config = HeartConfig::Instance();
    ar & *p_config;
    ar & p_config;

    ::new(t)MonodomainTissue<ELEMENT_DIM, SPACE_DIM>(p_mesh);
}
}
} // namespace ...


#endif /*MONODOMAINTISSUE_HPP_*/
