/*

Copyright (C) University of Oxford, 2005-2012

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
