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


#ifndef BIDOMAINTISSUE_HPP_
#define BIDOMAINTISSUE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

#include "AbstractCardiacTissue.hpp"
#include "AbstractConductivityTensors.hpp"




/**
 * BidomainTissue class.
 *
 * See documentation for AbstractCardiacTissue. This class also has extracellular
 * conductivity tensors.
 *
 */
template <unsigned SPACE_DIM>
class BidomainTissue : public virtual AbstractCardiacTissue<SPACE_DIM>
{
private:
    friend class TestBidomainTissue; // for testing.

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
        archive & boost::serialization::base_object<AbstractCardiacTissue<SPACE_DIM> >(*this);
        // Conductivity tensors are dealt with by HeartConfig, and the caches get regenerated.
    }

    /** Extracellular conductivity tensors. */
    AbstractConductivityTensors<SPACE_DIM,SPACE_DIM> *mpExtracellularConductivityTensors;

    /**
     * Convenience method for extracellular conductivity tensors creation
     */
    void CreateExtracellularConductivityTensors();

public:
    /**
     * Constructor sets up extracellular conductivity tensors.
     * @param pCellFactory factory to pass on to the base class constructor
     * @param exchangeHalos used in state-variable interpolation.  Defaults to false.
     */
    BidomainTissue(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, bool exchangeHalos=false);

    /**
     * Archiving constructor
     *
     * @param pMesh  a pointer to the AbstractTetrahedral mesh (recovered from archive).
     */
    BidomainTissue(AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh);

    /**
     * Destructor
     */
    ~BidomainTissue();

    /**
     * Get the extracellular conductivity tensor for the given element
     * @param elementIndex  index of the element of interest
     */
     const c_matrix<double, SPACE_DIM, SPACE_DIM>& rGetExtracellularConductivityTensor(unsigned elementIndex);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BidomainTissue)

namespace boost
{
namespace serialization
{

template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const BidomainTissue<SPACE_DIM> * t, const unsigned int file_version)
{
    const AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* p_mesh = t->pGetMesh();
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
template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, BidomainTissue<SPACE_DIM> * t, const unsigned int file_version)
{
    AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* p_mesh;
    ar & p_mesh;

    // CreateIntracellularConductivityTensor() is called by AbstractCardiacTissue constructor and uses HeartConfig.
    // (as does CreateExtracellularConductivityTensor). So make sure that it is
    // archived too (needs doing before construction so appears here instead of usual archive location).
    HeartConfig* p_config = HeartConfig::Instance();
    ar & *p_config;
    ar & p_config;

    ::new(t)BidomainTissue<SPACE_DIM>(p_mesh);
}
}
} // namespace ...


#endif /*BIDOMAINTISSUE_HPP_*/
