/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ABSTRACTCHASTEREGION_HPP_
#define ABSTRACTCHASTEREGION_HPP_

#include <cassert>
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "ChasteSerializationVersion.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>

#include "ChastePoint.hpp"

/**
 * Abstract base class for Chaste regions.
 */

template <unsigned SPACE_DIM>
class AbstractChasteRegion
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
        // No member variables, but this is here so boost is happy serializing these classes.
    }
public:

    /**
     * Constructor
     */
    AbstractChasteRegion()
    {}

    /**
     * \todo Proper memory management
     * Cleans any data which the concrete class may have created (archiving)
     */
    virtual void Destroy()
    {
    }

    /**
     * Virtual functions, so virtual destructor.
     */
    virtual ~AbstractChasteRegion()
    {
    }

    /**
     * Checks whether the Chaste point is contained in the region. implemented in the concrete classes
     *
     * @param rPointToCheck Point to be checked to be contained in the region
     * @return true if the point is contained, false otherwise
     */

    virtual bool DoesContain(const ChastePoint<SPACE_DIM>& rPointToCheck) const = 0;

};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractChasteRegion)

namespace boost {
namespace serialization {
/**
 * Specify a version number for archive backwards compatibility.
 *
 * This is how to do BOOST_CLASS_VERSION(AbstractChasteRegion, 1)
 * with a templated class.
 */
template <unsigned SPACE_DIM>
struct version<AbstractChasteRegion<SPACE_DIM> >
{
    /** Version number */
    CHASTE_VERSION_CONTENT(1);
};
} // namespace serialization
} // namespace boost

#endif /*ABSTRACTCHASTEREGION_HPP_*/
