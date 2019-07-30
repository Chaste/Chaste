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
