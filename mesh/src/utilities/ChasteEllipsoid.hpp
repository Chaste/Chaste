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


#ifndef CHASTEELLIPSOID_HPP_
#define CHASTEELLIPSOID_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractChasteRegion.hpp"
#include "ChastePoint.hpp"


/**
 * This class defines a 3D ellipsoid and provides a method to check
 * if a given point is contained in the volume.
 */
template <unsigned SPACE_DIM>
class ChasteEllipsoid : public AbstractChasteRegion<SPACE_DIM>
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
    /** Centre of the ellipsoid. */
    ChastePoint<SPACE_DIM> mCentre;

    /** Radii of the ellipsoid. */
    ChastePoint<SPACE_DIM> mRadii;

public:
    /**
     * The (axis aligned) ellipsoid is defined by its centre and its radii in the x, y and z directions.
     *
     * @param rCentre Centre of the ellipsoid.
     * @param rRadii Radii of the ellipsoid.
     */
    ChasteEllipsoid(ChastePoint<SPACE_DIM>& rCentre, ChastePoint<SPACE_DIM>& rRadii);

    /**
     * @return true if a given point is contained in the ellipsoid.
     *
     * @param rPointToCheck Point to be checked to be contained in the ellipsoid.
     */
    bool DoesContain(const ChastePoint<SPACE_DIM>& rPointToCheck) const;

    /** @return centre of the ellipsoid */
    const ChastePoint<SPACE_DIM>& rGetCentre() const;

    /** @return radii of the ellipsoid */
    const ChastePoint<SPACE_DIM>& rGetRadii() const;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChasteEllipsoid)

namespace boost
{
namespace serialization
{

template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const ChasteEllipsoid<SPACE_DIM> * t, const unsigned int file_version)
{
    const ChastePoint<SPACE_DIM>* p_centre =  &(t->rGetCentre());
    const ChastePoint<SPACE_DIM>* p_radii =  &(t->rGetRadii());
    ar & p_centre;
    ar & p_radii;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor)
 */
template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, ChasteEllipsoid<SPACE_DIM> * t, const unsigned int file_version)
{
    ChastePoint<SPACE_DIM>* p_centre;
    ChastePoint<SPACE_DIM>* p_radii;

    ar & p_centre;
    ar & p_radii;

    ::new(t)ChasteEllipsoid<SPACE_DIM>((*p_centre), (*p_radii));

    delete p_centre;
    delete p_radii;
}
}
} // namespace ...

#endif /*CHASTEELLIPSOID_HPP_*/
