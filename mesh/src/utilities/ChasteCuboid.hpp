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


#ifndef CHASTECUBOID_HPP_
#define CHASTECUBOID_HPP_
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractChasteRegion.hpp"
#include "ChastePoint.hpp"


/**
 * This class defines a 3D cuboid and provides a method to check
 * if a given point is contained in the volume.
 */

template <unsigned SPACE_DIM>
class ChasteCuboid : public AbstractChasteRegion<SPACE_DIM>
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
    /** Lower vertex of the cuboid. */
    ChastePoint<SPACE_DIM> mLowerCorner;

    /** Upper vertex of the cuboid.  The space-diagonal opposite corner of mLowerCorner. */
    ChastePoint<SPACE_DIM> mUpperCorner;

public:
    /**
     * The cuboid is defined by any of its two space-diagonal opposite corners.
     *
     * @param rLowerPoint Lower vertex of the cuboid.
     * @param rUpperPoint Upper vertex of the cuboid.
     */
    ChasteCuboid(ChastePoint<SPACE_DIM>& rLowerPoint, ChastePoint<SPACE_DIM>& rUpperPoint);

    /**
     * @return true if a given point is contained in the cuboid.
     *
     * @param rPointToCheck Point to be checked to be contained in the cuboid.
     */
    bool DoesContain(const ChastePoint<SPACE_DIM>& rPointToCheck) const;

    /** @return the upper vertex of the cuboid */
    const ChastePoint<SPACE_DIM>& rGetUpperCorner() const;

    /** @return the lower vertex of the cuboid */
    const ChastePoint<SPACE_DIM>& rGetLowerCorner() const;

    /**
     * @param rDimension dimension
     * @return the width in a particular dimension */
    double GetWidth(unsigned rDimension) const;

    /**
     * @return the longest axis of the cuboid (<SPACE_DIM).
     * In the case where there's a tie, then any of the longest
     * axes are returned.
     */
     unsigned GetLongestAxis() const;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChasteCuboid)

namespace boost
{
namespace serialization
{

template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const ChasteCuboid<SPACE_DIM> * t, const unsigned int file_version)
{
    const ChastePoint<SPACE_DIM>* p_upper_corner =  &(t->rGetUpperCorner());
    const ChastePoint<SPACE_DIM>* p_lower_corner =  &(t->rGetLowerCorner());
    ar & p_upper_corner;
    ar & p_lower_corner;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor)
 */
template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, ChasteCuboid<SPACE_DIM> * t, const unsigned int file_version)
{
    ChastePoint<SPACE_DIM>* p_upper_corner;
    ChastePoint<SPACE_DIM>* p_lower_corner;

    ar & p_upper_corner;
    ar & p_lower_corner;

    ::new(t)ChasteCuboid<SPACE_DIM>((*p_lower_corner), (*p_upper_corner));

    delete p_upper_corner;
    delete p_lower_corner;
}
}
} // namespace ...

#endif /*CHASTECUBOID_HPP_*/
