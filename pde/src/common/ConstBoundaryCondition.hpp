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

#ifndef _CONSTBOUNDARYCONDITION_HPP_
#define _CONSTBOUNDARYCONDITION_HPP_

#include "AbstractBoundaryCondition.hpp"

#include <boost/serialization/base_object.hpp>

/**
 * Boundary condition that takes a constant value wherever it is applied.
 */
template<unsigned SPACE_DIM>
class ConstBoundaryCondition : public AbstractBoundaryCondition<SPACE_DIM>
{
private:

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
        archive & boost::serialization::base_object<AbstractBoundaryCondition<SPACE_DIM> >(*this);
    }

    /** The constant value of the boundary condition. */
    double mValue;

public:

    /**
     * Create a new boundary condition object.
     *
     * @param value The value of this boundary condition at all points where it
     *    is applied.
     */
    ConstBoundaryCondition(const double value);

    /**
     * @param rX The point at which this boundary condition is to be evaluated.
     * @return The constant value given in the constructor.
     */
    double GetValue(const ChastePoint<SPACE_DIM>& rX) const;
};

#include "SerializationExportWrapper.hpp" // Must be last
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ConstBoundaryCondition)

namespace boost
{
namespace serialization
{

template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const ConstBoundaryCondition<SPACE_DIM> * t, const unsigned int file_version)
{
    const ChastePoint<SPACE_DIM> p;
    const double value = t->GetValue(p);

    ar & value;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor)
 */
template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, ConstBoundaryCondition<SPACE_DIM> * t, const unsigned int file_version)
{
    double value;
    ar & value;

    ::new(t)ConstBoundaryCondition<SPACE_DIM>(value);
}
}
} // namespace ...

#endif //_CONSTBOUNDARYCONDITION_HPP_
