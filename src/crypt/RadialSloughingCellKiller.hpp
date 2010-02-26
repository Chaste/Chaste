/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef RADIALSLOUGHINGCELLKILLER_HPP_
#define RADIALSLOUGHINGCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

/**
 *  Radial sloughing cell killer for use with the crypt projection model.
 *
 *  Kills cells if they are outside a circle whose centre and radius can be
 *  passed in but are take default values.
 */
class RadialSloughingCellKiller : public AbstractCellKiller<2>
{
private:

    /** Centre of death. */
    c_vector<double,2> mCentre;

    /** Radius of death. */
    double mRadius;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pTissue pointer to the tissue.
     * @param centre the centre of death.
     * @param radius the radius of death.
     */
    RadialSloughingCellKiller(AbstractTissue<2>* pTissue,
                              c_vector<double,2> centre,
                              double radius);

    /**
     * @return mCentre.
     */
    c_vector<double,2> GetCentre() const;

    /**
     * @return mRadius.
     */
    double GetRadius() const;

    /**
     *  Loop over cells and kills cells outside boundary.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath();

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(RadialSloughingCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a RadialSloughingCellKiller.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const RadialSloughingCellKiller * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<2>* const p_tissue = t->GetTissue();
    ar << p_tissue;
    c_vector<double,2> centre = t->GetCentre();
    ar << centre[0];
    ar << centre[1];
    double radius = t->GetRadius();
    ar << radius;
}

/**
 * De-serialize constructor parameters and initialise a RadialSloughingCellKiller.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, RadialSloughingCellKiller * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<2>* p_tissue;
    ar >> p_tissue;
    c_vector<double,2> centre;
    ar >> centre[0];
    ar >> centre[1];
    double radius;
    ar >> radius;

    // Invoke inplace constructor to initialise instance
    ::new(t)RadialSloughingCellKiller(p_tissue, centre, radius);
}
}
} // namespace ...


#endif /*RADIALSLOUGHINGCELLKILLER_HPP_*/

