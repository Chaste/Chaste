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

#ifndef SLOUGHINGCELLKILLER_HPP_
#define SLOUGHINGCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 *  A cell killer that kills cells if they are outside the domain.
 *  The domain is assumed to start at x=0 and y=0. By default only cells
 *  are sloughed if y>mSloughLength. To slough the sides call the constructor
 *  with the appropriate parameter.
 */
template<unsigned DIM>
class SloughingCellKiller : public AbstractCellKiller<DIM>
{
private:

    /** Whether cells should be sloughed from the sides of the crypt. */
    bool mSloughSides;

    /**
     * The height of the domain, non-dimensionalised with cell length.
     * This parameter determines when cells are sloughed from the domain.
     */
    double mSloughHeight;

    /**
    * The width of the domain, non-dimensionalised with cell length.
    * This determines when cells are sloughed from sides of the domain in 2D.
    */
    double mSloughWidth;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     *
     * @param pCrypt pointer to a cell population
     * @param sloughHeight the height at which to slough from the domain
     * @param sloughSides whether to slough cells at the side of the domain
     * @param sloughWidth the width of the domain (note slough on left and right)
     */
    SloughingCellKiller(AbstractCellPopulation<DIM>* pCrypt, double sloughHeight, bool sloughSides=false, double sloughWidth = 10.0);

    /**
     * @return mSloughSides.
     */
    bool GetSloughSides() const;

    /**
     * @return mSloughHeight.
     */
    double GetSloughHeight() const;

    /**
     * @return mSloughWidth.
     */
    double GetSloughWidth() const;

    /**
     * Loops over cells and kills cells outside boundary.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath();

    /**
     * Outputs cell killer parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SloughingCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SloughingCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SloughingCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_crypt = t->GetCellPopulation();
    ar << p_crypt;
    bool slough_sides = t->GetSloughSides();
    ar << slough_sides;
    double slough_height = t->GetSloughHeight();
    ar << slough_height;
    double slough_width = t->GetSloughWidth();
    ar << slough_width;
}

/**
 * De-serialize constructor parameters and initialise a SloughingCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SloughingCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_crypt;
    ar >> p_crypt;
    bool slough_sides;
    ar >> slough_sides;
    double slough_height;
    ar >> slough_height;
    double slough_width;
    ar >> slough_width;

    // Invoke inplace constructor to initialise instance
    ::new(t)SloughingCellKiller<DIM>(p_crypt, slough_height, slough_sides, slough_width);
}
}
} // namespace ...

#endif /*SLOUGHINGCELLKILLER_HPP_*/
