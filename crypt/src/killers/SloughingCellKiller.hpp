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
    SloughingCellKiller(AbstractCellPopulation<DIM>* pCrypt,
                        double sloughHeight,
                        bool sloughSides = false,
                        double sloughWidth = 10.0);

    /**
     * Destructor
     */
    virtual ~SloughingCellKiller(){};

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
    virtual void CheckAndLabelCellsForApoptosisOrDeath();

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
    Archive & ar, const SloughingCellKiller<DIM> * t, const unsigned int file_version)
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
