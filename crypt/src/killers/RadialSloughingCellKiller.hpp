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

#ifndef RADIALSLOUGHINGCELLKILLER_HPP_
#define RADIALSLOUGHINGCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
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
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population.
     * @param centre the centre of death.
     * @param radius the radius of death.
     */
    RadialSloughingCellKiller(AbstractCellPopulation<2>* pCellPopulation,
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
     * Loop over cells and kills cells outside boundary.
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
    Archive & ar, const RadialSloughingCellKiller * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
    ar & p_cell_population;
    c_vector<double,2> centre = t->GetCentre();
    ar & centre[0];
    ar & centre[1];
    double radius = t->GetRadius();
    ar & radius;
}

/**
 * De-serialize constructor parameters and initialise a RadialSloughingCellKiller.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, RadialSloughingCellKiller * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar & p_cell_population;
    c_vector<double,2> centre;
    ar & centre[0];
    ar & centre[1];
    double radius;
    ar & radius;

    // Invoke inplace constructor to initialise instance
    ::new(t)RadialSloughingCellKiller(p_cell_population, centre, radius);
}
}
} // namespace ...


#endif /*RADIALSLOUGHINGCELLKILLER_HPP_*/

