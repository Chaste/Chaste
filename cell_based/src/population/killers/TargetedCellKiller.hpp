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

#ifndef TARGETEDCELLKILLER_HPP_
#define TARGETEDCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Simple cell killer, which at the first timestep kills any cell
 * whose corresponding location index is a given number.
 */
template<unsigned DIM>
class TargetedCellKiller : public AbstractCellKiller<DIM>
{
private:

    /**
     * The index of the cell to kill
     */
    unsigned mTargetIndex;

    /**
     * Variable to reack when the cell has been killed.
     * Once the cell has been called mBloodLust will stop the killer killing more cells.
     */
    bool mBloodLust;

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
        // archive & mTargetIndex; // done in load_construct_data
        // archive & mBloodlust; // done in load_construct_data
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param targetedIndex The index of the cell to kill
     * @param bloodLust Wether to kill cells or not defaults to true (used by load methods)
     */
    TargetedCellKiller(AbstractCellPopulation<DIM>* pCellPopulation, unsigned targetedIndex, bool bloodLust = true);

    /**
     * @return mTargetIndex.
     */
    unsigned GetTargetIndex() const;

    /**
     * @return mBloodLust.
     */
    unsigned GetBloodLust() const;

    /**
     * Loop over cells and start apoptosis randomly, based on the user-set
     * probability.
     */
    void TestAndLabelCellsForApoptosisOrDeath();

    /**
     * Overridden OutputCellKillerParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TargetedCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a RandomCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const TargetedCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
    unsigned targeted_index = t->GetTargetIndex();
    ar << targeted_index;
    bool blood_lust = t->GetBloodLust();
    ar << blood_lust;
}

/**
 * De-serialize constructor parameters and initialise a RandomCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, TargetedCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
    unsigned targeted_index;
    ar >> targeted_index;
    bool blood_lust;
    ar >> blood_lust;

    // Invoke inplace constructor to initialise instance
    ::new(t)TargetedCellKiller<DIM>(p_cell_population, targeted_index, blood_lust);
}
}
} // namespace ...

#endif /*TARGETEDCELLKILLER_HPP_*/
