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
    void CheckAndLabelCellsForApoptosisOrDeath();

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
 * Serialize information required to construct a TargetedCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const TargetedCellKiller<DIM> * t, const unsigned int file_version)
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
 * De-serialize constructor parameters and initialise a TargetedCellKiller.
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
