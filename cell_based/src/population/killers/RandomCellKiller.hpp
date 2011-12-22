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

#ifndef RANDOMCELLKILLER_HPP_
#define RANDOMCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "RandomNumberGenerator.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


/**
 * A cell killer that randomly kills cells based on the user set probability.
 *
 * The probability passed into the constructor will be the probability
 * of any cell dying whenever TestAndLabelCellsForApoptosis() is called.
 *
 * Note this does take into account timesteps - the input probability is the
 * probability that in an hour's worth of trying, the cell killer will have
 * successfully killed a given cell. In the method TestAndLabelSingleCellForApoptosis()
 * this probability is used to calculate the probability that the cell is killed
 * at a given time step.
 */
template<unsigned DIM>
class RandomCellKiller : public AbstractCellKiller<DIM>
{
private:

    /**
     * Probability that in an hour's worth of trying, the cell killer
     * will have successfully killed a given cell.
      */
     double mProbabilityOfDeathInAnHour;

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

        // Make sure the random number generator is also archived
        SerializableSingleton<RandomNumberGenerator>* p_rng_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_rng_wrapper;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param probabilityOfDeathInAnHour probability that a cell is labelled for apoptosis in one hour's worth of trying
     */
    RandomCellKiller(AbstractCellPopulation<DIM>* pCellPopulation, double probabilityOfDeathInAnHour);

    /**
     * @return mProbabilityOfDeathInAnHour.
     */
    double GetDeathProbabilityInAnHour() const;

    /**
     * Overridden method to test a given cell for apoptosis.
     *
     * @param pCell the cell to test for apoptosis
     */
    void TestAndLabelSingleCellForApoptosis(CellPtr pCell);

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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a RandomCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const RandomCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
    double prob = t->GetDeathProbabilityInAnHour();
    ar << prob;
}

/**
 * De-serialize constructor parameters and initialise a RandomCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, RandomCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
    double prob;
    ar >> prob;

    // Invoke inplace constructor to initialise instance
    ::new(t)RandomCellKiller<DIM>(p_cell_population, prob);
}
}
} // namespace ...

#endif /*RANDOMCELLKILLER_HPP_*/
