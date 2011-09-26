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

#ifndef CABASEDSIMULATION_HPP_
#define CABASEDSIMULATION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "AbstractCaUpdateRule.hpp"


/**
 * A lattice-based cell-based simulation object.
 */
 template<unsigned DIM>
class CaBasedSimulation : public AbstractCellBasedSimulation<DIM>
{
    // Allow tests to access private members, in order to test computation of private functions e.g. DoCellBirth()
    friend class TestCaBasedSimulationWithCaBasedCellPopulation;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulation<DIM> >(*this);

        archive & mUpdateRuleCollection;
        archive & mIterateRandomlyOverUpdateRuleCollection;
        archive & mIterateRandomlyOverCells;
    }

    /** Helper member that is a static cast of the cell population. */
    CaBasedCellPopulation<DIM>* mpStaticCastCellPopulation;

    /** The rules used to determine the new locations of cells. */
    std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > > mUpdateRuleCollection;

    /**
     * Whether to iterate randomly over mUpdateRuleCollection when updating cell locations.
     * Defaults to false in the constructor.
     */
    bool mIterateRandomlyOverUpdateRuleCollection;

    /**
     * Whether to iterate randomly over cells when updating cell locations.
     * Defaults to false in the constructor.
     */
    bool mIterateRandomlyOverCells;

    /**
     * Overridden CalculateCellDivisionVector() method.
     *
     * @param pParentCell the parent cell
     */
    c_vector<double, DIM> CalculateCellDivisionVector(CellPtr pParentCell);

    /**
     * Move each cell to a new lattice site for this timestep by calling the
     * CellPopulation method MoveCell().
     */
    void UpdateCellLocations();

    /**
     * Overridden UpdateCellPopulation() method.
     * Calls the deaths, births and (if mUpdateCellPopulation is true) CellPopulation::Update() methods.
     * Does nothing if at the start of a simulation that has just been loaded, to ensure consistency
     * in random number generation.
     */
    void UpdateCellPopulation();

    /**
     * Overridden UpdateCellLocationsAndTopology() method.
     */
    void UpdateCellLocationsAndTopology();

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells whether to initialise cells (defaults to true, set to
     *        false when loading from an archive)
     */
    CaBasedSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                      bool deleteCellPopulationInDestructor=false,
                      bool initialiseCells=true);

    /**
     * @return const reference to mUpdateRuleCollection (used in archiving).
     */
    const std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > > rGetUpdateRuleCollection() const;

    /**
     * Outputs simulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses. Currently not used for Ca Based Simulations.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);

    /**
     * Add an update rule to be used in this simulation.
     *
     * @param pUpdateRule shared pointer to a CA update rule law
     */
    void AddUpdateRule(boost::shared_ptr<AbstractCaUpdateRule<DIM> > pUpdateRule);

    /**
     * Set mIterateRandomlyOverUpdateRuleCollection.
     * 
     * @param iterateRandomly whether to iterate randomly over mUpdateRuleCollection
     */
    void SetIterateRandomlyOverUpdateRuleCollection(bool iterateRandomly);

    /**
     * Set mIterateRandomlyOverUpdateRuleCollection.
     * 
     * @param iterateRandomly whether to iterate randomly over cells
     */
    void SetIterateRandomlyOverCells(bool iterateRandomly);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CaBasedSimulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CaBasedSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CaBasedSimulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM> * p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a CaBasedSimulation.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CaBasedSimulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)CaBasedSimulation<DIM>(*p_cell_population, true, false);
}
}
} // namespace

#endif /*CABASEDSIMULATION_HPP_*/
