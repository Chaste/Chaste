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

#ifndef LATTICEBASEDCELLBASEDSIMULATION_HPP_
#define LATTICEBASEDCELLBASEDSIMULATION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "CellBasedSimulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "AbstractCaUpdateRule.hpp"


/**
 * A lattice-based cell-based simulation object.
 */
 template<unsigned DIM>
class CaBasedCellBasedSimulation : public CellBasedSimulation<DIM>
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCellBasedSimulationWithCaBasedCellPopulation;

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
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<CellBasedSimulation<DIM> >(*this);

        archive & mUpdateRuleCollection;
        archive & mIterateRandomlyOverUpdateRuleCollection;
        archive & mIterateRandomlyOverCells;
    }

    /** Helper member that is a static cast of the cell population. */
    CaBasedCellPopulation<DIM>* mpStaticCastCellPopulation;

    /** The rules used to determine the new locations of cells. */
    std::vector<AbstractCaUpdateRule<DIM>*> mUpdateRuleCollection;

    /** Whether delete the collection of update rules in the destructor. */
    bool mAllocatedMemoryForUpdateRuleCollection;

    /** Whether to iterate randomly over mUpdateRuleCollection when updating cell locations. */
    bool mIterateRandomlyOverUpdateRuleCollection;

    /** Whether to iterate randomly over cells when updating cell locations. */
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

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param updateRuleCollection The cell movement rules to use in the simulation
     * @param iterateRandomlyOverUpdateRuleCollection whether to iterate randomly over
     *        mUpdateRuleCollection when updating cell locations (defaults to false)
     * @param iterateRandomlyOverCells whether to iterate randomly over cells when
     *        updating cell locations (defaults to false)
     * @param deleteCellPopulationAndForceCollection Whether to delete the cell population and force
     *        collection on destruction to free up memory (defaults to false)
     * @param initialiseCells whether to initialise cells (defaults to true, set to
     *        false when loading from an archive)
     */
    CaBasedCellBasedSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                      std::vector<AbstractCaUpdateRule<DIM>*> updateRuleCollection,
                      bool iterateRandomlyOverUpdateRuleCollection=false,
                      bool iterateRandomlyOverCells=false,
                      bool deleteCellPopulationAndForceCollection=false,
                      bool initialiseCells=true);

    /**
     * Destructor.
     */
    ~CaBasedCellBasedSimulation();

    /**
     * Main solve method.
     *
     * This method sets up the simulation time, creates output files, and initialises the
     * cell population. It then iterates through a time loop. At each time step, first any cell death
     * or birth is implemented, then the cell population topology is updated, then the update rule is
     * recalculated and the cell population evolved according to whatever update rules are present in
     * the simulation, and finally the results for that time step are output to file. At the
     * end of the time loop, the method closes any output files.
     */
    void Solve();

    /**
     * @return const reference to mUpdateRuleCollection (used in archiving).
     */
    const std::vector<AbstractCaUpdateRule<DIM>*> rGetUpdateRuleCollection() const;

    /**
     * Outputs simulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses. Currently not used for Ca Based Simulations.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CaBasedCellBasedSimulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CaBasedCellBasedSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CaBasedCellBasedSimulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM> * p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
    const std::vector<AbstractCaUpdateRule<DIM>*> update_rule_collection = t->rGetUpdateRuleCollection();
    ar & update_rule_collection;
}

/**
 * De-serialize constructor parameters and initialise a CaBasedCellBasedSimulation.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CaBasedCellBasedSimulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
    std::vector<AbstractCaUpdateRule<DIM>*> update_rule_collection;
    ar >> update_rule_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)CaBasedCellBasedSimulation<DIM>(*p_cell_population, update_rule_collection, false, false, true, false);
}
}
} // namespace

#endif /*LATTICEBASEDCELLBASEDSIMULATION_HPP_*/
