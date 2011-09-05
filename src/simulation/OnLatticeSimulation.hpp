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

#ifndef ONLATTICESIMULATION_HPP_
#define ONLATTICESIMULATION_HPP_

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractPottsUpdateRule.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

/**
 * Run an on-lattice 2D or 3D cell-based simulation using a Potts Model
 *
 * The OnLatticeSimulation is constructed with a CellPopulation, which
 * updates the correspondence between each Cell and its spatial representation
 * and handles cell division (governed by the CellCycleModel associated
 * with each cell). Once constructed, one or more Update rules may be passed
 * to the OnLatticeSimulation object, to define the processes which update
 * cells in the CellPopulation. Similarly, one or more CellKillers may be passed
 * to the OnLatticeSimulation object to specify conditions in which Cells
 * may die.
 */
template<unsigned DIM>
class OnLatticeSimulation : public AbstractCellBasedSimulation<DIM>
{

protected:

    /** Helper member that is a static cast of the cell population. */
    PottsBasedCellPopulation<DIM>* mpStaticCastCellPopulation;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<AbstractCellBasedSimulation<DIM> >(*this);
    }

    /**
     * Overridden this UpdateCellLocationsAndTopology method.
     *
     * In On lattice simulations this method performs the monte carlo sampling.
     */
    virtual void UpdateCellLocationsAndTopology();

public:
    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationAndCellKillersInDestructor Whether to delete the cell population and cell killer collection on destruction to free up memory
     * @param initialiseCells Whether to initialise cells (set to false when loading from an archive)
     */
    OnLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                         bool deleteCellPopulationAndCellKillersInDestructor=false,
                         bool initialiseCells=true);

    /**
     * Destructor.
     *
     * This frees the Forces and Boundary Conditions, if they were created by de-serialization.
     */
    virtual ~OnLatticeSimulation()
    {
    }

    /**
     * Add an update rule to be used in this simulation (use this to set the Hamiltonian).
     *
     * @param pUpdateRule pointer to a potts update rule law
     */
    void AddUpdateRule(AbstractPottsUpdateRule<DIM>* pUpdateRule);

    /**
     * Overridden OutputAdditionalSimulationSetup method to output the force and cell
     * population boundary condition information
     */
    void OutputAdditionalSimulationSetup(out_stream& rParamsFile);

    /**
     * Outputs simulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OnLatticeSimulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an OnLatticeSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const OnLatticeSimulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise an OnLatticeSimulation.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, OnLatticeSimulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance, last two variables set extra
    // member variables to be deleted as they are loaded from archive and to not initialise sells.
    ::new(t)OnLatticeSimulation<DIM>(*p_cell_population, true, false);
}
}
} // namespace

#endif /*ONLATTICESIMULATION_HPP_*/
