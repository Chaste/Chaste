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

    /**
     * Whether to write the cell velocities to a file.
     * Initialised to false in constuctor.
     */
    bool mOutputCellVelocities;

    /** Results file cell velocities. */
    out_stream mpCellVelocitiesFile;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulation<DIM> >(*this);
        archive & mOutputCellVelocities;
    }

    /**
     * Overridden UpdateCellLocationsAndTopology() method.
     *
     * In on-lattice simulations this method performs Monte Carlo sampling.
     */
    void UpdateCellLocationsAndTopology();

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells Whether to initialise cells (defaults to true, set to false when loading
     * from an archive)
     */
    OnLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                        bool deleteCellPopulationInDestructor=false,
                        bool initialiseCells=true);

    /**
     * Add an update rule to be used in this simulation (use this to set the Hamiltonian).
     *
     * @param pUpdateRule shared pointer to a Potts update rule law
     */
    void AddUpdateRule(boost::shared_ptr<AbstractPottsUpdateRule<DIM> > pUpdateRule);

    /**
     * Overridden OutputAdditionalSimulationSetup() method.
     * Outputs the update rule information.
     * @return mOutputCellVelocities
     */
    bool GetOutputCellVelocities();

    /**
     * Set mOutputCellVelocities.
     *
     * @param outputCellVelocities the new value of mOutputCellVelocities
     */
    void SetOutputCellVelocities(bool outputCellVelocities);

    /**
     * Overridden SetupSolve() method to setup the cell velocities file.
     */
    virtual void SetupSolve();

    /**
     * Overridden AfterSolve() method to close the cell velocities file.
     */
    virtual void AfterSolve();

    /**
     * Overridden OutputAdditionalSimulationSetup method to output the force and cell
     * population boundary condition information.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputAdditionalSimulationSetup(out_stream& rParamsFile);

    /**
     * Outputs simulation parameters to file.
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
