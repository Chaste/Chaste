/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef VERTEXCRYPTSIMULATION2D_HPP_
#define VERTEXCRYPTSIMULATION2D_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "CellBasedSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"


/**
 * A 2D vertex-based crypt simulation object.
 */
class VertexCryptSimulation2d : public CellBasedSimulation<2>
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestVertexCryptSimulation2d;

private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<CellBasedSimulation<2> >(*this);
        archive & mUseJiggledBottomCells;
    }

    /** Helper member that is a static cast of the cell population. */
    VertexBasedCellPopulation<2>* mpStaticCastCellPopulation;

    /** Whether to use a flat bottom surface or to jiggle the cells on the bottom surface */
    bool mUseJiggledBottomCells;

    /**
     * Overridden WriteVisualizerSetupFile() method.
     *
     * Writes out special information about the mesh to the visualizer.
     */
    void WriteVisualizerSetupFile();

    /**
     * Overridden CalculateCellDivisionVector() method.
     *
     * By default this method returns the zero vector. If the parent cell
     * is a stem cell, then this method returns the vector (0,1). This is
     * then used by the VertexBasedCellPopulation method AddCell() as the axis along
     * which the cell divides.
     *
     * @param pParentCell the parent cell
     */
    c_vector<double, 2> CalculateCellDivisionVector(CellPtr pParentCell);

public :

    /**
     *  Constructor.
     *
     *  @param rCellPopulation A cell population object
     *  @param forceCollection The mechanics to use in the simulation
     *  @param deleteCellPopulationAndForceCollection Whether to delete the cell population and force collection on destruction to free up memory
     *  @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     */
    VertexCryptSimulation2d(AbstractCellPopulation<2>& rCellPopulation,
                      std::vector<AbstractForce<2>*> forceCollection,
                      bool deleteCellPopulationAndForceCollection=false,
                      bool initialiseCells=true);

    /**
     * Set method for mUseJiggledBottomCells.
     */
    void UseJiggledBottomCells();

    /**
     * Overridden ApplyCellPopulationBoundaryConditions() method.
     *
     * If an instance of WntConcentration is not set up, then stem cells at the
     * bottom of the crypt are pinned. Any cell that has moved below the bottom
     * of the crypt is moved back up.
     *
     * @param rOldLocations the node locations at the previous time step
     */
    void ApplyCellPopulationBoundaryConditions(const std::vector<c_vector<double,2> >& rOldLocations);

    /**
     * Outputs Simulation Parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(VertexCryptSimulation2d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VertexCryptSimulation2d.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const VertexCryptSimulation2d * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2> * p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
    const std::vector<AbstractForce<2>*> force_collection = t->rGetForceCollection();
    ar & force_collection;
}

/**
 * De-serialize constructor parameters and initialise a VertexCryptSimulation2d.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, VertexCryptSimulation2d * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar & p_cell_population;
    std::vector<AbstractForce<2>*> force_collection;
    ar & force_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)VertexCryptSimulation2d(*p_cell_population, force_collection, true, false);
}
}
} // namespace

#endif /*VERTEXCRYPTSIMULATION2D_HPP_*/
