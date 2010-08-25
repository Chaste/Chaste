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

#ifndef LATTICEBASEDTISSUESIMULATION_HPP_
#define LATTICEBASEDTISSUESIMULATION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "TissueSimulation.hpp"
#include "LatticeBasedTissue.hpp"
#include "AbstractUpdateRule.hpp"


/**
 * A lattice-based tissue simulation object.
 */
 template<unsigned DIM>
class LatticeBasedTissueSimulation : public TissueSimulation<DIM>
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestLatticeBasedTissueSimulation;

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
        archive & boost::serialization::base_object<TissueSimulation<DIM> >(*this);

        archive & mUpdateRuleCollection;
        archive & mIterateRandomlyOverUpdateRuleCollection;
        archive & mIterateRandomlyOverCells;
    }

    /** Helper member that is a static cast of the tissue. */
    LatticeBasedTissue<DIM>* mpStaticCastTissue;

    /** The rules used to determine the new locations of cells. */
    std::vector<AbstractUpdateRule<DIM>*> mUpdateRuleCollection;

    /** Whether delete the collection of update rules in the destructor. */
    bool mAllocatedMemoryForUpdateRuleCollection;

    /** Whether to iterate randomly over mUpdateRuleCollection when updating cell locations. */
    bool mIterateRandomlyOverUpdateRuleCollection;

    /** Whether to iterate randomly over cells when updating cell locations. */
    bool mIterateRandomlyOverCells;

    /**
     * Overridden CalculateCellDivisionVector() method.
     *
     * By default this method returns the zero vector. If the parent cell
     * is a stem cell, then this method returns the vector (0,1). This is
     * then used by the VertexBasedTissue method AddCell() as the axis along
     * which the cell divides.
     *
     * @param pParentCell the parent cell
     */
    c_vector<double, DIM> CalculateCellDivisionVector(TissueCellPtr pParentCell);

    /**
     * Move each cell to a new lattice site for this timestep by calling the
     * tissue method MoveCell().
     */
    void UpdateCellLocations();

public:

    /**
     * Constructor.
     *
     * @param rTissue A tissue facade class (contains a mesh and cells)
     * @param updateRuleCollection The cell movement rules to use in the simulation
     * @param iterateRandomlyOverUpdateRuleCollection whether to iterate randomly over 
     *        mUpdateRuleCollection when updating cell locations (defaults to false)
     * @param iterateRandomlyOverCells whether to iterate randomly over cells when 
     *        updating cell locations (defaults to false)
     * @param deleteTissueAndForceCollection Whether to delete the tissue and force 
     *        collection on destruction to free up memory (defaults to false)
     * @param initialiseCells whether to initialise cells (defaults to true, set to 
     *        false when loading from an archive)
     */
    LatticeBasedTissueSimulation(AbstractTissue<DIM>& rTissue,
                      std::vector<AbstractUpdateRule<DIM>*> updateRuleCollection,
                      bool iterateRandomlyOverUpdateRuleCollection=false,
                      bool iterateRandomlyOverCells=false,
                      bool deleteTissueAndForceCollection=false,
                      bool initialiseCells=true);

    /**
     * Destructor.
     */
    ~LatticeBasedTissueSimulation();

    /**
     * Main solve method.
     *
     * This method sets up the simulation time, creates output files, and initialises the
     * tissue. It then iterates through a time loop. At each time step, first any cell death
     * or birth is implemented, then the tissue topology is updated, then the update rule is
     * recalculated and the tissue evolved according to whatever update rules are present in
     * the simulation, and finally the results for that time step are output to file. At the
     * end of the time loop, the method closes any output files.
     */
    void Solve();

    /**
     * @return const reference to mUpdateRuleCollection (used in archiving).
     */
    const std::vector<AbstractUpdateRule<DIM>*> rGetUpdateRuleCollection() const;

    /**
     * Outputs Simulation Parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses. Currently not used for Lattice Based Simulations.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);

};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LatticeBasedTissueSimulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a LatticeBasedTissueSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const LatticeBasedTissueSimulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    const std::vector<AbstractUpdateRule<DIM>*> update_rule_collection = t->rGetUpdateRuleCollection();
    ar & update_rule_collection;
}

/**
 * De-serialize constructor parameters and initialise a LatticeBasedTissueSimulation.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, LatticeBasedTissueSimulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;
    ar >> p_tissue;
    std::vector<AbstractUpdateRule<DIM>*> update_rule_collection;
    ar >> update_rule_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)LatticeBasedTissueSimulation<DIM>(*p_tissue, update_rule_collection, false, false, true, false);
}
}
} // namespace

#endif /*LATTICEBASEDTISSUESIMULATION_HPP_*/
