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

#ifndef ABSTRACTONLATTICECELLPOPULATION_HPP_
#define ABSTRACTONLATTICECELLPOPULATION_HPP_

#include "AbstractCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"

/**
 * An abstract class for on-lattice cell populations.
 */
template<unsigned DIM>
class AbstractOnLatticeCellPopulation : public AbstractCellPopulation<DIM>
{
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
        archive & boost::serialization::base_object<AbstractCellPopulation<DIM> >(*this);
        archive & mUpdateNodesInRandomOrder;
        archive & mIterateRandomlyOverUpdateRuleCollection;
    }

protected:

    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this cell population has been de-serialized.
     */
    bool mDeleteMesh;

    /**
     * Whether to update nodes in random order.
     * Initialized to true in the constructor.
     */
    bool mUpdateNodesInRandomOrder;

    /**
     * Whether to iterate randomly over mUpdateRuleCollection when updating cell locations.
     * Initialized to false in the constructor.
     */
    bool mIterateRandomlyOverUpdateRuleCollection;

    /**
     * Constructor for use by archiving only.
     */
    AbstractOnLatticeCellPopulation();

public:

    /**
     * Default constructor.
     *
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     *            (defaults to false)
     */
    AbstractOnLatticeCellPopulation(std::vector<CellPtr>& rCells,
                                    const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                                    bool deleteMesh=false);

    /**
     * Destructor.
     */
    virtual ~AbstractOnLatticeCellPopulation();

    /**
     * Update cell locations over the course of a time step of specified length.
     * 
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param dt time step
     */
    virtual void UpdateCellLocations(double dt)
    {}

    /**
     * Get whether we update nodes in a random order.
     * 
     * @return mUpdateNodesInRandomOrder
     */
    bool GetUpdateNodesInRandomOrder();
    
    /**
     * Get whether we update nodes in a random order.
     * 
     * @param updateNodesInRandomOrder Whether to update nodes in a random order.
     */
    void SetUpdateNodesInRandomOrder(bool updateNodesInRandomOrder);

    /**
     * Set mIterateRandomlyOverUpdateRuleCollection.
     * 
     * @param iterateRandomly whether to iterate randomly over mUpdateRuleCollection
     */
    void SetIterateRandomlyOverUpdateRuleCollection(bool iterateRandomly);

    /**
     * @return mIterateRandomlyOverUpdateRuleCollection.
     */
    bool GetIterateRandomlyOverUpdateRuleCollection();

    /**
     * Outputs CellPopulation parameters to file
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellPopulationParameters(out_stream& rParamsFile);
};

#endif /*ABSTRACTONLATTICECELLPOPULATION_HPP_*/
