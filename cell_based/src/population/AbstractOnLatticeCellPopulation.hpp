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

#ifndef ABSTRACTONLATTICECELLPOPULATION_HPP_
#define ABSTRACTONLATTICECELLPOPULATION_HPP_

#include "AbstractCellPopulation.hpp"
#include "AbstractUpdateRule.hpp"

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
        archive & mUpdateRuleCollection;
        archive & mUpdateNodesInRandomOrder;
        archive & mIterateRandomlyOverUpdateRuleCollection;
    }

protected:

    /**
     * The update rules used to determine the new location of the cells.
     * These rules specify how individual cells move into free spaces.
     */
    std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > > mUpdateRuleCollection;

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
     * Constructor that just takes in a mesh.
     *
     * @param rMesh the mesh for the cell population.
     */
    AbstractOnLatticeCellPopulation(AbstractMesh<DIM, DIM>& rMesh);

public:

    /**
     * Default constructor.
     *
     * @param rMesh a refernce to the mesh underlying the cell population
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     *            (defaults to false)
     */
    AbstractOnLatticeCellPopulation(AbstractMesh<DIM, DIM>& rMesh,
                                    std::vector<CellPtr>& rCells,
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
    virtual void UpdateCellLocations(double dt)=0;

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
     * Overridden SetNode() method.
     *
     * This method throws an exception if called on a subclass of AbstractOnLatticeCellPopulation,
     * since in such classes the lattice is assumed to be fixed.
     *
     * @param index the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    void SetNode(unsigned index, ChastePoint<DIM>& rNewLocation);

    /**
     * Overridden GetNeighbouringNodeIndices() method.
     *
     * This method currently returns an exception as the two types of neighbourhood
     * (Moore and Von Neumann) are defined in the PottsMesh.
     *
     * @param index the node index
     * @return the set of neighbouring node indices.
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned index);

    /**
     * Outputs parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellPopulationParameters(out_stream& rParamsFile);

    /**
     * Overridden GetDefaultTimeStep() method.
     *
     * @return a default value for the time step to use
     * when simulating the cell population.
     *
     * A hard-coded value of 0.1 is returned. However, note that the time
     * step can be reset by calling SetDt() on the simulation object used to
     * simulate the cell population.
     */
    virtual double GetDefaultTimeStep();

    /**
     * Add an update rule to be used with this population.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pUpdateRule pointer to an update rule
     */
    virtual void AddUpdateRule(boost::shared_ptr<AbstractUpdateRule<DIM> > pUpdateRule)=0;

    /**
     * Remove any update rules previously passed to this population.
     */
    virtual void RemoveAllUpdateRules();

    /**
     * Get the collection of update rules to be used with this population.
     *
     * @return the update rule collection
     */
    virtual const std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > > GetUpdateRuleCollection() const;
};

#endif /*ABSTRACTONLATTICECELLPOPULATION_HPP_*/
