/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef ABSTRACTOFFLATTICECELLPOPULATION_HPP_
#define ABSTRACTOFFLATTICECELLPOPULATION_HPP_

#include "AbstractCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"

/**
 * An abstract facade class encapsulating an off-lattice (centre- or
 * vertex-based) cell population.
 */
template<unsigned DIM>
class AbstractOffLatticeCellPopulation : public AbstractCellPopulation<DIM>
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
        archive & mDampingConstantNormal;
        archive & mDampingConstantMutant;
    }

protected:

    /**
     * Damping constant for normal cells has units of kg s^-1
     * Represented by the parameter eta in the model by Meineke et al (2001) in
     * their off-lattice model of the intestinal crypt
     * (doi:10.1046/j.0960-7722.2001.00216.x).
     */
    double mDampingConstantNormal;

    /**
     * Damping constant for mutant cells has units of kg s^-1.
     */
    double mDampingConstantMutant;

    /**
     * Constructor for use by archiving only.
     */
    AbstractOffLatticeCellPopulation();

public:

    /**
     * Default constructor.
     *
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    AbstractOffLatticeCellPopulation(std::vector<CellPtr>& rCells,
                                    const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Add a new node to the cell population.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in cell population.
     */
    virtual unsigned AddNode(Node<DIM>* pNewNode)=0;

    /**
     * Move the node with a given index to a new point in space.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param nodeIndex the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    virtual void SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)=0;

    /**
     * Update the location of each node in the cell population given
     * a vector of forces on nodes and a time step over which to
     * integrate the equations of motion.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rNodeForces  forces on nodes
     * @param dt time step
     */
    virtual void UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)=0;

    /**
     * Get the damping constant for this node - ie d in drdt = F/d.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param nodeIndex the global index of this node
     * @return the damping constant at the node.
     */
    virtual double GetDampingConstant(unsigned nodeIndex)=0;

    /**
      * Set mDampingConstantNormal.
      *
      * @param dampingConstantNormal  the new value of mDampingConstantNormal
     */
    void SetDampingConstantNormal(double dampingConstantNormal);

    /**
     * Set mDampingConstantMutant.
     *
     * @param dampingConstantMutant  the new value of mDampingConstantMutant
     */
    void SetDampingConstantMutant(double dampingConstantMutant);

    /**
     * @return mDampingConstantNormal
     */
    double GetDampingConstantNormal();

    /**
     * @return mDampingConstantMutant
     */
    double GetDampingConstantMutant();

    /**
     * Overridden OutputCellPopulationParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellPopulationParameters(out_stream& rParamsFile);
};

#endif /*ABSTRACTOFFLATTICECELLPOPULATION_HPP_*/
