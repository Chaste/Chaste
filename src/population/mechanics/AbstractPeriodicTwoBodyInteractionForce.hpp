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
#ifndef ABSTRACTPERIODICTWOBODYINTERACTIONFORCE_HPP_
#define ABSTRACTPERIODICTWOBODYINTERACTIONFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"

/**
 * An abstract class for two-body force laws.
 */
template<unsigned DIM>
class AbstractPeriodicTwoBodyInteractionForce : public AbstractTwoBodyInteractionForce<DIM>
{
private:

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
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mInitialWidth;
    }

protected:

    /** Initial width of the mesh for use in periodic boundaries. */
    double mInitialWidth;
    
    /** An extended mesh, used to implement periodicity. */
    MutableMesh<DIM,DIM>* mpExtendedMesh;

    /** A map from node indices in mpExtendedMesh to node indices in the cell population. */
    std::map<unsigned, unsigned> mExtendedMeshNodeIndexMap;

public:

    /**
     * Constructor.
     */
    AbstractPeriodicTwoBodyInteractionForce();

    /**
     * Destructor.
     */
    ~AbstractPeriodicTwoBodyInteractionForce();

    /**
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by rCalculateVelocitiesOfEachNode()
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     *
     * @return The force exerted on Node A by Node B.
     */
    virtual c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM>& rCellPopulation)=0;

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Returns the initial width.
     */
    double GetInitialWidth();

    /**
     * Returns the initial width.
     *
     * @param initialWidth the initial width
     */
    void SetInitialWidth(double initialWidth);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#endif /* ABSTRACTPERIODICTWOBODYINTERACTIONFORCE_HPP_ */
