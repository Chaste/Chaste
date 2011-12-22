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

#ifndef ABSTRACTTWOBODYINTERACTIONFORCE_HPP_
#define ABSTRACTTWOBODYINTERACTIONFORCE_HPP_

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
/**
 * An abstract class for two-body force laws.
 */
template<unsigned DIM>
class AbstractTwoBodyInteractionForce : public AbstractForce<DIM>
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
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mUseCutOffLength;
        archive & mMechanicsCutOffLength;
    }

protected:

    /** Whether to have zero force if the cells are far enough apart. */
    bool mUseCutOffLength;

    /** Mechanics cut off length. */
    double mMechanicsCutOffLength;

public:

    /**
     * Constructor.
     */
    AbstractTwoBodyInteractionForce();

    /**
     * Whether the force is using a cut of length.
     *
     * @return mUseCutOffLength
     */
    bool GetUseCutOffLength();

    /**
     * Use a cutoff point, ie specify zero force if two cells are greater
     * than the cutoff distance apart.
     *
     * @param cutOffLength the cutoff to use
     */
    void SetCutOffLength(double cutOffLength);

    /**
     * Get the cutoff length.
     *
     * @return mCutOffLenght
     */
    double GetCutOffLength();

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
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#endif /*ABSTRACTTWOBODYINTERACTIONFORCE_HPP_*/
