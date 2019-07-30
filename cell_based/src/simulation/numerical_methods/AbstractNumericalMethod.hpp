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

#ifndef ABSTRACTNUMERICALMETHOD_HPP_
#define ABSTRACTNUMERICALMETHOD_HPP_

#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "Identifiable.hpp"

#include "AbstractOffLatticeCellPopulation.hpp"
#include "AbstractForce.hpp"

/**
 * An abstract class representing a numerical method for off lattice cell based simulations.
 *
 * Numerical methods have access to the cell population and the force collection.
 * The method is then responsible for evaluating forces at whatever times and positions
 * are required, then updating all node positions.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AbstractNumericalMethod : public Identifiable
{
    friend class TestNumericalMethods;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Save or restore the simulation.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mUseAdaptiveTimestep;
        archive & mUseUpdateNodeLocation;
        archive & mGhostNodeForcesEnabled;
    }

protected:

    /** Pointer to the cell population being updated by this method*/
    AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>* mpCellPopulation;

    /** Pointer to the force collection to apply*/
    std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >* mpForceCollection;

    /**
     * Whether the numerical method uses an adaptive time step.
     * Initialized to false in the AbstractNumericalMethod constructor.
     */
    bool mUseAdaptiveTimestep;

    /**
     * A boolean flag indicating whether non-forward Euler methods are supported for this type of cell
     * population. Allows us to fall back to the old method of updating node positions for populations
     * that require it (This is only for NodeBasedCellPopulationWithBuskeUpdates).
     *
     * Initialized to false in the AbstractNumericalMethod constructor.
     *
     * \todo #2087 Consider replacing this with static_casts
     */
    bool mUseUpdateNodeLocation;

    /**
     * A boolean indicating whether this cell population type contains ghost nodes.
     * Initialized to true in the AbstractNumericalMethod constructor.
     */
    bool mGhostNodeForcesEnabled;

    /**
     * Computes and returns the force on each node, including the damping factor
     * @return A vector of applied forces
     */
    std::vector<c_vector<double, SPACE_DIM> > ComputeForcesIncludingDamping();

    /**
     * Saves the current location of each cell in the population in a vector.
     * @return A vector of cell positions
     */
    std::vector<c_vector<double, SPACE_DIM> > SaveCurrentLocations();

    /**
     * Updates a single node's position, taking into account periodic boundary conditions
     *
     * @param nodeIndex Index of the node to update
     * @param newPosition C_vector holding the new node position
     */
     void SafeNodePositionUpdate(unsigned nodeIndex, c_vector<double, SPACE_DIM> newPosition);

    /**
     * Detects whether a node has exceeded the acceptable displacement for one timestep.
     * If a step size exception has occurred, it either causes the simulation to terminate or,
     * in adaptive simulations, the exception is caught and the step size is reduced in response.
     *
     * @param nodeIndex Index of the node being examined
     * @param displacement Displacement of the node this step
     * @param dt Time step size
     */
    void DetectStepSizeExceptions(unsigned nodeIndex, c_vector<double,SPACE_DIM>& displacement, double dt);

public:

    /**
     * Constructor.
     * No input parameters are required, allowing the numerical method to be created first, then passed to
     * the simulation. The cell population and force collection pointers are then set by the simulation in
     * its constructor.
     */
    AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>();

    /**
     * Destructor.
     */
    virtual ~AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>();

    /**
     * Sets the pointer to the cell population updated by this method
     *
     * @param pPopulation Pointer to an AbstractOffLattice cell population
     */
    void SetCellPopulation(AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>* pPopulation);

    /**
     * Sets the pointer to the force collection applied by this method
     *
     * @param pForces Pointer to the simulation's force collection
     */
    void SetForceCollection(std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >* pForces);

    /**
     * Set mUseAdaptiveTimestep.
     *
     * @param useAdaptiveTimestep whether to use an adaptive time step
     */
    void SetUseAdaptiveTimestep(bool useAdaptiveTimestep);

    /**
     * Set mUseUpdateNodeLocation.
     *
     * @param useUpdateNodeLocation whether to use the population method UpdateNodeLocations()
     */
    void SetUseUpdateNodeLocation(bool useUpdateNodeLocation);

    /**
     * Get mUseUpdateNodeLocation.
     *
     * @return whether the population method UpdateNodeLocations() is being used
     */
    bool GetUseUpdateNodeLocation();

    /**
     * @return whether the numerical method uses an adaptive time step.
     */
    bool HasAdaptiveTimestep();

    /**
     * Updates node positions according to Newton's 2nd law with overdamping.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param dt Time step size
     */
    virtual void UpdateAllNodePositions(double dt)=0;

    /**
     * Saves the name of the numerical method to the parameters file
     *
     * @param rParamsFile Reference to the parameter output filestream
     */
    void OutputNumericalMethodInfo(out_stream& rParamsFile);

    /**
     * Saves any additional numerical method details to the parameters file.
     *
     * @param rParamsFile Reference to the parameter output filestream
     */
    virtual void OutputNumericalMethodParameters(out_stream& rParamsFile);
};

TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED(AbstractNumericalMethod)

#endif /*ABSTRACTNUMERICALMETHOD_HPP_*/
