/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef ABSTRACTVARIABLESIZETWOBODYINTERACTIONFORCE_HPP_
#define ABSTRACTVARIABLESIZETWOBODYINTERACTIONFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * An abstract two-body interaction force law for variable-size cells.
 *
 * This class contains common mechanics used by spring-like interactions where
 * the rest length can vary with node radii, spring marking after cell division,
 * and apoptosis. Concrete subclasses define the final force law via
 * CalculateLinkInteraction().
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AbstractVariableSizeTwoBodyInteractionForce : public AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM>
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
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mSpringStiffness;
        archive & mDivisionRestingSpringLength;
        archive & mSpringGrowthDuration;
    }

protected:

    /** Spring stiffness. */
    double mSpringStiffness;

    /** Initial resting spring length after cell division. */
    double mDivisionRestingSpringLength;

    /** Time taken for spring rest length to regrow after division. */
    double mSpringGrowthDuration;

public:

    /**
     * Constructor.
     */
    AbstractVariableSizeTwoBodyInteractionForce();

    /**
     * Destructor.
     */
    virtual ~AbstractVariableSizeTwoBodyInteractionForce();

    /**
     * Return a multiplication factor for the spring constant.
     *
     * Subclasses may override this method.
     * 
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @param isCloserThanRestLength whether the nodes are currently closer than the rest length
     * @return the multiplication factor for the spring constant
     * 
     */
    virtual double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                              unsigned nodeBGlobalIndex,
                                                              AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                              bool isCloserThanRestLength);

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * This computes the geometry/rest length logic shared by variable-size
     * interactions and then delegates the final force-law expression to
     * CalculateLinkInteraction().
     * 
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @return the force exerted on Node A by Node B    
     */
    c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                           unsigned nodeBGlobalIndex,
                                                           AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Calculate the force between two linked nodes after overlap/rest length
     * have been computed.
     * 
     * @param overlap the amount by which the distance between nodes is less than the rest length
     * @param restLength the rest length of the spring between the nodes
     * @param rUnitDifference the unit vector pointing from one node to the other
     * @param multiplicationFactor a multiplication factor for the spring constant
     * @return the force vector between the two nodes
     * 
     * As this method is pure virtual, it must be overridden in subclasses.
     */
    virtual c_vector<double, SPACE_DIM> CalculateLinkInteraction(double overlap,
                                                                 double restLength,
                                                                 const c_vector<double, SPACE_DIM>& rUnitDifference,
                                                                 double multiplicationFactor)=0;

    /**
     * @return mSpringStiffness
     */
    double GetSpringStiffness();

    /**
     * @return mDivisionRestingSpringLength
     */
    double GetDivisionRestingSpringLength();

    /**
     * @return mSpringGrowthDuration
     */
    double GetSpringGrowthDuration();

    /**
     * Set mSpringStiffness.
     * 
     * @param springStiffness the new value of mSpringStiffness
     */
    void SetSpringStiffness(double springStiffness);

    /**
     * Set mDivisionRestingSpringLength.
     * 
     * @param divisionRestingSpringLength the new value of mDivisionRestingSpringLength
     */
    void SetDivisionRestingSpringLength(double divisionRestingSpringLength);

    /**
     * Set mSpringGrowthDuration.
     * 
     * @param springGrowthDuration the new value of mSpringGrowthDuration   
     */
    void SetSpringGrowthDuration(double springGrowthDuration);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#endif /*ABSTRACTVARIABLESIZETWOBODYINTERACTIONFORCE_HPP_*/