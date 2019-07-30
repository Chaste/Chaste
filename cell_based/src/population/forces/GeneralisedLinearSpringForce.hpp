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

#ifndef GENERALISEDLINEARSPRINGFORCE_HPP_
#define GENERALISEDLINEARSPRINGFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A force law employed by Meineke et al (2001) in their off-lattice
 * model of the intestinal crypt (doi:10.1046/j.0960-7722.2001.00216.x).
 *
 * Each pair of neighbouring nodes are assumed to be connected by a linear
 * spring. The force of node \f$i\f$ is given
 * by
 *
 * \f[
 * \mathbf{F}_{i}(t) = \sum_{j} \mu_{i,j} ( || \mathbf{r}_{i,j} || - s_{i,j}(t) ) \hat{\mathbf{r}}_{i,j}.
 * \f]
 *
 * Here \f$\mu_{i,j}\f$ is the spring constant for the spring between nodes
 * \f$i\f$ and \f$j\f$, \f$s_{i,j}(t)\f$ is its natural length at time \f$t\f$,
 * \f$\mathbf{r}_{i,j}\f$ is their relative displacement and a hat (\f$\hat{}\f$)
 * denotes a unit vector.
 *
 * Length is scaled by natural length.
 * Time is in hours.
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class GeneralisedLinearSpringForce : public AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForces;

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
        archive & mMeinekeSpringStiffness;
        archive & mMeinekeDivisionRestingSpringLength;
        archive & mMeinekeSpringGrowthDuration;
    }

protected:

    /**
     * Spring stiffness.
     *
     * Represented by the parameter mu in the model by Meineke et al (2001) in
     * their off-lattice model of the intestinal crypt
     * (doi:10.1046/j.0960-7722.2001.00216.x).
     */
    double mMeinekeSpringStiffness;

    /**
     * Initial resting spring length after cell division.
     * Has units of cell size at equilibrium rest length
     *
     * The value of this parameter should be larger than mDivisionSeparation,
     * because of pressure from neighbouring springs.
     */
    double mMeinekeDivisionRestingSpringLength;

    /**
     * The time it takes for the springs rest length to increase from
     * mMeinekeDivisionRestingSpringLength to its natural length.
     *
     * The value of this parameter is usually the same as the M Phase of the cell cycle and defaults to 1.
     */
    double mMeinekeSpringGrowthDuration;

public:

    /**
     * Constructor.
     */
    GeneralisedLinearSpringForce();

    /**
     * Destructor.
     */
    virtual ~GeneralisedLinearSpringForce();

    /**
     * Return a multiplication factor for the spring constant, which
     * returns a default value of 1.
     *
     * This method may be overridden in subclasses.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
     *
     * @return the multiplication factor.
     */
    virtual double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                              unsigned nodeBGlobalIndex,
                                                              AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                              bool isCloserThanRestLength);

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by AddForceContribution()
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                     unsigned nodeBGlobalIndex,
                                                     AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);
    /**
     * @return mMeinekeSpringStiffness
     */
    double GetMeinekeSpringStiffness();

    /**
     * @return mMeinekeDivisionRestingSpringLength
     */
    double GetMeinekeDivisionRestingSpringLength();

    /**
     * @return mMeinekeSpringGrowthDuration
     */
    double GetMeinekeSpringGrowthDuration();

    /**
     * Set mMeinekeSpringStiffness.
     *
     * @param springStiffness the new value of mMeinekeSpringStiffness
     */
    void SetMeinekeSpringStiffness(double springStiffness);

    /**
     * Set mMeinekeDivisionRestingSpringLength.
     *
     * @param divisionRestingSpringLength the new value of mMeinekeDivisionRestingSpringLength
     */
    void SetMeinekeDivisionRestingSpringLength(double divisionRestingSpringLength);

    /**
     * Set mMeinekeSpringGrowthDuration.
     *
     * @param springGrowthDuration the new value of mMeinekeSpringGrowthDuration
     */
    void SetMeinekeSpringGrowthDuration(double springGrowthDuration);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(GeneralisedLinearSpringForce)

#endif /*GENERALISEDLINEARSPRINGFORCE_HPP_*/
