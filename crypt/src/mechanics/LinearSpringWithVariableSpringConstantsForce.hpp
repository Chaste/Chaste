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

#ifndef LINEARSPRINGWITHVARIABLESPRINGCONSTANTSFORCE_HPP
#define LINEARSPRINGWITHVARIABLESPRINGCONSTANTSFORCE_HPP

#include "GeneralisedLinearSpringForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

// Needed here to avoid serialization errors (on Boost<1.37)
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

/**
 * A subclass of GeneralisedLinearSpringForce with variable spring constants.
 */
template<unsigned DIM>
class LinearSpringWithVariableSpringConstantsForce : public GeneralisedLinearSpringForce<DIM>
{
    friend class TestLinearSpringWithVariableSpringConstantsForce;

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
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<DIM> >(*this);
        archive & mUseEdgeBasedSpringConstant;
        archive & mUseMutantSprings;
        archive & mMutantMutantMultiplier;
        archive & mNormalMutantMultiplier;
        archive & mUseBCatSprings;
        archive & mUseApoptoticSprings;
        archive & mBetaCatSpringScaler;
        archive & mApoptoticSpringTensionStiffness;
        archive & mApoptoticSpringCompressionStiffness;
    }

protected:

    /** Whether to use spring constant proportional to cell-cell contact length/area (defaults to false). */
    bool mUseEdgeBasedSpringConstant;

    /** Whether to use different stiffnesses depending on whether either cell is a mutant. */
    bool mUseMutantSprings;

    /** Multiplier for spring stiffness if mutant. */
    double mMutantMutantMultiplier;

    /** Multiplier for spring stiffness if mutant. */
    double mNormalMutantMultiplier;

    /** Use springs which are dependent on beta-catenin levels. */
    bool mUseBCatSprings;

    /** Use springs which are dependent on whether cells are apoptotic. */
    bool mUseApoptoticSprings;

    /** Scaling factor for beta catenin to spring strength. */
    double mBetaCatSpringScaler;

    /** Non-dimensionalized 'stiffness' of a apoptotic cell under tension. */
    double mApoptoticSpringTensionStiffness;

    /** Non-dimensionalized 'stiffness' of a apoptotic cell under compression. */
    double mApoptoticSpringCompressionStiffness;

public:

    /**
     * Constructor.
     */
    LinearSpringWithVariableSpringConstantsForce();

    /**
     * Destructor.
     */
    ~LinearSpringWithVariableSpringConstantsForce();

    /**
     * Set whether to use an edge-based spring constant.
     *
     * @param useEdgeBasedSpringConstant whether to use an edge-based spring constant
     */
    void SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant);

    /**
     * Use different spring strengths depending on two cells:
     * Normal-normal, Normal-mutant, mutant-mutant
     *
     * @param useMutantSprings whether to use mutant springs
     * @param mutantMutantMultiplier the multiplier for springs connecting two mutant cells
     * @param normalMutantMultiplier the multiplier for springs connecting a mutant cell with a normal cell
     */
    void SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier=2, double normalMutantMultiplier=1.5);

    /**
     * Use the amount of beta-catenin on an edge to find spring constant.
     *
     * @param useBCatSprings whether to use beta-catenin-dependent spring stiffness
     */
    void SetBetaCateninSprings(bool useBCatSprings);

    /**
     * Set spring stiffness to be dependent on whether cells are apoptotic
     *
     * @param useApoptoticSprings whether to have apoptosis-dependent spring stiffness
     */
    void SetApoptoticSprings(bool useApoptoticSprings);

    /**
     * Return a multiplication factor for the spring constant, which
     * may depend on whether the given pair of neighbouring cells are
     * e.g. undergoing apoptosis, have mutations, or experience variable
     * levels of beta catenin.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
     *
     * @return the multiplication factor.
     */
    double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                      unsigned nodeBGlobalIndex,
                                                      AbstractCellPopulation<DIM>& rCellPopulation,
                                                      bool isCloserThanRestLength);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * @return mBetaCatSpringScaler
     */
    double GetBetaCatSpringScaler();

    /**
     * Set mBetaCatSpringScaler.
     *
     * @param betaCatSpringScaler the new value of mBetaCatSpringScaler
     */
    void SetBetaCatSpringScaler(double betaCatSpringScaler);

    /**
     * @return mApoptoticSpringTensionStiffness
     */
    double GetApoptoticSpringTensionStiffness();

    /**
     * Set mApoptoticSpringTensionStiffness.
     *
     * @param apoptoticSpringTensionStiffness the new value of mApoptoticSpringTensionStiffness
     */
    void SetApoptoticSpringTensionStiffness(double apoptoticSpringTensionStiffness);

    /**
     * @return mApoptoticSpringCompressionStiffness
     */
    double GetApoptoticSpringCompressionStiffness();

    /**
     * Set mApoptoticSpringCompressionStiffness.
     *
     * @param apoptoticSpringCompressionStiffness the new value of mApoptoticSpringCompressionStiffness
     */
    void SetApoptoticSpringCompressionStiffness(double apoptoticSpringCompressionStiffness);

    /**
     * Outputs force parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(LinearSpringWithVariableSpringConstantsForce)


#endif /*LINEARSPRINGWITHVARIABLESPRINGCONSTANTSFORCE_HPP*/
