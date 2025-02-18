/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef IMMERSEDBOUNDARYMORSEINTERACTIONFORCE_HPP_
#define IMMERSEDBOUNDARYMORSEINTERACTIONFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "Exception.hpp"

#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"

#include <iostream>

/**
 * A force class for use in immersed boundary simulations. This force implements
 * Morse-potential-like links between nodes in adjacent immersed boundaries
 * (https://en.wikipedia.org/wiki/Morse_potential). The well width is a constant
 * interaction strength, the rest length is an equilibrium bond distance, and
 * the well width is a parameter governing the profile of the curve.
 */
template <unsigned DIM>
class ImmersedBoundaryMorseInteractionForce : public AbstractImmersedBoundaryForce<DIM>
{
private:
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractImmersedBoundaryForce<DIM> >(*this);
        archive& mWellDepth;
        archive& mRestLength;
        archive& mLaminaWellDepthMult;
        archive& mLaminaRestLengthMult;
        archive& mWellWidth;
    }

    /**
     * The basic interaction strength. Initialised to 1e3 in constructor.
     */
    double mWellDepth;

    /**
     * The basic rest length associated with interactions, as a fraction of cell
     * population's interaction distance. Initialised to 0.25 in constructor.
     */
    double mRestLength;

    /**
     * Multiplicative factor to change strength of interactions involving a
     * lamina node. Initialised to 1.0 in constructor.
     */
    double mLaminaWellDepthMult;

    /**
     * Multiplicative factor to change equilibrium length of interactions
     * involving a lamina node. Initialised to 1.0 in constructor.
     */
    double mLaminaRestLengthMult;

    /**
     * The well width as a fraction of the cell population's interaction
     * distance. Initialised to 0.25 in constructor.
     */
    double mWellWidth;

public:
    /**
     * Constructor.
     */
    ImmersedBoundaryMorseInteractionForce();

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundaryMorseInteractionForce();

    /**
     * Overridden AddImmersedBoundaryForceContribution() method.
     *
     * Calculates the force on each node in the immersed boundary cell
     * population as a result of cell-cell interactions.
     *
     * @param rNodePairs reference to a vector set of node pairs between which
     *                   to contribute the force
     * @param rCellPopulation reference to the cell population
     */
    void AddImmersedBoundaryForceContribution(
        std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
        ImmersedBoundaryCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputImmersedBoundaryForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputImmersedBoundaryForceParameters(out_stream& rParamsFile);

    /** @return mWellDepth */
    double GetWellDepth() const;

    /**
     * Set mWellDepth.
     *
     * @param wellDepth the new value of mWellDepth
     */
    void SetWellDepth(double wellDepth);

    /** @return mRestLength */
    double GetRestLength() const;

    /**
     * Set mRestLength.
     *
     * @param restLength the new value of mRestLength
     */
    void SetRestLength(double restLength);

    /** @return mLaminaWellDepthMult */
    double GetLaminaWellDepthMult() const;

    /**
     * Set mLaminaWellDepthMult.
     *
     * @param laminaWellDepthMult the new value of mLaminaWellDepthMult
     */
    void SetLaminaWellDepthMult(double laminaWellDepthMult);

    /** @return mLaminaRestLengthMult */
    double GetLaminaRestLengthMult() const;

    /**Set mLaminaRestLengthMult.
     *
     * @param laminaRestLengthMult the new value of mLaminaRestLengthMult
     */
    void SetLaminaRestLengthMult(double laminaRestLengthMult);

    /** @return mWellWidth */
    double GetWellWidth() const;

    /**
     * Set mWellWidth.
     *
     * @param wellWidth the new value of mWellWidth
     */
    void SetWellWidth(double wellWidth);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryMorseInteractionForce)

#endif /*IMMERSEDBOUNDARYMORSEINTERACTIONFORCE_HPP_*/
