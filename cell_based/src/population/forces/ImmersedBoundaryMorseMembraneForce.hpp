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

#ifndef IMMERSEDBOUNDARYMORSEMEMBRANEFORCE_HPP_
#define IMMERSEDBOUNDARYMORSEMEMBRANEFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"

/**
 * A force class for use in immersed boundary simulations. This force implements
 * Morse-potential-like links between adjacent nodes in each immersed boundary
 * (https://en.wikipedia.org/wiki/Morse_potential). The well width is a constant
 * interaction strength, the rest length is an equilibrium bond distance, and
 * the well width is a parameter governing the profile of the curve.
 */
template <unsigned DIM>
class ImmersedBoundaryMorseMembraneForce : public AbstractImmersedBoundaryForce<DIM>
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
        archive& mElementWellDepth;
        archive& mElementRestLength;
        archive& mLaminaWellDepth;
        archive& mLaminaRestLength;
        archive& mWellWidth;
    }

    /**
     * The basic interaction strength. Initialised to 1e6 in constructor.
     */
    double mElementWellDepth;

    /**
     * The rest length associated with each element as a fraction of the average
     * node spacing. Initialised to 0.5 in constructor.
     */
    double mElementRestLength;

    /**
     * The spring constant associated with each lamina. Initialised to 1e6 in
     * constructor.
     */
    double mLaminaWellDepth;

    /**
     * The rest length associated with each lamina as a fraction of the average
     * node spacing. Initialised to 0.5 in constructor.
     */
    double mLaminaRestLength;

    /**
     * The well width as a fraction of the average node spacing in either an
     * element or lamina. Initialised to 0.25 in constructor.
     */
    double mWellWidth;

    /**
     * Helper method for AddImmersedBoundaryForceContribution.
     * Calculates forces, and can accept either an element or a lamina
     *
     * @tparam ELEMENT_DIM either DIM or DIM-1 depending on whether receiving an
     *     element or a lamina
     * @param rElement the element or lamina add forces to
     * @param rCellPopulation the immersed boundary cell population
     * @param intrinsicSpacingSquared the square of the intrinsic spacing
     */
    template <unsigned ELEMENT_DIM>
    void CalculateForcesOnElement(ImmersedBoundaryElement<ELEMENT_DIM, DIM>& rElement,
                                  ImmersedBoundaryCellPopulation<DIM>& rCellPopulation,
                                  double intrinsicSpacingSquared);

public:
    /** Constructor */
    ImmersedBoundaryMorseMembraneForce();

    /** Destructor */
    virtual ~ImmersedBoundaryMorseMembraneForce();

    /**
     * Overridden AddImmersedBoundaryForceContribution() method.
     * Calculates basic elasticity in the membrane of each immersed boundary as
     * a result of interactions.
     *
     * @param rNodePairs reference to a vector set of node pairs between which to
     *     contribute the force
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

    /** @return mElementWellDepth */
    double GetElementWellDepth() const;

    /**
     * Set mElementWellDepth.
     *
     * @param elementWellDepth the new value of mElementWellDepth
     */
    void SetElementWellDepth(double elementWellDepth);

    /** @return mElementRestLength */
    double GetElementRestLength() const;

    /**
     * Set mElementRestLength/
     *
     * @param elementRestLength the new value of mElementRestLength
     */
    void SetElementRestLength(double elementRestLength);

    /** @return mLaminaWellDepth */
    double GetLaminaWellDepth() const;

    /**
     * Set mLaminaWellDepth.
     *
     * @param laminaWellDepth the new value of mLaminaWellDepth
     */
    void SetLaminaWellDepth(double laminaWellDepth);

    /** @return mLaminaRestLength */
    double GetLaminaRestLength() const;

    /**
     * Set mLaminaRestLength/
     *
     * @param laminaRestLength the new value of mLaminaRestLength
     */
    void SetLaminaRestLength(double laminaRestLength);

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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryMorseMembraneForce)

#endif /*IMMERSEDBOUNDARYMORSEMEMBRANEFORCE_HPP_*/
