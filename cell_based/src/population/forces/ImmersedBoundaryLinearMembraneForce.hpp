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

#ifndef IMMERSEDBOUNDARYLINEARMEMBRANEFORCE_HPP_
#define IMMERSEDBOUNDARYLINEARMEMBRANEFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"

/**
 * A force class for use in immersed boundary simulations. This force implements
 * elastic links between adjacent nodes in each immersed boundary.
 */
template <unsigned DIM>
class ImmersedBoundaryLinearMembraneForce : public AbstractImmersedBoundaryForce<DIM>
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
        archive& mElementSpringConst;
        archive& mElementRestLength;
        archive& mLaminaSpringConst;
        archive& mLaminaRestLength;
    }

    /**
     * The spring constant associated with each element. Initialised to 1e6 in
     * constructor.
     */
    double mElementSpringConst;

    /**
     * The rest length associated with each element as a fraction of the average
     * node spacing. Initialised to 0.5 in constructor.
     */
    double mElementRestLength;

    /**
     * The spring constant associated with each lamina. Initialised to 1e6 in
     * constructor.
     */
    double mLaminaSpringConst;

    /**
     * The rest length associated with each lamina as a fraction of the average
     * node spacing. Initialised to 0.5 in constructor.
     */
    double mLaminaRestLength;

    /**
     * Helper method for AddImmersedBoundaryForceContribution.
     * Calculates forces, and can accept either an element or a lamina
     *
     * @tparam ELEMENT_DIM either DIM or DIM-1 depending on whether receiving an
     *         element or a lamina
     * @param rElement the element or lamina add forces to
     * @param rCellPopulation the immersed boundary cell population
     * @param intrinsicSpacingSquared the intrinsic node spacing
     */
    template <unsigned ELEMENT_DIM>
    void CalculateForcesOnElement(
        ImmersedBoundaryElement<ELEMENT_DIM, DIM>& rElement,
        ImmersedBoundaryCellPopulation<DIM>& rCellPopulation,
        double intrinsicSpacingSquared);

public:
    /** Constructor */
    ImmersedBoundaryLinearMembraneForce();

    /** Destructor */
    virtual ~ImmersedBoundaryLinearMembraneForce();

    /**
     * Overridden AddImmersedBoundaryForceContribution() method.
     * Calculates basic elasticity in the membrane of each immersed boundary as
     * a result of interactions.
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
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputImmersedBoundaryForceParameters(out_stream& rParamsFile);

    /** @return mElementSpringConst */
    double GetElementSpringConst() const;

    /**
     * Set mElementSpringConst.
     *
     * @param elementSpringConst the new value of mElementSpringConst
     */
    void SetElementSpringConst(double elementSpringConst);

    /** @return mElementRestLength */
    double GetElementRestLength() const;

    /**
     * Set mElementRestLength
     *
     * @param elementRestLength the new value of mElementRestLength
     */
    void SetElementRestLength(double elementRestLength);

    /** @return mLaminaSpringConst */
    double GetLaminaSpringConst() const;

    /**
     * Set mLaminaSpringConst.
     *
     * @param laminaSpringConst the new value of mLaminaSpringConst
     */
    void SetLaminaSpringConst(double laminaSpringConst);

    /** @return mLaminaRestLength */
    double GetLaminaRestLength() const;

    /**
     * Set mLaminaRestLength.
     *
     * @param laminaRestLength the new value of mLaminaRestLength
     */
    void SetLaminaRestLength(double laminaRestLength);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryLinearMembraneForce)

#endif /*IMMERSEDBOUNDARYLINEARMEMBRANEFORCE_HPP_*/
