/*

Copyright (c) 2005-2015, University of Oxford.
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

/*
 * ImmersedBoundarySingleCellMigration.hpp
 *
 *  Created on: 20 Jun 2016
 *      Author: Bartosz Jan Bartmanski
 */

#ifndef IMMERSEDBOUNDARYSINGLECELLMIGRATIONFORCE_HPP_
#define IMMERSEDBOUNDARYSINGLECELLMIGRATIONFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"

/**
 * A force class for use in immersed boundary simulations. This force implements a force in a specific
 * direction on the specified element
 */
template <unsigned DIM>
class SingleCellMigrationForce : public AbstractImmersedBoundaryForce<DIM>
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
        archive& mStrength;
        archive& mElementIndex;
    }

    /** Strength of the force */
	double mStrength;

	/** direction of the force */
	c_vector<double, DIM> mDirection;

	/** Index element that will have this force applied to it */
	unsigned mElementIndex;

    /**
     * Helper method for AddImmersedBoundaryForceContribution.
     * Calculates forces, and can accept either an element or a lamina
     *
     * @tparam ELEMENT_DIM either DIM or DIM-1 depending on whether receiving an element or a lamina
     * @param rElement the element or lamina add forces to
     * @param rCellPopulation the immersed boundary cell population
     */
    template <unsigned ELEMENT_DIM>
    void CalculateForcesOnElement(ImmersedBoundaryElement<ELEMENT_DIM, DIM>& rElement,
                                  ImmersedBoundaryCellPopulation<DIM>& rCellPopulation,
                                  double intrinsicSpacingSquared);

public:
    /** Constructor */
    SingleCellMigrationForce();

    /** Destructor */
   virtual  ~SingleCellMigrationForce();

    /**
     * Overridden AddImmersedBoundaryForceContribution() method.
     * Calculates basic elasticity in the membrane of each immersed boundary as a result of interactions.
     *
     * @param rNodePairs reference to a vector set of node pairs between which to contribute the force
     * @param rCellPopulation reference to the cell population
     */
    void AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                              ImmersedBoundaryCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputImmersedBoundaryForceParameters() method.
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputImmersedBoundaryForceParameters(out_stream& rParamsFile);

    /** Setting the index of the element to which the force will be applied to */
    void SetElementIndex(unsigned element_index);

    /** Getting the index of the element to which the force will be applied to */
    unsigned GetElementIndex() const;

    /** Setting the strength of the force */
    void SetStrength(double strength);

    /** Getting the strength of the force */
    double GetStrength() const;

    /** Setting the direction of the force - will be normalized */
    void SetDirection(c_vector<double,DIM> direction);

    /** Getting the direction of the force */
    c_vector<double, DIM> GetDirection() const;



};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SingleCellMigrationForce)

#endif /* IMMERSEDBOUNDARYSINGLECELLMIGRATIONFORCE_HPP_ */
