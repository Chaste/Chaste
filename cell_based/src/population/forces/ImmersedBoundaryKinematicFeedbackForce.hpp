/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef IMMERSEDBOUNDARYKINEMATICFEEDBACKFORCE_HPP_
#define IMMERSEDBOUNDARYKINEMATICFEEDBACKFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"

#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"

#include <iostream>

/**
 * A force class for use in immersed boundary simulations. This force implements kinematic feedback; amplification of
 * relative motion of cells sliding past one another.
 */
template<unsigned DIM>
class ImmersedBoundaryKinematicFeedbackForce : public AbstractImmersedBoundaryForce<DIM>
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
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractImmersedBoundaryForce<DIM> >(*this);
        archive & mSpringConst;
        archive & mRestLength;
    }

    /** The basic spring constant associated with interactions */
    double mSpringConst;

    /** The basic rest length associated with interactions, as a fraction of cell population's interaction distance */
    double mRestLength;

    /** Vector to contain the location of each node at the previous time step */
    std::vector<c_vector<double, DIM>> mPreviousLocations;

    /**
     * Helper function for AddImmersedBoundaryForceContribution().  Repopulate mPreviousLocations with new values.
     * @param rCellPopulation the cell population
     */
    void UpdatePreviousLocations(ImmersedBoundaryCellPopulation<DIM>& rCellPopulation);

public:

    /**
     * Constructor.
     */
    ImmersedBoundaryKinematicFeedbackForce();

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundaryKinematicFeedbackForce() = default;

    /**
     * Overridden AddImmersedBoundaryForceContribution() method.
     * Calculates the force on each node in the immersed boundary cell population as a result of kinematic feedback.
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

    /** @return mSpringConst */
    double GetSpringConst() const;

    /** @param springConst the new value of mSpringConst */
    void SetSpringConst(double springConst);

    /** @return mRestLength */
    double GetRestLength() const;

    /** @param restLength the new value of mRestLength */
    void SetRestLength(double restLength);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryKinematicFeedbackForce)

#endif /*IMMERSEDBOUNDARYKINEMATICFEEDBACKFORCE_HPP_*/
