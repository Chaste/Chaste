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

#ifndef POPULATIONTESTINGFORCE_HPP_
#define POPULATIONTESTINGFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractForce.hpp"

/**
 * A simple force law used to test node location updates across the off lattice population test files.
 */

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class PopulationTestingForce : public AbstractForce<ELEMENT_DIM, SPACE_DIM>
{
private:

    /**
     * Archiving.
     */
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
        archive & boost::serialization::base_object<AbstractForce<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mWithPositionDependence;
    }

    /**
     * Whether the applied force should depend on the position of the node, or only its index
     * Defaults to true, since that's better for testing a variety of numerical methods.
     * Crypt tests require a position independent test force though.
     */
    bool mWithPositionDependence;

public:

    /**
     * Constructor.
     *
     * @param hasPositionDependence specifys which testing force to use (defaults to true)
     */
    PopulationTestingForce(bool hasPositionDependence = true);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     *
     * THis has two forces depending on the member variable mWithPositionDependence
     *
     * If true then
     *
     * For node i
     * F[j] = 0.01 (j+1) i r[j]
     *
     * Where r is the position of node i.
     *
     * if false then
     *
     * F[j] = 0.01 (j+1) i
     *
     */
    void AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Helper method to return the expected step location for ForwardEulerNumericalMethod.
     *
     * @return the expected location after one step
     *
     * @param nodeIndex the index of the node
     * @param damping the damping constant
     * @param oldLocation the old location of the node
     * @param dt the step size
     */
    c_vector<double, SPACE_DIM> GetExpectedOneStepLocationFE(unsigned nodeIndex,
                                                           double damping,
                                                           c_vector<double, SPACE_DIM>& oldLocation,
                                                           double dt);

    /**
     * Helper method to return the expected step location for RK4NumericalMethod.
     *
     * @return the expected location after one step
     *
     * @param nodeIndex the index of the node
     * @param damping the damping constant
     * @param oldLocation the old location of the node
     * @param dt the step size
     */
    c_vector<double, SPACE_DIM> GetExpectedOneStepLocationRK4(unsigned nodeIndex,
                                                           double damping,
                                                           c_vector<double, SPACE_DIM>& oldLocation,
                                                           double dt);

    /**
     * Helper method to return the expected step location for AdamsMoultonNumericalMethod.
     *
     * @return the expected location after one step
     *
     * @param nodeIndex the index of the node
     * @param damping the damping constant
     * @param oldLocation the old location of the node
     * @param dt the step size
     */
    c_vector<double, SPACE_DIM> GetExpectedOneStepLocationAM2(unsigned nodeIndex,
                                                           double damping,
                                                           c_vector<double, SPACE_DIM>& oldLocation,
                                                           double dt);

    /**
     * Helper method to return the expected step location for BackwardEulerNumericalMethod.
     *
     * @return the expected location after one step
     *
     * @param nodeIndex the index of the node
     * @param damping the damping constant
     * @param oldLocation the old location of the node
     * @param dt the step size
     */
    c_vector<double, SPACE_DIM> GetExpectedOneStepLocationBE(unsigned nodeIndex,
                                                           double damping,
                                                           c_vector<double, SPACE_DIM>& oldLocation,
                                                           double dt);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PopulationTestingForce)

#endif /*POPULATIONTESTINGFORCE_HPP_*/
