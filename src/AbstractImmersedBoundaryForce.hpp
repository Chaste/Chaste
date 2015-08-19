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

#ifndef ABSTRACTIMMERSEDBOUNDARYFORCE_HPP_
#define ABSTRACTIMMERSEDBOUNDARYFORCE_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "ImmersedBoundaryCellPopulation.hpp"

/**
 * An abstract force class, for use in Immersed Boundary simulations.
 */
template<unsigned DIM>
class AbstractImmersedBoundaryForce : public Identifiable
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
    }

public:

    /**
     * Default constructor.
     */
    AbstractImmersedBoundaryForce();

    /**
     * Destructor.
     */
    virtual ~AbstractImmersedBoundaryForce();

    /**
     * Calculates the force on each node.
     *
     * This method must be overridden in concrete classes.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)=0;

    /**
     * Outputs force used in the simulation to file and then calls OutputForceParameters to output all relevant parameters.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceInfo(out_stream& rParamsFile);

    /**
     * Outputs force parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile)=0;
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractImmersedBoundaryForce)

#endif /*ABSTRACTIMMERSEDBOUNDARYFORCE_HPP_*/
