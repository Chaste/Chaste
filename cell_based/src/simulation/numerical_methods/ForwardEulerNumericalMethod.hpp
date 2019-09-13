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

#ifndef FORWARDEULERNUMERICALMETHOD_HPP_
#define FORWARDEULERNUMERICALMETHOD_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractNumericalMethod.hpp"

/**
 * Implements forward Euler time stepping.
 *
 * Solves the equations of motion dr/dt = F
 * Using the scheme
 *
 * r^(t+1) = r^t + dt F^t.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class ForwardEulerNumericalMethod : public AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM> {

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Save or restore the simulation.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM> >(*this);
    }

public:

    /**
     * Constructor.
     */
    ForwardEulerNumericalMethod();

    /**
     * Destructor.
     */
    virtual ~ForwardEulerNumericalMethod();

    /**
     * Overridden UpdateAllNodePositions() method.
     *
     * @param dt Time step size
     */
    void UpdateAllNodePositions(double dt);

    /**
     * Overridden OutputNumericalMethodParameters() method.
     *
     * @param rParamsFile Reference to the parameter output filestream
     */
    virtual void OutputNumericalMethodParameters(out_stream& rParamsFile);
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ForwardEulerNumericalMethod)

#endif /*FORWARDEULERNUMERICALMETHOD_HPP_*/
