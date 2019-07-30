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

#ifndef _MOCKEULERIVPODESOLVER_HPP_
#define _MOCKEULERIVPODESOLVER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "EulerIvpOdeSolver.hpp"

/**
 * This 'mock' class is only used in testing. It is the same
 * as the EulerIvpOdeSolver class, but also keeps a count of
 * how many times it has been called. This is useful to check
 * ODE solving has been parallelised correctly.
 */
class MockEulerIvpOdeSolver : public EulerIvpOdeSolver
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the abstract IVP Solver, never used directly - boost uses this.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<EulerIvpOdeSolver>(*this);
        archive & mCallCount;
    }

    /** How many times the ODE solver has been called. */
    unsigned mCallCount;

protected:

    /**
     * Method that actually performs the solving on behalf of the public Solve methods.
     *
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param rCurrentYValues  the current (initial) state; results will also be returned
     *                         in here
     * @param rWorkingMemory  working memory; same size as rCurrentYValues
     * @param startTime  initial time
     * @param endTime  time to solve to
     * @param timeStep  dt
     */
    virtual void InternalSolve(AbstractOdeSystem* pAbstractOdeSystem,
                               std::vector<double>& rCurrentYValues,
                               std::vector<double>& rWorkingMemory,
                               double startTime,
                               double endTime,
                               double timeStep);

public:

    /**
     * Constructor.
     */
    MockEulerIvpOdeSolver();

    /**
     * Get the number of times the ODE solver has been called.
     *
     * @return #mCallCount.
     */
    unsigned GetCallCount();

    /**
     * Destructor.
     */
    virtual ~MockEulerIvpOdeSolver()
    {}
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MockEulerIvpOdeSolver)

#endif //_MOCKEULERIVPODESOLVER_HPP_
