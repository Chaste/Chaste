/*

Copyright (C) University of Oxford, 2005-2011

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
     * @return mCallCount.
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
