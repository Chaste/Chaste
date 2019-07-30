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

#ifndef COMBINEDODESYSTEM_HPP_
#define COMBINEDODESYSTEM_HPP_

#include <vector>
#include <map>

#include "AbstractOdeSystem.hpp"

/**
 * An ODE system formed by combining several subsystems.
 *
 * Instances of this class are formed by passing other ODE system instances
 * to our constructor, and then calling Configure to set up any coupling
 * between subsystems.  This allows state variables from one system to be
 * used as parameters in other systems.
 *
 * This allows us to set up coupled systems such as
 *   \f$dy/dt = f(y, t)\f$
 * from subsystems
 *   \f$dy_1/dt = f_1(y_1, t; y_2)\f$
 * and
 *   \f$dy_2/dt = f_2(y_2, t; y_1)\f$
 * where
 *   \f$y = (y_1, y_2)\f$.
 *
 * The vector of state variables for the combined system is formed as the
 * concatenation of the state variable vectors of the subsystem.  This
 * class also makes use of the CombinedOdeSystemInformation class to
 * provide initial conditions, etc.
 */
class CombinedOdeSystem : public AbstractOdeSystem
{
private:

    /** The subsystems forming this combined system. */
    std::vector<AbstractOdeSystem*> mOdeSystems;

    /** Working memory for the Y vectors of the subsystems. */
    std::vector<std::vector<double> > mWorkingStateVars;
    /** Working memory for the DY vectors of the subsystems. */
    std::vector<std::vector<double> > mWorkingDerivs;
    /**
     * Keeps track of where the state variable vector for each subsystem
     * is located within the combined state variable vector.
     */
    std::vector<unsigned> mOffsets;

    /**
     * A convenience structure for recording the information passed
     * to Configure calls.
     */
    struct VariableParameterMap
    {
        /**
         * A map from state variable index to parameter index
         * (within the vectors of state variables and parameters, respectively).
         */
        std::map<unsigned, unsigned> theMap;
        /**
         * The index of the system whose state variables #theMap refers to
         * within our vector of subsystems.
         */
        unsigned pVariableOdeSystemIndex;
        /** The ODE system whose parameters #theMap refers to. */
        AbstractOdeSystem* pParameterOdeSystem;
    };

    /** Stores the information passed to Configure calls. */
    std::vector<struct VariableParameterMap> mVariableParameterMaps;

public:
    /**
     * Create a combined ODE system from a vector of subsystems.
     *
     * @param odeSystems  the subsystems.
     */
    CombinedOdeSystem(std::vector<AbstractOdeSystem*> odeSystems);

    /**
     * Configure a mapping between the state variables of one subsystem and
     * the parameters of another.
     *
     * @param rVariableParameterMap  a map specifying which state variables (keys)
     *    are mapped to which parameters (values).
     * @param pVariableOdeSystem  the ODE subsystem providing state variable values.
     * @param pParameterOdeSystem  the ODE subsystem whose parameters should be set.
     */
    void Configure(const std::map<unsigned, unsigned>& rVariableParameterMap,
                   AbstractOdeSystem* pVariableOdeSystem,
                   AbstractOdeSystem* pParameterOdeSystem);

    /**
     * Evaluate the right-hand side of the combined system.
     *
     * This calls EvaluateYDerivatives for each subsystem with the appropriate
     * portion of rY and rDY, having set parameters from values in rY according
     * to the configured maps.
     *
     * @param time  the current time
     * @param rY  the current values of the state variables
     * @param rDY  storage for the derivatives of the system; will be filled in on return
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);
};


#endif /*COMBINEDODESYSTEM_HPP_*/
