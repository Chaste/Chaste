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

#ifndef ABSTRACTPARAMETERISEDSYSTEM_HPP_
#define ABSTRACTPARAMETERISEDSYSTEM_HPP_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "AbstractUntemplatedParameterisedSystem.hpp"


/**
 * This class contains the state variable and parameter vectors for an ODE system,
 * along with methods to access these.
 *
 * Its main purpose is to be a common base class for both AbstractOdeSystem and
 * AbstractCvodeSystem, which require similar functionality but use different vector
 * types.
 */
template<typename VECTOR>
class AbstractParameterisedSystem : public AbstractUntemplatedParameterisedSystem
{
friend class TestAbstractCvodeSystem;

private:
    /**
     * Helper method to construct a string containing a dump of the vector
     *
     * @param rMessage  a string to prefix (e.g. an error or the name of the vector)
     * @param Y  a vector
     * @return  a string containing the contents of the vector.
     */
    std::string GetStateMessage(const std::string& rMessage, VECTOR Y);

protected:
    /** Vector containing the current values of the state variables. */
    VECTOR mStateVariables;

    /** Vector containing parameter values. */
    VECTOR mParameters;

    /**
     * Used to include extra debugging information in exception messages.
     * For example,
     *      EXCEPTION(DumpState("Gating variable out of range"));
     *
     * @param rMessage  the exception message
     * @return an augmented message which includes the values of the internal state variables
     */
    std::string DumpState(const std::string& rMessage);

    /**
     * Used to include extra debugging information in exception messages.
     * For example,
     *      EXCEPTION(DumpState("Gating variable out of range", state_variables));
     *
     * @param rMessage  the exception message
     * @param Y  the values of the state variables
     * @return an augmented message which includes the values of the state variables from Y
     */
    std::string DumpState(const std::string& rMessage,
                          VECTOR Y);

    /**
     * Used to include extra debugging information in exception messages.
     * For example,
     *      EXCEPTION(DumpState("Gating variable out of range", state_variables, time));
     *
     * @param rMessage  the exception message
     * @param Y  the values of the state variables
     * @param time  the independent variable (usually time).
     *
     * @return an augmented message which includes the values of the state variables from Y and the time.
     */
    std::string DumpState(const std::string& rMessage,
                          VECTOR Y,
                          double time);

    /**
     * This method is called by subclasses on completion of the load method of serialization.
     *
     * It checks that the parameters that were loaded match those that should be in the class, and fills in
     * with default values any that are missing. Hence this method updates mParameters.
     *
     * @param rParameters  the parameters that were loaded.
     * @param rParameterNames  the parameter names that were loaded.
     */
    void CheckParametersOnLoad(const std::vector<double>& rParameters,
                               const std::vector<std::string>& rParameterNames);

public:
    /**
     * Constructor.
     *
     * @param numberOfStateVariables  the number of state variables in the ODE system
     */
    AbstractParameterisedSystem(unsigned numberOfStateVariables);

    //
    // State variable methods
    //

    /**
     * @return the values of the state variables in the ODE system.
     */
    VECTOR& rGetStateVariables();

    /**
     * @return a copy of the state variable vector.
     * Caller takes responsibility for deleting the returned vector (if required for the VECTOR type).
     */
    VECTOR GetStateVariables();

    /**
     * Set the state variables equal to the values in the given vector, copying it.
     * Caller thus maintains responsibility for deleting the input vector (if required for the VECTOR type).
     *
     * @param rStateVariables  new values for the state variables
     */
    void SetStateVariables(const VECTOR& rStateVariables);

    /**
     * @return the value of a given state variable.
     *
     * @param index the index of the state variable
     */
    double GetStateVariable(unsigned index) const;

    /**
     * @return the value of a given state variable.
     *
     * @param rName the name of the state variable
     */
    double GetStateVariable(const std::string& rName) const;

    /**
     * Set the value of a single state variable in the ODE system.
     *
     * @param index index of the state variable to be set
     * @param newValue new value of the state variable
     */
    void SetStateVariable(unsigned index, double newValue);

    /**
     * Set the value of a single state variable in the ODE system.
     *
     * @param rName name of the state variable to be set
     * @param newValue new value of the state variable
     */
    void SetStateVariable(const std::string& rName, double newValue);

    /**
     * Empty method which can be over-ridden and used in solvers to
     * go through the current state vector and do range checking on the values
     * (e.g. check that concentrations are positive and probabilities are between
     * zero and one).
     *
     * This method is overridden with a currently commented out
     * method in AbstractCardiacCell which would be called by the
     * ComputeExceptVoltage method (in heart).
     *
     * This method is called by the AbstractCvodeSystem::Solve() method (in ode).
     */
    virtual void VerifyStateVariables()
    {}

    //
    // Initial condition methods
    //

    /**
     * Set the default initial conditions for the system.
     *
     * @note The default initial conditions are shared among all instances of the particular concrete
     *     system class.
     * @note This method DOES NOT change the state variables of the object on which it is called.
     *
     * @param rInitialConditions  vector containing initial values for the state variables
     */
    void SetDefaultInitialConditions(const VECTOR& rInitialConditions);

    /**
     * Set a single component of the default initial conditions for the system.
     *
     * @note The default initial conditions are shared among all instances of the particular concrete
     *     system class.
     * @note This method DOES NOT change the state variables of the object on which it is called.
     *
     * @param index  the index of the state variable in the system
     * @param initialCondition  the initial value for the state variable
     */
    void SetDefaultInitialCondition(unsigned index, double initialCondition);

    /**
     * @return the default initial conditions for this system.
     *
     * @note Returns a fresh vector, which the caller should delete if appropriate for the type.
     */
    VECTOR GetInitialConditions() const;

    /**
     * Reset the system's state variables to the default initial conditions.
     */
    void ResetToInitialConditions();

    //
    // Parameter methods
    //

    /**
     * @return the value of a given parameter.
     *
     * @param index the index of the parameter
     */
    double GetParameter(unsigned index) const;

    /**
     * @return the value of a given parameter.
     *
     * @param rName the name of the parameter
     */
    double GetParameter(const std::string& rName) const;

    /**
     * Set the value of a given parameter.
     *
     * @param rName the name of the parameter
     * @param value the value
     */
    void SetParameter(const std::string& rName, double value);

    /**
     * Set the value of a given parameter.
     *
     * @param index the index of the parameter
     * @param value the value
     */
    void SetParameter(unsigned index, double value);

    //
    // "Any variable" methods
    //

    /**
     * @return the value of a variable, whether a state variable, parameter,
     * or derived quantity.
     *
     * Note that if the variable is a derived quantity, this method will compute
     * all derived quantities, so may not be very efficient.  To avoid this, pass
     * a pre-computed vector of derived quantities as the optional third argument.
     *
     * @param index the index of the variable, as given by GetAnyVariableIndex.
     * @param time  the current simulation time, possibly needed if the variable
     *     is a derived quantity.
     * @param pDerivedQuantities  optional vector of pre-computed derived quantity values.
     */
    double GetAnyVariable(unsigned index, double time=0.0,
                          VECTOR* pDerivedQuantities=NULL);

    /**
     * @return the value of a variable, whether a state variable, parameter,
     * or derived quantity.
     *
     * Note that if the variable is a derived quantity, this method will compute
     * all derived quantities, so may not be very efficient.  To avoid this, pass
     * a pre-computed vector of derived quantities as the optional third argument.
     *
     * @param rName the name of the variable, (this method is the same as doing GetAnyVariableIndex(rName) and then calling the method above).
     * @param time  the current simulation time, possibly needed if the variable
     *     is a derived quantity.
     * @param pDerivedQuantities  optional vector of pre-computed derived quantity values.
     */
    double GetAnyVariable(const std::string& rName, double time=0.0,
                          VECTOR* pDerivedQuantities=NULL);

    /**
     * Set the value of a variable, whether a state variable or parameter.
     * Attempting to set the value of a derived quantity will raise an exception.
     *
     * @param index  the index of the variable, as given by GetAnyVariableIndex.
     * @param value  the value to give the variable.
     */
    void SetAnyVariable(unsigned index, double value);

    /**
     * Set the value of a variable, whether a state variable or parameter.
     * Attempting to set the value of a derived quantity will raise an exception.
     *
     * @param rName  the name of the variable.
     * @param value  the value to give the variable.
     */
    void SetAnyVariable(const std::string& rName, double value);

    //
    // Derived quantity methods
    //

    /**
     * Compute the derived quantities from the given system state.
     * Uses the current values for the parameters.
     *
     * @param time  the time at which to compute the derived quantities
     * @param rState  values for the state variables
     * @return a VECTOR of derived quantities
     */
    virtual VECTOR ComputeDerivedQuantities(double time,
                                            const VECTOR& rState);

    /**
     * Compute the derived quantities based on the current system state.
     *
     * @param time  the time at which to compute the derived quantities
     * @return a VECTOR of derived quantities
     */
    VECTOR ComputeDerivedQuantitiesFromCurrentState(double time);
};

#endif /*ABSTRACTPARAMETERISEDSYSTEM_HPP_*/
