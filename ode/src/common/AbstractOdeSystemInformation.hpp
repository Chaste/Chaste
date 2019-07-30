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


#ifndef _ABSTRACTODESYSTEMINFORMATION_HPP_
#define _ABSTRACTODESYSTEMINFORMATION_HPP_

#include <vector>
#include <string>
#include <map>

/**
 * An abstract class which provides access to information about a particular
 * ODE system *class* (as opposed to an instance).
 *
 * The information available includes:
 *  - a name for the system
 *  - name & units of the free variable
 *  - names & units of state variables
 *  - suggested initial conditions
 *  - names & units of (settable) parameters
 *  - names & units of derived quantities
 *
 * This class requires a subclass defining the Initialise method in order to set
 * up the information.  Developers may do this by defining their own subclass, but
 * the most convenient method is likely to be to use the OdeSystemInformation
 * class, which is a templated singleton subclass of this AbstractOdeSystemInformation
 * class.  See its documentation for details of how to use it.
 */
class AbstractOdeSystemInformation
{
    friend class TestAbstractOdeSystem;

protected:
    /** Human-friendly name for the ODE system */
    std::string mSystemName;

    /** The name of the free variable. */
    std::string mFreeVariableName;

    /** The units of the free variable. */
    std::string mFreeVariableUnits;

    /** State variable names */
    std::vector<std::string> mVariableNames;

    /** State variable units */
    std::vector<std::string> mVariableUnits;

    /** Parameter names */
    std::vector<std::string> mParameterNames;

    /** Parameter units */
    std::vector<std::string> mParameterUnits;

    /** Derived quantity names */
    std::vector<std::string> mDerivedQuantityNames;

    /** Derived quantity units */
    std::vector<std::string> mDerivedQuantityUnits;

    /** Named attributes */
    std::map<std::string, double> mAttributes;

    /** Suggested initial conditions */
    std::vector<double> mInitialConditions;

    /** Whether a 'real' Initialise method has been called */
    bool mInitialised;

    /**
     * Initialise the ODE system information.
     *
     * This must be provided by subclasses.
     */
    virtual void Initialise()=0;

public:

    /**
     * Constructor.
     */
    AbstractOdeSystemInformation();

    /**
     * Virtual destructor since we have virtual methods.
     */
    virtual ~AbstractOdeSystemInformation();

    /**
     * @return the name of this system of ODEs.
     */
    std::string GetSystemName() const;

    /**
     * @return the name of the free variable in this system of ODEs.
     */
    std::string GetFreeVariableName() const;

    /**
     * @return the units of the free variable in this system of ODEs.
     */
    std::string GetFreeVariableUnits() const;

    /**
     * Set the default initial conditions for the ODE system. This method DOES NOT change the
     * state variables of the ODE object on which it is called.
     *
     * @param rInitialConditions  vector containing initial values for the state variables
     */
    void SetDefaultInitialConditions(const std::vector<double>& rInitialConditions);

    /**
     * Set a single component of the default initial conditions for the ODE system. This method
     * DOES NOT change the state variables of the ODE object on which it is called.
     *
     * @param index  the index of the state variable in the system
     * @param initialCondition  the initial value for the state variable
     */
    void SetDefaultInitialCondition(unsigned index, double initialCondition);

    /**
     * @return a copy of the suggested initial conditions.
     */
    std::vector<double> GetInitialConditions() const;

    /**
     * @return the state variable names vector.
     */
    const std::vector<std::string>& rGetStateVariableNames() const;

    /**
     * @return the state variable units vector.
     */
    const std::vector<std::string>& rGetStateVariableUnits() const;

    /**
     * This method is used to establish a state varible's position within
     * the vector of state variables of an ODE system. This number can
     * then be used with the method GetStateVariableUnits.
     *
     * @param rName  the name of a state variable
     * @return the state variable's position within the vector of state
     *         variables associated with the ODE system.
     */
    unsigned GetStateVariableIndex(const std::string& rName) const;

    /**
     * This method is used to establish whether a state variable is in
     * an ODE system. You can then safely call GetStateVariableIndex
     * without a try...catch statement.
     *
     * @param rName  the name of a state variable
     * @return whether the state variable is in this ODE system
     */
    bool HasStateVariable(const std::string& rName) const;

    /**
     * @return the units of a state variable given its index in the ODE system.
     *
     * @param index  a state variable's position within the vector of
     *               state variables associated with the ODE system.
     * @return the units of the state variable.
     */
    std::string GetStateVariableUnits(unsigned index) const;

    /**
     * @return the vector of parameter names.
     */
    const std::vector<std::string>& rGetParameterNames() const;

    /**
     * @return the vector of parameter units.
     */
    const std::vector<std::string>& rGetParameterUnits() const;

    /**
     * This method is used to establish a parameter's position within
     * the vector of parameters of an ODE system. This number can
     * then be used with the method GetParameterUnits.
     *
     * @param rName  the name of a parameter
     * @return the parameter's position within the vector of parameters
     *         associated with the ODE system.
     */
    unsigned GetParameterIndex(const std::string& rName) const;

    /**
     * This method is used to establish whether a parameter is in
     * an ODE system. You can then safely call GetParameterIndex
     * without a try...catch statement.
     *
     * @param rName  the name of a parameter
     * @return whether the parameter is in this ODE system
     */
    bool HasParameter(const std::string& rName) const;

    /**
     * Get the units of a parameter given its index in the ODE system.
     *
     * @param index  a state variable's position within the vector of
     *               state variables associated with the ODE system.
     * @return the units of the state variable.
     */
    std::string GetParameterUnits(unsigned index) const;

    /**
     * @return The number of parameters in this parameterised system.
     */
    unsigned GetNumberOfParameters() const;

    /**
     * @return the index of a variable, whether a state variable, parameter,
     * or derived quantity, with the given name.  The returned index is
     * suitable for use with GetAnyVariableUnits.
     *
     * Indices go through Variables, Parameters then derived quantities.
     *
     * @param rName  the name of a variable
     */
    unsigned GetAnyVariableIndex(const std::string& rName) const;

    /**
     * This method is used to establish whether a variable is in
     * an ODE system's state vars, parameters or derived quantitites.
     * You can then safely call GetAnyVariableIndex
     * without a try...catch statement.
     *
     * @param rName  the name of a variable
     * @return whether the variable is in this ODE system
     */
    bool HasAnyVariable(const std::string& rName) const;

    /**
     * Get the units of a variable, whether a state variable, parameter, or
     * derived quantity, given its index as returned by GetAnyVariableIndex.
     *
     * @param index  an index from GetAnyVariableIndex.
     * @return the units of the variable.
     */
    std::string GetAnyVariableUnits(unsigned index) const;

    /**
     * @return the vector of derived quantity names.
     */
    const std::vector<std::string>& rGetDerivedQuantityNames() const;

    /**
     * @return the vector of derived quantity units.
     */
    const std::vector<std::string>& rGetDerivedQuantityUnits() const;

    /**
     * @return the index of a derived quantity, given its name.
     *
     * @param rName  the name of a derived quantity.
     */
    unsigned GetDerivedQuantityIndex(const std::string& rName) const;

    /**
     * This method is used to establish whether a derived quantity is in
     * an ODE system. You can then safely call GetDerivedQuantityIndex
     * without a try...catch statement.
     *
     * @param rName  the name of a derived quantity
     * @return whether the derived quantity is in this ODE system
     */
    bool HasDerivedQuantity(const std::string& rName) const;

    /**
     * Get the units of a derived quantity.
     *
     * @param index  an index from GetDerivedQuantityIndex.
     * @return the units of the variable.
     */
    std::string GetDerivedQuantityUnits(unsigned index) const;

    /**
     * @return the number of derived quantities in this system
     */
    unsigned GetNumberOfDerivedQuantities() const;

    //
    // Attribute methods
    //

    /**
     * @return the number of named attributes that this system has.
     */
    unsigned GetNumberOfAttributes() const;

    /**
     * @return true if this system has a particular named attribute.
     * @param rName  the attribute name.
     */
    bool HasAttribute(const std::string& rName) const;

    /**
     * @return the value of a named attribute.
     * @param rName  the attribute name.
     */
    double GetAttribute(const std::string& rName) const;
};

#endif /*_ABSTRACTODESYSTEMINFORMATION_HPP_*/
