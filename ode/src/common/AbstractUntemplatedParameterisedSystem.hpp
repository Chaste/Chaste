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

#ifndef ABSTRACTUNTEMPLATEDPARAMETERISEDSYSTEM_HPP_
#define ABSTRACTUNTEMPLATEDPARAMETERISEDSYSTEM_HPP_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "AbstractOdeSystemInformation.hpp"

/**
 * This class is an untemplated base class for AbstractParameterisedSystem, containing
 * those methods which don't require knowledge of the vector type, in order to make it
 * easier to move between templated and generic parts of the codebase.  In particular
 * it holds the AbstractOdeSystemInformation pointer, and methods to access this object
 * to provide information about the ODE system, such as state variable/parameter names
 * and units.
 */
class AbstractUntemplatedParameterisedSystem
{
public:
    /**
     * Constructor.
     *
     * @param numberOfStateVariables  the number of state variables in the ODE system
     */
    AbstractUntemplatedParameterisedSystem(unsigned numberOfStateVariables);

    /** Make this class polymorphic. */
    virtual ~AbstractUntemplatedParameterisedSystem();

    /**
     * @return the object which provides information about this ODE system.
     */
    boost::shared_ptr<const AbstractOdeSystemInformation> GetSystemInformation() const;

    /**
     * @return the name of this system.
     */
    std::string GetSystemName() const;

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

    //
    // State variable methods
    //

    /**
     * @return the number of state variables in the ODE system.
     *
     * @return #mNumberOfStateVariables
     */
    unsigned GetNumberOfStateVariables() const;

    /**
     * @return the names of the state variables in the ODE system.
     */
    const std::vector<std::string>& rGetStateVariableNames() const;

    /**
     * @return the units of the state variables in the ODE system.
     */
    const std::vector<std::string>& rGetStateVariableUnits() const;

    /**
     * This method is used to establish a state variable's position within
     * the vector of state variables of an ODE system.  This number can
     * then be used with the methods GetStateVariable and GetStateVariableUnits.
     *
     * @param rName  the name of a state variable.
     *
     * @return the state variable's position within the vector of state variables
     *         associated with the ODE system.
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
     * Reset the system's state variables to the default initial conditions.
     */
    virtual void ResetToInitialConditions()=0;

    //
    // Parameter methods
    //

    /**
     * @return the number of parameters.
     */
    unsigned GetNumberOfParameters() const;

    /**
     * @return the names of the parameters in the ODE system.
     */
    const std::vector<std::string>& rGetParameterNames() const;

    /**
     * @return the units of the parameters in the ODE system.
     */
    const std::vector<std::string>& rGetParameterUnits() const;

    /**
     * This method is used to establish a parameter's position within
     * the vector of parameters of an ODE system. This number can
     * then be used with the methods GetParameterUnits and GetParameter.
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
     * @return the units of a parameter given its index in the ODE system.
     *
     * @param index  a state variable's position within the vector of
     *               state variables associated with the ODE system.
     * @return the units of the state variable.
     */
    std::string GetParameterUnits(unsigned index) const;

    //
    // Derived quantity methods
    //

    /**
     * @return the number of derived quantities.
     */
    unsigned GetNumberOfDerivedQuantities() const;

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
     * @return the units of a derived quantity.
     *
     * @param index  an index from GetDerivedQuantityIndex.
     * @return the units of the variable.
     */
    std::string GetDerivedQuantityUnits(unsigned index) const;

    //
    // "Any variable" methods
    //

    /**
     * @return the index of a variable, whether a state variable, parameter,
     * or derived quantity, with the given name.
     * The returned index is suitable for use with GetAnyVariableUnits,
     * GetAnyVariable, etc.
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
     * @return the units of a variable, whether a state variable, parameter, or
     * derived quantity, given its index as returned by GetAnyVariableIndex.
     *
     * @param index  an index from GetAnyVariableIndex.
     * @return the units of the variable.
     */
    std::string GetAnyVariableUnits(unsigned index) const;

    /**
     * @return the units of a variable, whether a state variable, parameter, or
     * derived quantity, given its index as returned by GetAnyVariableIndex.
     *
     * @param rName  the name of any variable in the model.
     * @return the units of the variable.
     */
    std::string GetAnyVariableUnits(const std::string& rName) const;


protected:
    /** The number of state variables in the system. */
    unsigned mNumberOfStateVariables;

    /**
     * Information about the concrete ODE system class.
     *
     * Subclasses @b need to set this in their constructor to point to an instance
     * of a suitable class.  See for example the OdeSystemInformation class.
     */
    boost::shared_ptr<AbstractOdeSystemInformation> mpSystemInfo;
};

#endif /*ABSTRACTUNTEMPLATEDPARAMETERISEDSYSTEM_HPP_*/
