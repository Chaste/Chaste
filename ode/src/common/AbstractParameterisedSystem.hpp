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

#ifndef ABSTRACTPARAMETERISEDSYSTEM_HPP_
#define ABSTRACTPARAMETERISEDSYSTEM_HPP_

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
     * Get the object which provides information about this ODE system.
     */
    boost::shared_ptr<const AbstractOdeSystemInformation> GetSystemInformation() const;

    /**
     * Get the name of this system.
     */
    std::string GetSystemName() const;

    //
    // Attribute methods
    //

    /**
     * Return the number of named attributes that this system has.
     */
    unsigned GetNumberOfAttributes() const;

    /**
     * Test whether this system has a particular named attribute.
     * @param rName  the attribute name.
     */
    bool HasAttribute(const std::string& rName) const;

    /**
     * Get the value of a named attribute.
     * @param rName  the attribute name.
     */
    double GetAttribute(const std::string& rName) const;

    //
    // State variable methods
    //

    /**
     * Get the number of state variables in the ODE system.
     *
     * @return mNumberOfStateVariables
     */
    unsigned GetNumberOfStateVariables() const;

    /**
     * Get the names of the state variables in the ODE system.
     */
    const std::vector<std::string>& rGetStateVariableNames() const;

    /**
     * Get the units of the state variables in the ODE system.
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
     * Get the units of a state variable given its index in the ODE system.
     *
     * @param index  a state variable's position within the vector of
     *               state variables associated with the ODE system.
     * @return the units of the state variable.
     */
    std::string GetStateVariableUnits(unsigned index) const;

    //
    // Parameter methods
    //

    /**
     * Get the number of parameters.
     */
    unsigned GetNumberOfParameters() const;

    /**
     * Get the names of the parameters in the ODE system.
     */
    const std::vector<std::string>& rGetParameterNames() const;

    /**
     * Get the units of the parameters in the ODE system.
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
     * Get the units of a parameter given its index in the ODE system.
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
     * Get the number of derived quantities.
     */
    unsigned GetNumberOfDerivedQuantities() const;

    /**
     * Get the vector of derived quantity names.
     */
    const std::vector<std::string>& rGetDerivedQuantityNames() const;

    /**
     * Get the vector of derived quantity units.
     */
    const std::vector<std::string>& rGetDerivedQuantityUnits() const;

    /**
     * Get the index of a derived quantity, given its name.
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

    //
    // "Any variable" methods
    //

    /**
     * Get the index of a variable, whether a state variable, parameter,
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
     * Get the units of a variable, whether a state variable, parameter, or
     * derived quantity, given its index as returned by GetAnyVariableIndex.
     *
     * @param index  an index from GetAnyVariableIndex.
     * @return the units of the variable.
     */
    std::string GetAnyVariableUnits(unsigned index) const;

    /**
     * Get the units of a variable, whether a state variable, parameter, or
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
     * @param message  a string to prefix (e.g. an error or the name of the vector)
     * @param Y  a vector
     * @return  a string containing the contents of the vector.
     */
    std::string GetStateMessage(const std::string& message, VECTOR Y);

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
     */
    std::string DumpState(const std::string& rMessage);

    /**
     * Used to include extra debugging information in exception messages.
     * For example,
     *      EXCEPTION(DumpState("Gating variable out of range", state_variables));
     *
     * @param rMessage  the exception message
     * @param Y  the values of the state variables
     */
    std::string DumpState(const std::string& rMessage,
                          VECTOR Y);

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
     * Get the values of the state variables in the ODE system.
     */
    VECTOR& rGetStateVariables();

    /**
     * Get a copy of the state variable vector.
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
     * Get the value of a given state variable.
     *
     * @param index the index of the state variable
     */
    double GetStateVariable(unsigned index) const;

    /**
     * Get the value of a given state variable.
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
     * Get the default initial conditions for this system.
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
     * Get the value of a given parameter.
     *
     * @param index the index of the parameter
     */
    double GetParameter(unsigned index) const;

    /**
     * Get the value of a given parameter.
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
     * Get the value of a variable, whether a state variable, parameter,
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
     * Get the value of a variable, whether a state variable, parameter,
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
     */
    virtual VECTOR ComputeDerivedQuantities(double time,
                                            const VECTOR& rState);

    /**
     * Compute the derived quantities based on the current system state.
     *
     * @param time  the time at which to compute the derived quantities
     */
    VECTOR ComputeDerivedQuantitiesFromCurrentState(double time);
};



#endif /*ABSTRACTPARAMETERISEDSYSTEM_HPP_*/
