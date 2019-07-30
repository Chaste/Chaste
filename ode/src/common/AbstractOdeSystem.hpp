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

#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_

#include <vector>
#include <string>
#include <algorithm>


#include "ChasteSerialization.hpp"
#include "ChasteSerializationVersion.hpp"
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#include "ClassIsAbstract.hpp"

#include "AbstractParameterisedSystem.hpp"
#include "Exception.hpp"

/**
 * Abstract OdeSystem class.
 *
 * Sets up variables and functions for a general ODE system.
 *
 * ODE systems are specified primarily by the EvaluateYDerivatives() method,
 * which calculates the right-hand side of the system.
 *
 * Instances can store their state internally in the mStateVariables vector
 * in our base class AbstractParameterisedSystem (see also
 * GetNumberOfStateVariables(), SetStateVariables() and rGetStateVariables()),
 * although this is not essential - the vector may be empty, in which case
 * AbstractIvpOdeSolver::SolveAndUpdateStateVariable may not be used to
 * solve the system.
 *
 * ODE systems may also have a vector of parameters, which can be accessed
 * through the GetParameter() and SetParameter() methods of our base class.
 *
 * Information about what the parameters and state variables represent is
 * provided by a subclass of AbstractOdeSystemInformation.  Various wrapper
 * methods (e.g. rGetStateVariableNames()) are provided in our base class to
 * access this information.
 *
 * There are two more advanced facilities available for subclass authors.
 * An analytic form for the Jacobian matrix of the system may be provided,
 * in which case you must subclass AbstractOdeSystemWithAnalyticJacobian.
 * The GetUseAnalyticJacobian() method will test whether this is the case.
 *
 * Also, subclasses may define a condition at which ODE solvers should stop
 * prematurely.  For the Chaste solvers this is done by overriding
 * CalculateStoppingEvent(); if the more advanced CVODE solvers are being used
 * then implement CalculateRootFunction() instead to detect the stopping time
 * more accurately.
 */
class AbstractOdeSystem : public AbstractParameterisedSystem<std::vector<double> >
{
    friend class TestAbstractOdeSystem;

private:


    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        // Despite the fact that 3 of these variables actually live in our base class,
        // we still archive them here to maintain backwards compatibility.
        // Since the N_Vector version of mStateVariables and mParameters needs converting
        // to a standard vector before archiving, this doesn't hurt too much.
        archive & mNumberOfStateVariables;
        archive & mUseAnalyticJacobian;
        archive & mStateVariables;
        archive & mParameters;

        if (version > 0)
        {
            archive & rGetParameterNames();
        }

        // This is always set up by subclass constructors, and is essentially
        // 'static' data, so shouldn't go in the archive.
        //archive &mpSystemInfo;
    }
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & mNumberOfStateVariables;
        archive & mUseAnalyticJacobian;
        archive & mStateVariables;
        std::vector<double> parameters;
        archive & parameters;

        if (version > 0)
        {
            std::vector<std::string> param_names;
            archive & param_names;

            CheckParametersOnLoad(parameters,param_names);
        }
        else
        {
            mParameters = parameters;
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

protected:

    /** Whether to use an analytic Jacobian. */
    bool mUseAnalyticJacobian;

public:

    /**
     * Constructor.
     *
     * @param numberOfStateVariables  the number of state variables in the ODE system
     */
    AbstractOdeSystem(unsigned numberOfStateVariables);

    /**
     * Virtual destructor since we have virtual methods.
     */
    virtual ~AbstractOdeSystem();

    /**
     * Method to evaluate the derivatives of the system.
     *
     * @param time  the current time
     * @param rY  the current values of the state variables
     * @param rDY  storage for the derivatives of the system; will be filled in on return
     */
    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                                      std::vector<double>& rDY)=0;

    /**
     * CalculateStoppingEvent() - can be overloaded if the ODE is to be solved
     * only until a particular event (for example, only until the y value becomes
     * negative.
     *
     * After each timestep the solver will call this method on the ODE to see if
     * it should stop there.
     * @return true if the solver should stop now.  By default, false is returned here.
     *
     * @param time  the current time
     * @param rY  the current values of the state variables
     */
    virtual bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    /**
     * An alternative approach to stopping events; currently only useful with CVODE.
     * CVODE can search for roots (zeros) of this function while solving the ODE system,
     * and home in on them to find sign transitions to high precision.
     *
     * The default implementation here fakes a root function using CalculateStoppingEvent.
     *
     * @param time  the current time
     * @param rY  the current values of the state variables
     * @return value of the root function
     */
    virtual double CalculateRootFunction(double time, const std::vector<double>& rY);

    /**
     * Get whether an analytic Jacobian is used.
     *
     * @return #mUseAnalyticJacobian
     */
    bool GetUseAnalyticJacobian();

    /**
     * \todo move to AbstractParameterisedSystem? (1540)
     *
     * @return const reference to the state variables in the ODE system (used in archiving).
     */
    const std::vector<double>& rGetConstStateVariables() const;
};

CLASS_IS_ABSTRACT(AbstractOdeSystem)
BOOST_CLASS_VERSION(AbstractOdeSystem, 1u)

#endif //_ABSTRACTODESYSTEM_HPP_
