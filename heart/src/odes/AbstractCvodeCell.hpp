/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifdef CHASTE_CVODE
#ifndef _ABSTRACTCVODECELL_HPP_
#define _ABSTRACTCVODECELL_HPP_

#include <boost/shared_ptr.hpp>

// Chaste headers
#include "AbstractOdeSystemInformation.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractCvodeSystem.hpp"
#include "AbstractCardiacCellInterface.hpp"
#include "VectorHelperFunctions.hpp"


/**
 * A cardiac cell that is designed to be simulated using CVODE.
 * It uses CVODE's vector type natively.
 *
 * Functionality is similar to that provided by AbstractCardiacCell and AbstractOdeSystem,
 * but not identical.  It also includes a direct interface to the CVODE solver, via the
 * Solve methods, since the CvodeAdaptor class doesn't work for us.
 *
 * Assumes that it will be solving stiff systems, so uses BDF/Newton.
 *
 * Note that a call to Solve will initialise the CVODE solver, and free its
 * working memory when done.  There is thus a non-trivial overhead involved.
 *
 * \todo #890 Add an option to just initialise once, and assume subsequent Solve
 *   calls are continuing from where we left off.
 */
class AbstractCvodeCell : public AbstractCardiacCellInterface, public AbstractCvodeSystem
{
private:
    /** The maximum timestep to use. */
    double mMaxDt;

public:
    /**
     * Create a new cardiac cell.
     *
     * @note subclasses @b must call Init() in their constructor after setting #mpSystemInfo.
     *
     * @param pSolver  not used for these cells (they're always solved with CVODE); may be empty
     * @param numberOfStateVariables  the size of the ODE system modelling this cell
     * @param voltageIndex  the index of the transmembrane potential within the vector of state variables
     * @param pIntracellularStimulus  the intracellular stimulus current
     */
    AbstractCvodeCell(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                      unsigned numberOfStateVariables,
                      unsigned voltageIndex,
                      boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Free the state variables, if they have been set.
     */
    virtual ~AbstractCvodeCell();

    /**
     * Get the current value of the transmembrane potential, as given
     * in our state variable vector.
     */
    double GetVoltage();

    /**
     * Set the transmembrane potential
     * @param voltage  new value
     */
    void SetVoltage(double voltage);

    /**
     * Set the maximum timestep to use for simulating this cell.
     *
     * @param maxDt  the maximum timestep
     */
    void SetTimestep(double maxDt);

    /**
     * Simulate this cell's behaviour between the time interval [tStart, tEnd],
     * updating the internal state variable values.
     *
     * The maximum time step to use is given by #mMaxDt, which defaults to
     * HeartConfig::Instance()->GetPrintingTimeStep() if unset.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    virtual void SolveAndUpdateState(double tStart, double tEnd);

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * and return state variable values.
     *
     * The maximum time step to use will be taken as #mMaxDt.  If this is unset
     * it is the same as the sampling interval, which defaults to
     * HeartConfig::Instance()->GetPrintingTimeStep().
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     * @param tSamp  sampling interval for returned results (defaults to dt)
     */
    OdeSolution Compute(double tStart, double tEnd, double tSamp=0.0);

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * but does not update the voltage.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    void ComputeExceptVoltage(double tStart, double tEnd);

    /**
     * Set whether to clamp the voltage by setting its derivative to zero.
     * We need to ensure CVODE is re-initialised if this setting changes.
     * @param clamp  whether to clamp
     */
    void SetVoltageDerivativeToZero(bool clamp=true);

    /**
     * This just returns the number of state variables in the cell model.
     *
     * It is here because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface
     *
     * @return the number of state variables
     */
    unsigned GetNumberOfStateVariables() const;

    /**
     * This just returns the number of parameters in the cell model.
     *
     * It is here because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface
     *
     * @return the number of parameters
     */
    unsigned GetNumberOfParameters() const;

    /**
     * This just returns the state variables in the cell model.
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @return the state variables
     */
    std::vector<double> GetStdVecStateVariables();


    /**
     * This just sets the state variables in the cell model.
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @param rVariables  the state variables (to take a copy of).
     */
    void SetStateVariables(const std::vector<double>& rVariables);

    /**
     * This is also needed now just to show there is an alternative to the above method!
     * @param rVariables  the state variables (to take a copy of).
     */
    void SetStateVariables(const N_Vector& rVariables);

    /**
     * This just calls the method AbstractCvodeSystem::GetAnyVariable
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @param rName variable name
     * @param time  the time at which to evaluate variable (only needed for derived quantities).
     * @return value of the variable at that time
     */
    double GetAnyVariable(const std::string& rName, double time=0.0);

    /**
     * This just calls AbstractCvodeSystem::GetParameter
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @param rParameterName  the name of a parameter to get the value of,
     * @return  the parameter's value.
     */
    double GetParameter(const std::string& rParameterName);

    /**
     * This is just here to show there is an alternative to the above!
     *
     * It just calls the base class method.
     *
     * @param parameterIndex  the index of the parameter vector entry to return
     * @return the parameter value
     */
    double GetParameter(unsigned parameterIndex);

    /**
     * This just calls AbstractCvodeSystem::SetParameter
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @param rParameterName  the parameter name to set the value of,
     * @param value  value to set it to.
     */
    void SetParameter(const std::string& rParameterName, double value);

    /**
     * This is just here to show there is an alternative to the above!
     *
     * It just calls the base class method.
     *
     * @param parameterIndex  the index of the parameter vector to alter
     * @param value  the value the parameter should take
     */
    void SetParameter(unsigned parameterIndex, double value);

};


#endif // _ABSTRACTCVODECELL_HPP_
#endif // CHASTE_CVODE
