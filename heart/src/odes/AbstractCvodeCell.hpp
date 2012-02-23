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

};


#endif // _ABSTRACTCVODECELL_HPP_
#endif // CHASTE_CVODE
