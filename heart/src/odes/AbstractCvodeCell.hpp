/*

Copyright (C) University of Oxford, 2005-2012

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
