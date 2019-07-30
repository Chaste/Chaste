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


#ifndef ABSTRACTBACKWARDEULERCARDIACCELL_HPP_
#define ABSTRACTBACKWARDEULERCARDIACCELL_HPP_

#include <cassert>
#include <cmath>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "ClassIsAbstract.hpp"

#include "AbstractCardiacCell.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "TimeStepper.hpp"

/**
 * This is the base class for cardiac cells solved using a (decoupled) backward
 * Euler approach (see http://dx.doi.org/10.1109/TBME.2006.879425).
 *
 * The basic approach to solving such models is:
 *  \li Update the transmembrane potential, either from solving an external PDE,
 *      or using a forward Euler step.
 *  \li Update any gating variables (or similar) using a backward euler step.
 *      Suitable ODEs can be written in the form  \f$du/dt = g(V) + h(V)*u\f$.  The update
 *      expression is then \f$u_n = ( u_{n-1} + g(V_n)*dt ) / ( 1 - h(V_n)*dt )\f$.
 *  \li Update the remaining state variables using Newton's method to solve the
 *      nonlinear system \f$U_n - U_{n-1} = dt*F(U_n, V_n)\f$.
 *      The template parameter to the class specifies the size of this nonlinear system.
 */
template<unsigned SIZE>
class AbstractBackwardEulerCardiacCell : public AbstractCardiacCell
{
    private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractCardiacCell>(*this);
    }
public:

    /**
     * Standard constructor for a cell.
     *
     * @param numberOfStateVariables  the size of the ODE system
     * @param voltageIndex  the index of the variable representing the transmembrane
     *     potential within the state variable vector
     * @param pIntracellularStimulus  the intracellular stimulus function
     *
     * Some notes for future reference:
     *  \li We may want to remove the timestep from this class, and instead pass it to
     *      the Compute* methods, especially if variable timestepping is to be used.
     *  \li It's a pity that inheriting from AbstractCardiacCell forces us to store a
     *      null pointer (for the unused ODE solver) in every instance.  We may want
     *      to revisit this design decision at a later date.
     */
    AbstractBackwardEulerCardiacCell(
        unsigned numberOfStateVariables,
        unsigned voltageIndex,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /** Virtual destructor */
    virtual ~AbstractBackwardEulerCardiacCell();

    /**
     * Compute the residual of the nonlinear system portion of the cell model.
     *
     * @param time  the current time
     * @param rCurrentGuess  the current guess for \f$U_n\f$
     * @param rResidual  to be filled in with the residual vector
     */
    virtual void ComputeResidual(double time, const double rCurrentGuess[SIZE], double rResidual[SIZE])=0;

    /**
     * Compute the Jacobian matrix for the nonlinear system portion of the cell model.
     *
     * @param time  the current time
     * @param rCurrentGuess  the current guess for \f$U_n\f$
     * @param rJacobian  to be filled in with the Jacobian matrix
     */
    virtual void ComputeJacobian(double time, const double rCurrentGuess[SIZE], double rJacobian[SIZE][SIZE])=0;

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt.  Uses a forward Euler step to update the transmembrane
     * potential at each timestep.
     *
     * The length of the time interval must be a multiple of the timestep.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     * @param tSamp  sampling interval for returned results (defaults to #mDt)
     * @return  the values of each state variable, at intervals of tSamp.
     */
    OdeSolution Compute(double tStart, double tEnd, double tSamp=0.0);

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt.  The transmembrane potential is kept fixed throughout.
     *
     * The length of the time interval must be a multiple of the timestep.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    void ComputeExceptVoltage(double tStart, double tEnd);

    /**
     * Simulate this cell's behaviour between the time interval [tStart, tEnd],
     * with timestemp #mDt, updating the internal state variable values.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    void SolveAndUpdateState(double tStart, double tEnd);

private:
// LCOV_EXCL_START
    /**
     * This function should never be called - the cell class incorporates its own solver.
     *
     * @param time
     * @param rY
     * @param rDY
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
    {
        NEVER_REACHED;
    }
// LCOV_EXCL_STOP

protected:
    /**
     * Compute the values of all state variables, except the voltage, using backward Euler,
     * for one timestep from tStart.
     *
     * \note This method must be provided by subclasses.
     *
     * @param tStart  start of this timestep
     */
    virtual void ComputeOneStepExceptVoltage(double tStart)=0;

    /**
     * Perform a forward Euler step to update the transmembrane potential.
     *
     * \note This method must be provided by subclasses.
     *
     * @param time  start of this timestep
     */
    virtual void UpdateTransmembranePotential(double time)=0;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*
 * NOTE: Explicit instantiation is not used for this class, because the SIZE
 * template parameter could take arbitrary values.
 */


template <unsigned SIZE>
AbstractBackwardEulerCardiacCell<SIZE>::AbstractBackwardEulerCardiacCell(
    unsigned numberOfStateVariables,
    unsigned voltageIndex,
    boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCardiacCell(boost::shared_ptr<AbstractIvpOdeSolver>(),
                              numberOfStateVariables,
                              voltageIndex,
                              pIntracellularStimulus)
{}

template <unsigned SIZE>
AbstractBackwardEulerCardiacCell<SIZE>::~AbstractBackwardEulerCardiacCell()
{}

template <unsigned SIZE>
OdeSolution AbstractBackwardEulerCardiacCell<SIZE>::Compute(double tStart, double tEnd, double tSamp)
{
    // In this method, we iterate over timesteps, doing the following for each:
    //   - update V using a forward Euler step
    //   - call ComputeExceptVoltage(t) to update the remaining state variables
    //     using backward Euler

    // Check length of time interval
    if (tSamp < mDt)
    {
        tSamp = mDt;
    }
    double _n_steps = (tEnd - tStart) / tSamp;
    const unsigned n_steps = (unsigned) floor(_n_steps+0.5);
    assert(fabs(tStart+n_steps*tSamp - tEnd) < 1e-12);
    const unsigned n_small_steps = (unsigned) floor(tSamp/mDt+0.5);
    assert(fabs(mDt*n_small_steps - tSamp) < 1e-12);

    // Initialise solution store
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(n_steps);
    solutions.rGetSolutions().push_back(rGetStateVariables());
    solutions.rGetTimes().push_back(tStart);
    solutions.SetOdeSystemInformation(this->mpSystemInfo);

    // Loop over time
    double curr_time = tStart;
    for (unsigned i=0; i<n_steps; i++)
    {
        for (unsigned j=0; j<n_small_steps; j++)
        {
            curr_time = tStart + i*tSamp + j*mDt;

            // Compute next value of V
            UpdateTransmembranePotential(curr_time);

            // Compute other state variables
            ComputeOneStepExceptVoltage(curr_time);

            // check gating variables are still in range
            VerifyStateVariables();
        }

        // Update solutions
        solutions.rGetSolutions().push_back(rGetStateVariables());
        solutions.rGetTimes().push_back(curr_time+mDt);
    }

    return solutions;
}

template <unsigned SIZE>
void AbstractBackwardEulerCardiacCell<SIZE>::ComputeExceptVoltage(double tStart, double tEnd)
{
    // This method iterates over timesteps, calling ComputeExceptVoltage(t) at
    // each one, to update all state variables except for V, using backward Euler.

    // Check length of time interval
    unsigned n_steps = (unsigned)((tEnd - tStart) / mDt + 0.5);
    assert(fabs(tStart + n_steps*mDt - tEnd) < 1e-12);

    // Loop over time
    double curr_time;
    for (unsigned i=0; i<n_steps; i++)
    {
        curr_time = tStart + i*mDt;

        // Compute other state variables
        ComputeOneStepExceptVoltage(curr_time);

#ifndef NDEBUG
        // Check gating variables are still in range
        VerifyStateVariables();
#endif // NDEBUG
    }
}

template<unsigned SIZE>
void AbstractBackwardEulerCardiacCell<SIZE>::SolveAndUpdateState(double tStart, double tEnd)
{
    TimeStepper stepper(tStart, tEnd, mDt);

    while (!stepper.IsTimeAtEnd())
    {
        double time = stepper.GetTime();

        // Compute next value of V
        UpdateTransmembranePotential(time);

        // Compute other state variables
        ComputeOneStepExceptVoltage(time);

        // check gating variables are still in range
        VerifyStateVariables();

        stepper.AdvanceOneTimeStep();
    }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * Specialization for the case where there are no non-linear ODEs in the model.
 */
template<>
class AbstractBackwardEulerCardiacCell<0u> : public AbstractCardiacCell
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractCardiacCell>(*this);
    }

public:
    /**
     * Standard constructor for a cell.
     *
     * @param numberOfStateVariables  the size of the ODE system
     * @param voltageIndex  the index of the variable representing the transmembrane
     *     potential within the state variable vector
     * @param pIntracellularStimulus  the intracellular stimulus function
     */
    AbstractBackwardEulerCardiacCell(unsigned numberOfStateVariables,
                                     unsigned voltageIndex,
                                     boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCardiacCell(boost::shared_ptr<AbstractIvpOdeSolver>(),
                              numberOfStateVariables,
                              voltageIndex,
                              pIntracellularStimulus)
    {}

    /** Virtual destructor */
    virtual ~AbstractBackwardEulerCardiacCell()
    {}

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt.  Uses a forward Euler step to update the transmembrane
     * potential at each timestep.
     *
     * The length of the time interval must be a multiple of the timestep.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     * @param tSamp  sampling interval for returned results (defaults to #mDt)
     * @return  the values of each state variable, at intervals of tSamp.
     */
    OdeSolution Compute(double tStart, double tEnd, double tSamp=0.0)
    {
        // In this method, we iterate over timesteps, doing the following for each:
        //   - update V using a forward Euler step
        //   - call ComputeExceptVoltage(t) to update the remaining state variables
        //     using backward Euler

        // Check length of time interval
        if (tSamp < mDt)
        {
            tSamp = mDt;
        }
        double _n_steps = (tEnd - tStart) / tSamp;
        const unsigned n_steps = (unsigned) floor(_n_steps+0.5);
        assert(fabs(tStart+n_steps*tSamp - tEnd) < 1e-12);
        const unsigned n_small_steps = (unsigned) floor(tSamp/mDt+0.5);
        assert(fabs(mDt*n_small_steps - tSamp) < 1e-12);

        // Initialise solution store
        OdeSolution solutions;
        solutions.SetNumberOfTimeSteps(n_steps);
        solutions.rGetSolutions().push_back(rGetStateVariables());
        solutions.rGetTimes().push_back(tStart);
        solutions.SetOdeSystemInformation(this->mpSystemInfo);

        // Loop over time
        double curr_time = tStart;
        for (unsigned i=0; i<n_steps; i++)
        {
            for (unsigned j=0; j<n_small_steps; j++)
            {
                curr_time = tStart + i*tSamp + j*mDt;

                // Compute next value of V
                UpdateTransmembranePotential(curr_time);

                // Compute other state variables
                ComputeOneStepExceptVoltage(curr_time);

                // check gating variables are still in range
                VerifyStateVariables();
            }

            // Update solutions
            solutions.rGetSolutions().push_back(rGetStateVariables());
            solutions.rGetTimes().push_back(curr_time+mDt);
        }

        return solutions;
    }

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt.  The transmembrane potential is kept fixed throughout.
     *
     * The length of the time interval must be a multiple of the timestep.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    void ComputeExceptVoltage(double tStart, double tEnd)
    {
        // This method iterates over timesteps, calling ComputeExceptVoltage(t) at
        // each one, to update all state variables except for V, using backward Euler.

        // Check length of time interval
        unsigned n_steps = (unsigned)((tEnd - tStart) / mDt + 0.5);
        assert(fabs(tStart + n_steps*mDt - tEnd) < 1e-12);

        // Loop over time
        double curr_time;
        for (unsigned i=0; i<n_steps; i++)
        {
            curr_time = tStart + i*mDt;

            // Compute other state variables
            ComputeOneStepExceptVoltage(curr_time);

#ifndef NDEBUG
            // Check gating variables are still in range
            VerifyStateVariables();
#endif // NDEBUG
        }
    }

    /**
     * Simulate this cell's behaviour between the time interval [tStart, tEnd],
     * with timestemp #mDt, updating the internal state variable values.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    void SolveAndUpdateState(double tStart, double tEnd)
    {
        TimeStepper stepper(tStart, tEnd, mDt);

        while (!stepper.IsTimeAtEnd())
        {
            double time = stepper.GetTime();

            // Compute next value of V
            UpdateTransmembranePotential(time);

            // Compute other state variables
            ComputeOneStepExceptVoltage(time);

            // Check gating variables are still in range
            VerifyStateVariables();

            stepper.AdvanceOneTimeStep();
        }
    }

private:
// LCOV_EXCL_START
    /**
     * This function should never be called - the cell class incorporates its own solver.
     *
     * @param time
     * @param rY
     * @param rDY
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
    {
        NEVER_REACHED;
    }
// LCOV_EXCL_STOP

protected:
    /**
     * Compute the values of all state variables, except the voltage, using backward Euler,
     * for one timestep from tStart.
     *
     * \note This method must be provided by subclasses.
     *
     * @param tStart  start of this timestep
     */
    virtual void ComputeOneStepExceptVoltage(double tStart)=0;

    /**
     * Perform a forward Euler step to update the transmembrane potential.
     *
     * \note This method must be provided by subclasses.
     *
     * @param time  start of this timestep
     */
    virtual void UpdateTransmembranePotential(double time)=0;
};


// Debugging
#ifndef NDEBUG
#include "OutputFileHandler.hpp"
#include <boost/foreach.hpp>
template<unsigned SIZE>
void DumpJacobianToFile(double time, const double rCurrentGuess[SIZE], double rJacobian[SIZE][SIZE],
                        const std::vector<double>& rY)
{
    OutputFileHandler handler("DumpJacs", false);
    out_stream p_file = handler.OpenOutputFile("J.txt", std::ios::app);
    (*p_file) << "At " << time << " " << SIZE << std::endl;
    (*p_file) << "rY";
    BOOST_FOREACH(double y, rY)
    {
        (*p_file) << " " << y;
    }
    (*p_file) << std::endl;
    (*p_file) << "rCurrentGuess";
    for (unsigned i=0; i<SIZE; i++)
    {
        (*p_file) << " " << rCurrentGuess[i];
    }
    (*p_file) << std::endl;
    (*p_file) << "rJacobian";
    for (unsigned i=0; i<SIZE; i++)
    {
        for (unsigned j=0; j<SIZE; j++)
        {
            (*p_file) << " " << rJacobian[i][j];
        }
    }
    (*p_file) << std::endl;
}
#endif

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractBackwardEulerCardiacCell)

#endif /*ABSTRACTBACKWARDEULERCARDIACCELL_HPP_*/
