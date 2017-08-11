/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef ABSTRACTRUSHLARSEN21CARDIACCELL_HPP_
#define ABSTRACTRUSHLARSEN21CARDIACCELL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "ClassIsAbstract.hpp"

#include "AbstractCardiacCell.hpp"
#include "PetscTools.hpp"

/**

 * This is the base class for cardiac cells solved using the Rush-Larsen21 method.
 * written by Wenxian Guo (w.guo@usask.ca), University of Saskatchewan
 * with support from Dr. Raymond Spiteri (spiteri@cs.usask.ca) 
 * It is based on code contributed by Megan Lewis, University of Saskatchewan
 *
 * The basic idea is to split the ODE system with gates (gating variables) and non-gates
 * Gates are solved using exponential integrator as in Rush-Larsen method
 * Non-gates are solved using RKC21 (2-stage Runge-Kutta-Chebyshev method of order 1)
 * See ode/src/solver/RKC21IvpOdeSolver.hpp or Chaste ticket#2901 for implementation details
 *
 *  \li Compute alpha & beta values for gating variables, and derivatives (rDY1) for RKC stage 1
 *      (EvaluateEquations)
 *  \li Advance gating variabless using AdvanceGatingVars
 *  \li Advance non-gating variables using ComputeOneStepExceptVoltage
 *  \li Advance V using either UpdateTransmembranePotential or external equations
 *  
 *
 *  ABOVE WORKFLOW MIGHT NOT BE A GOOD DESIGN PATTERN
 *
 *
 */

class AbstractRushLarsen21CardiacCell : public AbstractCardiacCell
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
     *  \li It's a pity that inheriting from AbstractCardiacCell forces us to store a
     *      null pointer (for the unused ODE solver) in every instance.  We may want
     *      to revisit this design decision at a later date.
     */
    AbstractRushLarsen21CardiacCell(
        unsigned numberOfStateVariables,
        unsigned voltageIndex,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /** Virtual destructor */
    virtual ~AbstractRushLarsen21CardiacCell();

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt.  Uses a RKC21 step to update the transmembrane
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
     * with timestep #mDt.  The transmembrane potential is kept fixed throughout
     * (and is updated through external models), but the other state variables are 
     * updated (using RKC21 step or exponential integrator).
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
    * This function does the following:
    * 1. Update gating variables using exponential integrator
    * 2. Store new gating variables and old non-gating variables in rState (in this class)
    * 3. Time integration of non-gating variables (including V) using one step RKC21, get result
    *    at t = mu1t * dt
    * 4. Linear interpolation of gating variables, get result at t = mu1t * dt
    * 5. Override state vector in cell class with results from 3 and 4
    * 
    * This function is implemented in cell class
    *
    **/
    virtual void AdvanceGatingVars(const std::vector<double> &rDY1, // Advance non-gating variables
                                   std::vector<double> &rState, // To protect states for 2 step method
                                   const std::vector<double> &rAlphaOrTau,
                                   const std::vector<double> &rBetaOrInf)=0;

    /**
     * Update the values of non-gating variables using the RKC21 method for a single timestep.
     *
     * \note This method must be provided by subclasses.
     *
     * @param rState  vector containing state variables. Gating variables are at time t = dt,
     *              Non-gating variables are at time t = 0
     * @param rDY1  vector containing RHS at t = 0
     * @param rDY2  vector containing RHS at t = mu1t * dt
     */
    virtual void ComputeOneStepForNonGatingVarsExceptVoltage(const std::vector<double> &rState,
                                             const std::vector<double> &rDY1,
                                             const std::vector<double> &rDY2)=0;


     /**
     * Compute RHS and alpha and beta values.
     *
     * \note This method must be provided by subclasses.
     *
     * @param time  start of this timestep
     * @param rDY  vector to fill in with dy/dt values
     * @param rAlphaOrTau  vector to fill in with alpha or tau values, depending on the formulation
     * @param rBetaOrInf  vector to fill in with beta or inf values, depending on the formulation
     */
    virtual void EvaluateEquations(double time,
                                   std::vector<double> &rDY,
                                   std::vector<double> &rAlphaOrTau,
                                   std::vector<double> &rBetaOrInf)=0;

    /*
     * Transmembrane potential is updated using RKC21
     */
    void UpdateTransmembranePotential(const std::vector<double> &rDY1, const std::vector<double> &rDY2);





};

CLASS_IS_ABSTRACT(AbstractRushLarsen21CardiacCell)

#endif // ABSTRACTRUSHLARSEN21CARDIACCELL_HPP_
