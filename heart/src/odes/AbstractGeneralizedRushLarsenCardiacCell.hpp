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

#ifndef ABSTRACTGENERALIZEDRUSHLARSENCARDIACCELL_HPP_
#define ABSTRACTGENERALIZEDRUSHLARSENCARDIACCELL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "ClassIsAbstract.hpp"
#include "AbstractCardiacCell.hpp"
#include "PetscTools.hpp"

/*
Megan E. Marsh, Raymond J. Spiteri
Numerical Simulation Laboratory
University of Saskatchewan
December 2011
Partial support provided by research grants from the National
Science and Engineering Research Council (NSERC) of Canada
and the MITACS/Mprime Canadian Network of Centres of Excellence.
*/

/**
 * This is the base class for cardiac cells solved using the GRL methods (GRL1 and GRL2).
 * Modified from AbstractRushLarsenCardiacCell.hpp
 */
class AbstractGeneralizedRushLarsenCardiacCell : public AbstractCardiacCell
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
    AbstractGeneralizedRushLarsenCardiacCell(
            unsigned numberOfStateVariables,
            unsigned voltageIndex,
            boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /** Virtual destructor */
    virtual ~AbstractGeneralizedRushLarsenCardiacCell();

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt.
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
     * with timestep #mDt.  The transmembrane potential is kept fixed throughout,
     * but the other state variables are updated.
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


    /**
     * @return whether the ODE system has an analytic Jacobian (#mHasAnalyticJacobian).
     */
    bool HasAnalyticJacobian() const;

    /**
     * Force the use of a numerical Jacobian, even if an analytic form is provided.
     * This is needed for a handful of troublesome models.
     *
     * @param useNumericalJacobian  Whether to use a numerical instead of the analytic Jacobian.
     */
    void ForceUseOfNumericalJacobian(bool useNumericalJacobian = true);


private:
// LCOV_EXCL_START
    /**
     * This function should never be called - the cell class incorporates its own solver.
     *
     * @param time
     * @param rY
     * @param rDY
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        NEVER_REACHED;
    }
// LCOV_EXCL_STOP

protected:
    /**
     * Update the values of all variables except the transmembrane potential using a GRL method.
     *
     * @param time  the current simulation time
     */
    virtual void ComputeOneStepExceptVoltage(double time)=0;

    /**
     * Perform a forward Euler step to update the transmembrane potential.
     *
     * @param time  the current simulation time
     */
    virtual void UpdateTransmembranePotential(double time)=0;

    /** The diagonal of the Jacobian, working memory for use by subclasses. */
    std::vector<double> mPartialF;

    /** The derivatives, working memory for use by subclasses. */
    std::vector<double> mEvalF;

    /** The state at the beginning of the current step, working memory for use by subclasses. */
    std::vector<double> mYInit;

    /** Whether we have an analytic Jacobian. */
    bool mHasAnalyticJacobian;
};

CLASS_IS_ABSTRACT(AbstractGeneralizedRushLarsenCardiacCell)

#endif // ABSTRACTGENERALIZEDRUSHLARSENCARDIACCELL_HPP_
