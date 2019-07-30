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

#ifndef _ABSTRACTSTEADYSTATERUNNER_HPP_
#define _ABSTRACTSTEADYSTATERUNNER_HPP_

#ifdef CHASTE_CVODE

#include <boost/shared_ptr.hpp>
#include "AbstractCvodeCell.hpp"
#include "Exception.hpp"
#include "RegularStimulus.hpp"
#include "Warnings.hpp"

/**
 * This class is abstract to investigate different ways of getting a cell model to steady state.
 *
 * A simple concrete example is SteadyStateRunner.
 */
class AbstractSteadyStateRunner
{
private:
    /**
     * Default constructor is private to ensure it is never used.
     */
    AbstractSteadyStateRunner(){};

protected:
    /** The cell model to run to steady state */
    boost::shared_ptr<AbstractCvodeCell> mpModel;

    /** The number of evaluations that it took */
    unsigned mNumEvaluations;

    /** The maximum number of paces that are allowed */
    unsigned mMaxNumPaces;

    /** Whether to suppress output from this class (defaults to false) */
    bool mSuppressOutput;

    /**
     * Run the cell model to steady state
     * (defined by a small change in state variables between beats)
     *
     * Must be overridden in subclasses for different approximations.
     *
     * Subclasses must check they don't do more than mMaxNumPaces
     * and should increment mNumEvaluations to reflect how many paces
     * were simulated
     */
    virtual void RunToSteadyStateImplementation() = 0;

public:
    /**
     * Constructor
     *
     * @param pModel  The cell model to run to steady state.
     */
    AbstractSteadyStateRunner(boost::shared_ptr<AbstractCvodeCell> pModel);

    /**
     * Destructor (empty)
     */
    virtual ~AbstractSteadyStateRunner(){};

    /**
     * Run the cell model to steady state
     * (defined by a small change in state variables between beats)
     *
     * Requires a method called RunToSteadyStateImplementation which must be
     * overridden in subclasses for different approximations.
     *
     * @return whether the model reached steady state within the maximum number of evaluations
     */
    bool RunToSteadyState();

    /**
     * Stop class, and subclasses printing output to std::cout
     *
     * @param  suppress  whether to suppress output.
     */
    void SuppressOutput(bool suppress = true);

    /**
     * @return The number of cell model 'paces/beats' that had to be evaluated.
     */
    unsigned GetNumEvaluations();

    /**
     * @param numPaces  The maximum number of paces to do on the way to steady state
     */
    void SetMaxNumPaces(unsigned numPaces);
};
#endif // CHASTE_CVODE

#endif // _ABSTRACTSTEADYSTATERUNNER_HPP_
