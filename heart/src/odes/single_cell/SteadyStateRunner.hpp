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

#ifndef _STEADYSTATERUNNER_HPP_
#define _STEADYSTATERUNNER_HPP_

#ifdef CHASTE_CVODE

#include "AbstractSteadyStateRunner.hpp"
#include "VectorHelperFunctions.hpp"

/**
 * This class is to get a cell model to (approximately) steady state.
 *
 * It does it by simply pacing and looking for a small change in the norm of the
 * state variables between subsequent paces.
 *
 * Note - because it looks for a small change in state variables between solves
 * this method requires quite an accurate CVODE solution to do this efficiently
 * (otherwise you just have to wait to get lucky on subsequent solves).
 *
 * So although running with stricter tolerances takes longer per pace
 * it can mean you detect a steady state in up to 3x fewer paces, so it is
 * recommended to use a model with SetTolerances(1e-6, 1e-8); or stricter.
 */
class SteadyStateRunner: public AbstractSteadyStateRunner
{
private:
    /** whether we should do two paces at once (should detect steady alternans as well as single paces) */
    bool mTwoPaceScan;

protected:
    /**
     * Run the cell model to steady state
     *
     * Here we don't do anything clever - we just gradually drift to the steady state,
     * defined by < 1e-6 change in the norm of state variables between 1 (or 2 - see #mTwoPaceScan) beats.
     */
    virtual void RunToSteadyStateImplementation();

public:
    /**
     * Constructor of a helper class for getting action potential models to steady state
     *
     * @param pModel  The cell model to run to steady state.
     * @param twoPaces  Whether to run two paces at once, for detection of steady state alternans.
     */
    SteadyStateRunner(boost::shared_ptr<AbstractCvodeCell> pModel, bool twoPaces=false)
     : AbstractSteadyStateRunner(pModel),
       mTwoPaceScan(twoPaces)
    {};
};

#endif // CHASTE_CVODE

#endif // _STEADYSTATERUNNER_HPP_

