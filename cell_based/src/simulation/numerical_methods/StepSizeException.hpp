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

#ifndef _STEPSIZEEXCEPTION_HPP_
#define _STEPSIZEEXCEPTION_HPP_

/**
 * Exception class to handle the adaptive time stepping in off-lattice simulations.
 *
 * We use this object since we are using it to signal to parent code that it needs
 * to adapt the time step, rather than being an error for presentation to the user
 * (for which we would use an Exception object).
 */
class StepSizeException : public std::exception
{
private:

    /** A suggested new time step. */
    double mSuggestedNewStep;

    /** Error message. */
    const std::string mMessage;

    /** Whether or not to terminate the simulation. */
    bool mIsTerminal;

public:

    /**
     * Construct an exception with a message string.
     *
     * @param suggestedNewStep a suggestion for an updated timestep
     * @param message the error message to display
     * @param isTerminal whether the error is terminal if true the the simulation stops
     */
    StepSizeException(double suggestedNewStep, const std::string message, bool isTerminal):
        std::exception(),
        mSuggestedNewStep(suggestedNewStep),
        mMessage(message),
        mIsTerminal(isTerminal)
    {
    }

    /**
     * Overidden what() method to return the message specific to this error handler.
     *
     * @return the exception message
     */
    virtual const char* what() const throw() {
        return mMessage.c_str();
    }

    /** Destructor. */
    ~StepSizeException() throw() {}

    /**
     * @return mSuggestedNewStep.
     */
    double GetSuggestedNewStep()
    {
        return mSuggestedNewStep;
    }

    /**
     * @return mIsTerminal.
     */
    bool IsTerminal()
    {
        return mIsTerminal;
    }
};

#endif // _STEPSIZEEXCEPTION_HPP_
