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

#ifndef _ABSTRACTSTIMULUSFUNCTION_HPP_
#define _ABSTRACTSTIMULUSFUNCTION_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include <cfloat>

#include "Exception.hpp"


/**
 * Represents an abstract stimulus function. Sub-classes will implement the
 * GetStimulus() function to represent the various type of stimuli to the cardiac
 * cell.
 */
class AbstractStimulusFunction
{
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
        // No member variables, but this is here so boost is happy serializing these classes.
    }
public:

    /**
     * @return the stimulus at a given time.
     *
     * @param time  time at which to return the stimulus
     */
    virtual double GetStimulus(double time) = 0;

    /**
     * Destructor.
     */
    virtual ~AbstractStimulusFunction();

    /**
     * Clear is used to managed memory in subclasses where the destructor may or may not need to clean up.
     */
    virtual void Clear();
};

CLASS_IS_ABSTRACT(AbstractStimulusFunction)

#endif //_ABSTRACTSTIMULUSFUNCTION_HPP_

