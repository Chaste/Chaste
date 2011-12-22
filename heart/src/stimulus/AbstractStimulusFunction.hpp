/*

Copyright (C) University of Oxford, 2005-2011

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
     * Returns the stimulus at a given time.
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

