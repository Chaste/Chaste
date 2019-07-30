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

#ifndef _ELEMENTATTRIBUTES_HPP_
#define _ELEMENTATTRIBUTES_HPP_

#include "UblasVectorInclude.hpp"
#include "ChasteSerialization.hpp"
#include "Exception.hpp"
#include <boost/serialization/vector.hpp>

/**
 * A container for attributes associated with the element classes.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ElementAttributes
{
private:

    /** Arbitrary attributes that a user gives meaning to */
    std::vector<double> mAttributes;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mAttributes;
    }

public:

    /**
     * Defaults all variables.
     */
    ElementAttributes();

    /**
     * @return mAttributes
     */
    std::vector<double>& rGetAttributes();

    /**
     * Push an attribute back onto mAttributes
     *
     * @param attribute the value of the attribute.
     */
    void AddAttribute(double attribute);

    /**
     * Either push the first attribute onto the vector or replace its value if it already exists.
     * This method gives back-compatibility to use-cases where there was only one attribute
     * @param attribute the value of the first attribute.
     */
    void SetFirstAttribute(double attribute);

    /**
     * @return the first (zero-indexed) attribute value or zero if none exist
     * \todo #2739 This should throw an exception if there are no attributes
     */
    double GetFirstAttribute();
};
#endif //_ELEMENTATTRIBUTES_HPP_
