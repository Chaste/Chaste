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

#ifndef IDENTIFIABLE_HPP_
#define IDENTIFIABLE_HPP_

#include <string>
#include "ChasteSerialization.hpp"

/**
 * "Mix-in" base class for any class that needs to provide a unique ID specifying
 * what type the derived class is.
 *
 * All you need to do to use this is inherit from this class at the base of
 * your hierarchy, and provide a Boost Serialization export key for every
 * concrete class.  Also, any test in which the GetIdentifier method is used,
 * even via the main cell_based code, \b must include CheckpointArchiveTypes.hpp
 * or CellBasedSimulationArchiver.hpp as the first Chaste header included.
 * Failure to do so will result in a segmentation fault.
 */
class Identifiable
{
public:

    /**
     * Virtual destructor to make this class polymorphic.
     */
    virtual ~Identifiable();

    /**
     * @return the unique identifier of the concrete class.
     *
     * This method uses Boost's serialization's extended_type_info and returns
     * the identifier of the derived class (this is defined when the macro
     * CHASTE_CLASS_EXPORT is invoked in each derived class, and is usually just
     * the name of the class).
     *
     * Note that you must include the header CheckpointArchiveTypes.hpp in any
     * test suite that calls this method.
     */
    std::string GetIdentifier() const;

private:

    /**
     * @return a name which is suitable for use as an XML element name
     * Templated classes get Boost Serialization export keys that look like
     * "pack<void (NameOfDerivedType< DIM >)>::type".
     * This method converts it to a nice name suitable for use as an XML element name,
     * i.e. of the form "NameOfDerivedType-DIM".
     * Works for classes templated over any number of parameters, providing the
     * values are allowable in XML element names.
     *
     * @param identifier  the identifier to tidy
     */
    std::string TidyTemplatedExportIdentifier(std::string identifier) const;
};

#endif // IDENTIFIABLE_HPP_
