/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef IDENTIFIABLE_HPP_
#define IDENTIFIABLE_HPP_

#include <string>
#include "ChasteSerialization.hpp"

#include <typeinfo>
#include <boost/serialization/extended_type_info.hpp>
#include <boost/serialization/extended_type_info_typeid.hpp>
#include <boost/serialization/extended_type_info_no_rtti.hpp>
#include <boost/serialization/type_info_implementation.hpp>

/**
 * "Mix-in" base class for any class that needs to provide a unique ID specifying
 * what type the derived class is.
 *
 * All you need to do to use this is inherit from this class at the base of
 * your hierarchy, and provide a Boost Serialization export key for every
 * concrete class.  Also, any test in which the GetIdentifier method is used,
 * even via the main cell_based code, \b must include CheckpointArchiveTypes.hpp
 * or CellBasedSimulationArchiver.hpp as the first Chaste header included.
 * Failure to do so will result in a seg fault.
 */
class Identifiable
{
public:

    /**
     * Virtual destructor to make this class polymorphic.
     */
    virtual ~Identifiable();

    /**
     * Return the unique identifier of the concrete class.
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
