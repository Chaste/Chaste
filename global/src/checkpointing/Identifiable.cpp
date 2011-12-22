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

#include "Identifiable.hpp"

#include <algorithm>
#include <typeinfo>

#include <boost/serialization/extended_type_info.hpp>
#include <boost/serialization/extended_type_info_typeid.hpp>
#include <boost/serialization/extended_type_info_no_rtti.hpp>
#include <boost/serialization/type_info_implementation.hpp>

std::string Identifiable::TidyTemplatedExportIdentifier(std::string identifier) const
{
    // First remove spaces, so identifier now takes the form "pack<void(NameOfDerivedType<DIM>)>::type"
    std::string::iterator end_pos = std::remove(identifier.begin(), identifier.end(), ' ');
    identifier.erase(end_pos, identifier.end());

    // Then remove "pack<void(", so identifier now takes the form "NameOfDerivedType<DIM>)>::type"
    const std::string s_pack = "pack<void(";
    std::string::size_type i = identifier.find(s_pack);
    if (i != identifier.npos)
    {
        identifier.erase(i, s_pack.length());
    }

    // Then replace "<" with "-", so identifier now takes the form "NameOfDerivedType-DIM>)>::type"
    const std::string s_open = "<";
    const std::string s_dash = "-";
    i = identifier.find(s_open);
    if (i != identifier.npos)
    {
        identifier.replace(i, s_open.length(), s_dash);
    }

    // Then replace "," with "-" to account for multiple template parameters
    const std::string s_comma = ",";
    i = identifier.find(s_comma);
    assert( i == identifier.npos ); /// #1453 There are currently no identifiables with multiple template parameters...

    /**\todo #1453 - implement the following if needed
    if (i != identifier.npos)
    {
        identifier.replace(i, s_comma.length(), s_dash);
    }
    */

    // Finally remove ">)>::type", so that identifier now takes the form "NameOfDerivedType-DIM"
    const std::string s_end = ">)>::type";
    i = identifier.find(s_end);
    if (i != identifier.npos)
    {
        identifier.erase(i, s_end.length());
    }

    return identifier;
}

Identifiable::~Identifiable()
{
}

std::string Identifiable::GetIdentifier() const
{
    std::string id;
#if BOOST_VERSION >= 103700
    id = boost::serialization::type_info_implementation<Identifiable>::type::get_const_instance().get_derived_extended_type_info(*this)->get_key();
#else
    id = boost::serialization::type_info_implementation<Identifiable>::type::get_derived_extended_type_info(*this)->get_key();
#endif
    return TidyTemplatedExportIdentifier(id);
}
