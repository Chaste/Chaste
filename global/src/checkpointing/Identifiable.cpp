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

#include "Identifiable.hpp"

#include <algorithm>
#include <typeinfo>

#include <boost/serialization/extended_type_info.hpp>
#include <boost/serialization/extended_type_info_typeid.hpp>
#include <boost/serialization/extended_type_info_no_rtti.hpp>
#include <boost/serialization/type_info_implementation.hpp>
#include "Warnings.hpp"

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
    while (i != identifier.npos)
    {
        identifier.replace(i, s_comma.length(), s_dash);
        i = identifier.find(s_comma, i);
    }

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
    const boost::serialization::extended_type_info* p_type_info =
            boost::serialization::type_info_implementation<Identifiable>::type::get_const_instance().get_derived_extended_type_info(*this);
    if(p_type_info!=nullptr)
    {
        id = p_type_info->get_key();
    }
    else
    {
        WARN_ONCE_ONLY("Unable to determine the class type. If you are using C++ serialization may not be set up correctly. \n If you are using Python this behaviour is expected.");
        id = "UnknownClass-ReflectionFailed";
    }
#else
    id = boost::serialization::type_info_implementation<Identifiable>::type::get_derived_extended_type_info(*this)->get_key();
#endif
    return TidyTemplatedExportIdentifier(id);
}
