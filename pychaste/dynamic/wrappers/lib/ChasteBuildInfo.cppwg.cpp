/*

Copyright (c) 2005-2024, University of Oxford.
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

// This file is auto-generated; manual changes will be overwritten.
// To make enduring changes, see pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PythonUblasObjectConverters.hpp"
#include "Version.hpp"

#include "ChasteBuildInfo.cppwg.hpp"

namespace py = pybind11;
typedef ChasteBuildInfo ChasteBuildInfo;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_ChasteBuildInfo_class(py::module &m)
{
    py::class_<ChasteBuildInfo, boost::shared_ptr<ChasteBuildInfo>>(m, "ChasteBuildInfo")
        .def(py::init<>())
        .def_static("GetLicenceText",
            (::std::string(*)()) &ChasteBuildInfo::GetLicenceText,
            " ")
        .def_static("GetRootDir",
            (char const *(*)()) &ChasteBuildInfo::GetRootDir,
            " ", py::return_value_policy::reference)
        .def_static("GetVersionString",
            (::std::string(*)()) &ChasteBuildInfo::GetVersionString,
            " ")
        .def_static("GetMajorReleaseNumber",
            (unsigned int(*)()) &ChasteBuildInfo::GetMajorReleaseNumber,
            " ")
        .def_static("GetMinorReleaseNumber",
            (unsigned int(*)()) &ChasteBuildInfo::GetMinorReleaseNumber,
            " ")
        .def_static("GetRevisionNumber",
            (long long unsigned int(*)()) &ChasteBuildInfo::GetRevisionNumber,
            " ")
        .def_static("IsWorkingCopyModified",
            (bool(*)()) &ChasteBuildInfo::IsWorkingCopyModified,
            " ")
        .def_static("GetBuildTime",
            (char const *(*)()) &ChasteBuildInfo::GetBuildTime,
            " ", py::return_value_policy::reference)
        .def_static("GetCurrentTime",
            (char const *(*)()) &ChasteBuildInfo::GetCurrentTime,
            " ", py::return_value_policy::reference)
        .def_static("GetBuilderUnameInfo",
            (char const *(*)()) &ChasteBuildInfo::GetBuilderUnameInfo,
            " ", py::return_value_policy::reference)
        .def_static("GetBuildInformation",
            (char const *(*)()) &ChasteBuildInfo::GetBuildInformation,
            " ", py::return_value_policy::reference)
        .def_static("GetCompilerType",
            (char const *(*)()) &ChasteBuildInfo::GetCompilerType,
            " ", py::return_value_policy::reference)
        .def_static("GetCompilerVersion",
            (char const *(*)()) &ChasteBuildInfo::GetCompilerVersion,
            " ", py::return_value_policy::reference)
        .def_static("GetCompilerFlags",
            (char const *(*)()) &ChasteBuildInfo::GetCompilerFlags,
            " ", py::return_value_policy::reference)
        .def_static("GetXsdVersion",
            (char const *(*)()) &ChasteBuildInfo::GetXsdVersion,
            " ", py::return_value_policy::reference)
        .def_static("rGetProjectVersions",
            (::std::map<std::basic_string<char>, std::basic_string<char>> const &(*)()) &ChasteBuildInfo::rGetProjectVersions,
            " ", py::return_value_policy::reference_internal)
        .def_static("rGetIfProjectsModified",
            (::std::map<std::basic_string<char>, std::basic_string<char>> const &(*)()) &ChasteBuildInfo::rGetIfProjectsModified,
            " ", py::return_value_policy::reference_internal)
        .def_static("GetProvenanceString",
            (::std::string(*)()) &ChasteBuildInfo::GetProvenanceString,
            " ")
        .def_static("GetChasteCodegenVersion",
            (::std::string(*)()) &ChasteBuildInfo::GetChasteCodegenVersion,
            " ")
    ;
}
