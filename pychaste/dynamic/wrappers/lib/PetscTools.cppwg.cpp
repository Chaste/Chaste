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
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include "PythonPetscObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PythonUblasObjectConverters.hpp"
#include "PetscTools.hpp"

#include "PetscTools.cppwg.hpp"

namespace py = pybind11;
typedef PetscTools PetscTools;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_PetscTools_class(py::module &m)
{
    py::class_<PetscTools, boost::shared_ptr<PetscTools>>(m, "PetscTools")
        .def(py::init<>())
        .def_static("ResetCache",
            (void(*)()) &PetscTools::ResetCache,
            " ")
        .def_static("IsInitialised",
            (bool(*)()) &PetscTools::IsInitialised,
            " ")
        .def_static("IsSequential",
            (bool(*)()) &PetscTools::IsSequential,
            " ")
        .def_static("IsParallel",
            (bool(*)()) &PetscTools::IsParallel,
            " ")
        .def_static("IsIsolated",
            (bool(*)()) &PetscTools::IsIsolated,
            " ")
        .def_static("GetNumProcs",
            (unsigned int(*)()) &PetscTools::GetNumProcs,
            " ")
        .def_static("GetMyRank",
            (unsigned int(*)()) &PetscTools::GetMyRank,
            " ")
        .def_static("AmMaster",
            (bool(*)()) &PetscTools::AmMaster,
            " ")
        .def_static("AmTopMost",
            (bool(*)()) &PetscTools::AmTopMost,
            " ")
        .def_static("Barrier",
            (void(*)(::std::string const)) &PetscTools::Barrier,
            " ", py::arg("callerId") = "")
        .def_static("BeginRoundRobin",
            (void(*)()) &PetscTools::BeginRoundRobin,
            " ")
        .def_static("EndRoundRobin",
            (void(*)()) &PetscTools::EndRoundRobin,
            " ")
        .def_static("IsolateProcesses",
            (void(*)(bool)) &PetscTools::IsolateProcesses,
            " ", py::arg("isolate") = true)
        .def_static("CreateVec",
            (::Vec(*)(int, int, bool)) &PetscTools::CreateVec,
            " ", py::arg("size"), py::arg("localSize") = -1, py::arg("ignoreOffProcEntries") = true, py::return_value_policy::reference)
        .def_static("CreateVec",
            (::Vec(*)(::std::vector<double>)) &PetscTools::CreateVec,
            " ", py::arg("data"), py::return_value_policy::reference)
        .def_static("CreateAndSetVec",
            (::Vec(*)(int, double)) &PetscTools::CreateAndSetVec,
            " ", py::arg("size"), py::arg("value"), py::return_value_policy::reference)
        .def_static("ReplicateBool",
            (bool(*)(bool)) &PetscTools::ReplicateBool,
            " ", py::arg("flag"))
        .def_static("ReplicateException",
            (void(*)(bool)) &PetscTools::ReplicateException,
            " ", py::arg("flag"))
        .def_static("DumpPetscObject",
            (void(*)(::Mat const &, ::std::string const &)) &PetscTools::DumpPetscObject,
            " ", py::arg("rMat"), py::arg("rOutputFileFullPath"))
        .def_static("DumpPetscObject",
            (void(*)(::Vec const &, ::std::string const &)) &PetscTools::DumpPetscObject,
            " ", py::arg("rVec"), py::arg("rOutputFileFullPath"))
        .def_static("HasParMetis",
            (bool(*)()) &PetscTools::HasParMetis,
            " ")
        .def_static("SetOption",
            (void(*)(char const *, char const *)) &PetscTools::SetOption,
            " ", py::arg("pOptionName"), py::arg("pOptionValue"))
        .def_static("ChasteMatCopy",
            (::PetscErrorCode(*)(::Mat, ::Mat, ::MatStructure)) &PetscTools::ChasteMatCopy,
            " ", py::arg("A"), py::arg("B"), py::arg("str"))
    ;
}
