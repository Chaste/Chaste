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
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"

#include "ImmersedBoundaryPalisadeMeshGenerator.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryPalisadeMeshGenerator ImmersedBoundaryPalisadeMeshGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_ImmersedBoundaryPalisadeMeshGenerator_class(py::module &m)
{
    py::class_<ImmersedBoundaryPalisadeMeshGenerator, boost::shared_ptr<ImmersedBoundaryPalisadeMeshGenerator>>(m, "ImmersedBoundaryPalisadeMeshGenerator")
        .def(py::init<unsigned int, unsigned int, double, double, double, bool, bool, bool, unsigned int, double>(), py::arg("numCellsWide"), py::arg("numNodesPerCell") = 100, py::arg("ellipseExponent") = 0.20000000000000001, py::arg("cellAspectRatio") = 2., py::arg("randomYMult") = 0., py::arg("basalLamina") = false, py::arg("apicalLamina") = false, py::arg("leakyLaminas") = false, py::arg("numFluidMeshPoints") = (2147483647 * 2U + 1U), py::arg("absoluteGap") = DOUBLE_UNSET)
        .def(py::init<>())
        .def("GetMesh",
            (::ImmersedBoundaryMesh<2, 2> *(ImmersedBoundaryPalisadeMeshGenerator::*)()) &ImmersedBoundaryPalisadeMeshGenerator::GetMesh,
            " ", py::return_value_policy::reference)
    ;
}
