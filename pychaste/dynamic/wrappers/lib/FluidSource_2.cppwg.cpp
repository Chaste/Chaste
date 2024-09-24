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
#include "FluidSource.hpp"

#include "FluidSource_2.cppwg.hpp"

namespace py = pybind11;
typedef FluidSource<2> FluidSource_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_FluidSource_2_class(py::module &m)
{
    py::class_<FluidSource_2, boost::shared_ptr<FluidSource_2>>(m, "FluidSource_2")
        .def(py::init<unsigned int, ::ChastePoint<2>>(), py::arg("index"), py::arg("point"))
        .def(py::init<unsigned int, ::boost::numeric::ublas::c_vector<double, 2>>(), py::arg("index"), py::arg("location"))
        .def(py::init<unsigned int, double, double, double>(), py::arg("index"), py::arg("v1") = 0., py::arg("v2") = 0., py::arg("v3") = 0.)
        .def("GetIndex",
            (unsigned int(FluidSource_2::*)() const) &FluidSource_2::GetIndex,
            " ")
        .def("SetIndex",
            (void(FluidSource_2::*)(unsigned int)) &FluidSource_2::SetIndex,
            " ", py::arg("index"))
        .def("GetPoint",
            (::ChastePoint<2>(FluidSource_2::*)() const) &FluidSource_2::GetPoint,
            " ")
        .def("rGetLocation",
            (::boost::numeric::ublas::c_vector<double, 2> const &(FluidSource_2::*)() const) &FluidSource_2::rGetLocation,
            " ", py::return_value_policy::reference_internal)
        .def("rGetModifiableLocation",
            (::boost::numeric::ublas::c_vector<double, 2> &(FluidSource_2::*)()) &FluidSource_2::rGetModifiableLocation,
            " ", py::return_value_policy::reference_internal)
        .def("GetStrength",
            (double(FluidSource_2::*)() const) &FluidSource_2::GetStrength,
            " ")
        .def("SetStrength",
            (void(FluidSource_2::*)(double)) &FluidSource_2::SetStrength,
            " ", py::arg("strength"))
        .def("SetIfSourceIsAssociatedWithElement",
            (void(FluidSource_2::*)(bool)) &FluidSource_2::SetIfSourceIsAssociatedWithElement,
            " ", py::arg("associated"))
        .def("IsSourceAssociatedWithElement",
            (bool(FluidSource_2::*)()) &FluidSource_2::IsSourceAssociatedWithElement,
            " ")
        .def("GetAssociatedElementIndex",
            (unsigned int(FluidSource_2::*)() const) &FluidSource_2::GetAssociatedElementIndex,
            " ")
        .def("SetAssociatedElementIndex",
            (void(FluidSource_2::*)(unsigned int)) &FluidSource_2::SetAssociatedElementIndex,
            " ", py::arg("associatedElementIndex"))
    ;
}
