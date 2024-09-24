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
#include "BuskeCompressionForce.hpp"

#include "BuskeCompressionForce_3.cppwg.hpp"

namespace py = pybind11;
typedef BuskeCompressionForce<3> BuskeCompressionForce_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class BuskeCompressionForce_3_Overrides : public BuskeCompressionForce_3
{
public:
    using BuskeCompressionForce_3::BuskeCompressionForce;
    void AddForceContribution(::AbstractCellPopulation<3> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            void,
            BuskeCompressionForce_3,
            AddForceContribution,
            rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            BuskeCompressionForce_3,
            OutputForceParameters,
            rParamsFile);
    }
};

void register_BuskeCompressionForce_3_class(py::module &m)
{
    py::class_<BuskeCompressionForce_3, BuskeCompressionForce_3_Overrides, boost::shared_ptr<BuskeCompressionForce_3>, AbstractForce<3>>(m, "BuskeCompressionForce_3")
        .def(py::init<>())
        .def("GetCompressionEnergyParameter",
            (double(BuskeCompressionForce_3::*)()) &BuskeCompressionForce_3::GetCompressionEnergyParameter,
            " ")
        .def("SetCompressionEnergyParameter",
            (void(BuskeCompressionForce_3::*)(double)) &BuskeCompressionForce_3::SetCompressionEnergyParameter,
            " ", py::arg("compressionEnergyParameter"))
        .def("AddForceContribution",
            (void(BuskeCompressionForce_3::*)(::AbstractCellPopulation<3> &)) &BuskeCompressionForce_3::AddForceContribution,
            " ", py::arg("rCellPopulation"))
        .def("OutputForceParameters",
            (void(BuskeCompressionForce_3::*)(::out_stream &)) &BuskeCompressionForce_3::OutputForceParameters,
            " ", py::arg("rParamsFile"))
    ;
}
