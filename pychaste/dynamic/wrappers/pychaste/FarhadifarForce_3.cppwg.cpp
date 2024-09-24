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

// This file is auto-generated. Manual changes will be overwritten. For changes
// to persist, update the configuration in pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FarhadifarForce.hpp"

#include "FarhadifarForce_3.cppwg.hpp"

namespace py = pybind11;
typedef FarhadifarForce<3> FarhadifarForce_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class FarhadifarForce_3_Overrides : public FarhadifarForce_3
{
public:
    using FarhadifarForce_3::FarhadifarForce;
    void AddForceContribution(::AbstractCellPopulation<3> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            void,
            FarhadifarForce_3,
            AddForceContribution,
            rCellPopulation);
    }
    double GetLineTensionParameter(::Node<3> * pNodeA, ::Node<3> * pNodeB, ::VertexBasedCellPopulation<3> & rVertexCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            double,
            FarhadifarForce_3,
            GetLineTensionParameter,
            pNodeA,
            pNodeB,
            rVertexCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            FarhadifarForce_3,
            OutputForceParameters,
            rParamsFile);
    }
};

void register_FarhadifarForce_3_class(py::module &m)
{
    py::class_<FarhadifarForce_3, FarhadifarForce_3_Overrides, boost::shared_ptr<FarhadifarForce_3>, AbstractForce<3>>(m, "FarhadifarForce_3")
        .def(py::init<>())
        .def("AddForceContribution",
            (void(FarhadifarForce_3::*)(::AbstractCellPopulation<3> &)) &FarhadifarForce_3::AddForceContribution,
            " ", py::arg("rCellPopulation"))
        .def("GetLineTensionParameter",
            (double(FarhadifarForce_3::*)(::Node<3> *, ::Node<3> *, ::VertexBasedCellPopulation<3> &)) &FarhadifarForce_3::GetLineTensionParameter,
            " ", py::arg("pNodeA"), py::arg("pNodeB"), py::arg("rVertexCellPopulation"))
        .def("GetAreaElasticityParameter",
            (double(FarhadifarForce_3::*)()) &FarhadifarForce_3::GetAreaElasticityParameter,
            " ")
        .def("GetPerimeterContractilityParameter",
            (double(FarhadifarForce_3::*)()) &FarhadifarForce_3::GetPerimeterContractilityParameter,
            " ")
        .def("GetLineTensionParameter",
            (double(FarhadifarForce_3::*)()) &FarhadifarForce_3::GetLineTensionParameter,
            " ")
        .def("GetBoundaryLineTensionParameter",
            (double(FarhadifarForce_3::*)()) &FarhadifarForce_3::GetBoundaryLineTensionParameter,
            " ")
        .def("GetTargetAreaParameter",
            (double(FarhadifarForce_3::*)()) &FarhadifarForce_3::GetTargetAreaParameter,
            " ")
        .def("SetAreaElasticityParameter",
            (void(FarhadifarForce_3::*)(double)) &FarhadifarForce_3::SetAreaElasticityParameter,
            " ", py::arg("areaElasticityParameter"))
        .def("SetPerimeterContractilityParameter",
            (void(FarhadifarForce_3::*)(double)) &FarhadifarForce_3::SetPerimeterContractilityParameter,
            " ", py::arg("perimeterContractilityParameter"))
        .def("SetLineTensionParameter",
            (void(FarhadifarForce_3::*)(double)) &FarhadifarForce_3::SetLineTensionParameter,
            " ", py::arg("lineTensionParameter"))
        .def("SetBoundaryLineTensionParameter",
            (void(FarhadifarForce_3::*)(double)) &FarhadifarForce_3::SetBoundaryLineTensionParameter,
            " ", py::arg("boundaryLineTensionParameter"))
        .def("SetTargetAreaParameter",
            (void(FarhadifarForce_3::*)(double)) &FarhadifarForce_3::SetTargetAreaParameter,
            " ", py::arg("targetAreaParameter"))
        .def("OutputForceParameters",
            (void(FarhadifarForce_3::*)(::out_stream &)) &FarhadifarForce_3::OutputForceParameters,
            " ", py::arg("rParamsFile"))
    ;
}
