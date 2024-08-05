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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FarhadifarForce.hpp"

#include "FarhadifarForce2.cppwg.hpp"

namespace py = pybind11;
typedef FarhadifarForce<2 > FarhadifarForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class FarhadifarForce2_Overrides : public FarhadifarForce2{
    public:
    using FarhadifarForce2::FarhadifarForce;
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            FarhadifarForce2,
            AddForceContribution,
                    rCellPopulation);
    }
    double GetLineTensionParameter(::Node<2> * pNodeA, ::Node<2> * pNodeB, ::VertexBasedCellPopulation<2> & rVertexCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            FarhadifarForce2,
            GetLineTensionParameter,
                    pNodeA,
        pNodeB,
        rVertexCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            FarhadifarForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_FarhadifarForce2_class(py::module &m){
py::class_<FarhadifarForce2 , FarhadifarForce2_Overrides , boost::shared_ptr<FarhadifarForce2 >  , AbstractForce<2>  >(m, "FarhadifarForce2")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(FarhadifarForce2::*)(::AbstractCellPopulation<2> &)) &FarhadifarForce2::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "GetLineTensionParameter",
            (double(FarhadifarForce2::*)(::Node<2> *, ::Node<2> *, ::VertexBasedCellPopulation<2> &)) &FarhadifarForce2::GetLineTensionParameter,
            " " , py::arg("pNodeA"), py::arg("pNodeB"), py::arg("rVertexCellPopulation") )
        .def(
            "GetAreaElasticityParameter",
            (double(FarhadifarForce2::*)()) &FarhadifarForce2::GetAreaElasticityParameter,
            " "  )
        .def(
            "GetPerimeterContractilityParameter",
            (double(FarhadifarForce2::*)()) &FarhadifarForce2::GetPerimeterContractilityParameter,
            " "  )
        .def(
            "GetLineTensionParameter",
            (double(FarhadifarForce2::*)()) &FarhadifarForce2::GetLineTensionParameter,
            " "  )
        .def(
            "GetBoundaryLineTensionParameter",
            (double(FarhadifarForce2::*)()) &FarhadifarForce2::GetBoundaryLineTensionParameter,
            " "  )
        .def(
            "GetTargetAreaParameter",
            (double(FarhadifarForce2::*)()) &FarhadifarForce2::GetTargetAreaParameter,
            " "  )
        .def(
            "SetAreaElasticityParameter",
            (void(FarhadifarForce2::*)(double)) &FarhadifarForce2::SetAreaElasticityParameter,
            " " , py::arg("areaElasticityParameter") )
        .def(
            "SetPerimeterContractilityParameter",
            (void(FarhadifarForce2::*)(double)) &FarhadifarForce2::SetPerimeterContractilityParameter,
            " " , py::arg("perimeterContractilityParameter") )
        .def(
            "SetLineTensionParameter",
            (void(FarhadifarForce2::*)(double)) &FarhadifarForce2::SetLineTensionParameter,
            " " , py::arg("lineTensionParameter") )
        .def(
            "SetBoundaryLineTensionParameter",
            (void(FarhadifarForce2::*)(double)) &FarhadifarForce2::SetBoundaryLineTensionParameter,
            " " , py::arg("boundaryLineTensionParameter") )
        .def(
            "SetTargetAreaParameter",
            (void(FarhadifarForce2::*)(double)) &FarhadifarForce2::SetTargetAreaParameter,
            " " , py::arg("targetAreaParameter") )
        .def(
            "OutputForceParameters",
            (void(FarhadifarForce2::*)(::out_stream &)) &FarhadifarForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
