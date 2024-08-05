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
#include "WelikyOsterForce.hpp"

#include "WelikyOsterForce3.cppwg.hpp"

namespace py = pybind11;
typedef WelikyOsterForce<3 > WelikyOsterForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class WelikyOsterForce3_Overrides : public WelikyOsterForce3{
    public:
    using WelikyOsterForce3::WelikyOsterForce;
    void AddForceContribution(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            WelikyOsterForce3,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            WelikyOsterForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_WelikyOsterForce3_class(py::module &m){
py::class_<WelikyOsterForce3 , WelikyOsterForce3_Overrides , boost::shared_ptr<WelikyOsterForce3 >  , AbstractForce<3>  >(m, "WelikyOsterForce3")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(WelikyOsterForce3::*)(::AbstractCellPopulation<3> &)) &WelikyOsterForce3::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "GetWelikyOsterAreaParameter",
            (double(WelikyOsterForce3::*)()) &WelikyOsterForce3::GetWelikyOsterAreaParameter,
            " "  )
        .def(
            "GetWelikyOsterPerimeterParameter",
            (double(WelikyOsterForce3::*)()) &WelikyOsterForce3::GetWelikyOsterPerimeterParameter,
            " "  )
        .def(
            "SetWelikyOsterAreaParameter",
            (void(WelikyOsterForce3::*)(double)) &WelikyOsterForce3::SetWelikyOsterAreaParameter,
            " " , py::arg("welikyOsterAreaParameter") )
        .def(
            "SetWelikyOsterPerimeterParameter",
            (void(WelikyOsterForce3::*)(double)) &WelikyOsterForce3::SetWelikyOsterPerimeterParameter,
            " " , py::arg("welikyOsterPerimeterParameter") )
        .def(
            "OutputForceParameters",
            (void(WelikyOsterForce3::*)(::out_stream &)) &WelikyOsterForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
