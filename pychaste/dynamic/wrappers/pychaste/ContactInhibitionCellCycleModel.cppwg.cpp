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
#include "ContactInhibitionCellCycleModel.hpp"

#include "ContactInhibitionCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef ContactInhibitionCellCycleModel ContactInhibitionCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class ContactInhibitionCellCycleModel_Overrides : public ContactInhibitionCellCycleModel
{
public:
    using ContactInhibitionCellCycleModel::ContactInhibitionCellCycleModel;
    void UpdateCellCyclePhase() override
    {
        PYBIND11_OVERRIDE(
            void,
            ContactInhibitionCellCycleModel,
            UpdateCellCyclePhase,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override
    {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            ContactInhibitionCellCycleModel,
            CreateCellCycleModel,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            ContactInhibitionCellCycleModel,
            OutputCellCycleModelParameters,
            rParamsFile);
    }
};

void register_ContactInhibitionCellCycleModel_class(py::module &m)
{
    py::class_<ContactInhibitionCellCycleModel, ContactInhibitionCellCycleModel_Overrides, boost::shared_ptr<ContactInhibitionCellCycleModel>, AbstractSimplePhaseBasedCellCycleModel>(m, "ContactInhibitionCellCycleModel")
        .def(py::init<>())
        .def("UpdateCellCyclePhase",
            (void(ContactInhibitionCellCycleModel::*)()) &ContactInhibitionCellCycleModel::UpdateCellCyclePhase,
            " ")
        .def("CreateCellCycleModel",
            (::AbstractCellCycleModel *(ContactInhibitionCellCycleModel::*)()) &ContactInhibitionCellCycleModel::CreateCellCycleModel,
            " ", py::return_value_policy::reference)
        .def("SetQuiescentVolumeFraction",
            (void(ContactInhibitionCellCycleModel::*)(double)) &ContactInhibitionCellCycleModel::SetQuiescentVolumeFraction,
            " ", py::arg("quiescentVolumeFraction"))
        .def("GetQuiescentVolumeFraction",
            (double(ContactInhibitionCellCycleModel::*)() const) &ContactInhibitionCellCycleModel::GetQuiescentVolumeFraction,
            " ")
        .def("SetEquilibriumVolume",
            (void(ContactInhibitionCellCycleModel::*)(double)) &ContactInhibitionCellCycleModel::SetEquilibriumVolume,
            " ", py::arg("equilibriumVolume"))
        .def("GetEquilibriumVolume",
            (double(ContactInhibitionCellCycleModel::*)() const) &ContactInhibitionCellCycleModel::GetEquilibriumVolume,
            " ")
        .def("GetCurrentQuiescentDuration",
            (double(ContactInhibitionCellCycleModel::*)() const) &ContactInhibitionCellCycleModel::GetCurrentQuiescentDuration,
            " ")
        .def("GetCurrentQuiescentOnsetTime",
            (double(ContactInhibitionCellCycleModel::*)() const) &ContactInhibitionCellCycleModel::GetCurrentQuiescentOnsetTime,
            " ")
        .def("OutputCellCycleModelParameters",
            (void(ContactInhibitionCellCycleModel::*)(::out_stream &)) &ContactInhibitionCellCycleModel::OutputCellCycleModelParameters,
            " ", py::arg("rParamsFile"))
    ;
}
