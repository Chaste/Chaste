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
#include "CellSrnModel.hpp"

#include "CellSrnModel.cppwg.hpp"

namespace py = pybind11;
typedef CellSrnModel CellSrnModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractSrnModel * _AbstractSrnModelPtr;

class CellSrnModel_Overrides : public CellSrnModel
{
public:
    using CellSrnModel::CellSrnModel;
    void Initialise() override
    {
        PYBIND11_OVERRIDE(
            void,
            CellSrnModel,
            Initialise,
            );
    }
    void ResetForDivision() override
    {
        PYBIND11_OVERRIDE(
            void,
            CellSrnModel,
            ResetForDivision,
            );
    }
    void SimulateToCurrentTime() override
    {
        PYBIND11_OVERRIDE(
            void,
            CellSrnModel,
            SimulateToCurrentTime,
            );
    }
    ::AbstractSrnModel * CreateSrnModel() override
    {
        PYBIND11_OVERRIDE(
            _AbstractSrnModelPtr,
            CellSrnModel,
            CreateSrnModel,
            );
    }
    void SetCell(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            void,
            CellSrnModel,
            SetCell,
            pCell);
    }
};

void register_CellSrnModel_class(py::module &m)
{
    py::class_<CellSrnModel, CellSrnModel_Overrides, boost::shared_ptr<CellSrnModel>, AbstractSrnModel>(m, "CellSrnModel")
        .def(py::init<>())
        .def("begin",
            (::CellSrnModel::iterator(CellSrnModel::*)()) &CellSrnModel::begin,
            " ")
        .def("end",
            (::CellSrnModel::iterator(CellSrnModel::*)()) &CellSrnModel::end,
            " ")
        .def("begin",
            (::CellSrnModel::const_iterator(CellSrnModel::*)() const) &CellSrnModel::begin,
            " ")
        .def("end",
            (::CellSrnModel::const_iterator(CellSrnModel::*)() const) &CellSrnModel::end,
            " ")
        .def("cbegin",
            (::CellSrnModel::const_iterator(CellSrnModel::*)() const) &CellSrnModel::cbegin,
            " ")
        .def("cend",
            (::CellSrnModel::const_iterator(CellSrnModel::*)() const) &CellSrnModel::cend,
            " ")
        .def("Initialise",
            (void(CellSrnModel::*)()) &CellSrnModel::Initialise,
            " ")
        .def("ResetForDivision",
            (void(CellSrnModel::*)()) &CellSrnModel::ResetForDivision,
            " ")
        .def("SimulateToCurrentTime",
            (void(CellSrnModel::*)()) &CellSrnModel::SimulateToCurrentTime,
            " ")
        .def("CreateSrnModel",
            (::AbstractSrnModel *(CellSrnModel::*)()) &CellSrnModel::CreateSrnModel,
            " ", py::return_value_policy::reference)
        .def("AddEdgeSrn",
            (void(CellSrnModel::*)(::std::vector<boost::shared_ptr<AbstractSrnModel>>)) &CellSrnModel::AddEdgeSrn,
            " ", py::arg("edgeSrns"))
        .def("AddEdgeSrnModel",
            (void(CellSrnModel::*)(::AbstractSrnModelPtr)) &CellSrnModel::AddEdgeSrnModel,
            " ", py::arg("pEdgeSrn"))
        .def("GetNumEdgeSrn",
            (unsigned int(CellSrnModel::*)() const) &CellSrnModel::GetNumEdgeSrn,
            " ")
        .def("GetEdgeSrn",
            (::AbstractSrnModelPtr(CellSrnModel::*)(unsigned int) const) &CellSrnModel::GetEdgeSrn,
            " ", py::arg("index"))
        .def("GetEdges",
            (::std::vector<boost::shared_ptr<AbstractSrnModel>> const &(CellSrnModel::*)() const) &CellSrnModel::GetEdges,
            " ", py::return_value_policy::reference_internal)
        .def("SetInteriorSrnModel",
            (void(CellSrnModel::*)(::AbstractSrnModelPtr)) &CellSrnModel::SetInteriorSrnModel,
            " ", py::arg("pInteriorSrn"))
        .def("GetInteriorSrn",
            (::AbstractSrnModelPtr(CellSrnModel::*)() const) &CellSrnModel::GetInteriorSrn,
            " ")
        .def("SetCell",
            (void(CellSrnModel::*)(::CellPtr)) &CellSrnModel::SetCell,
            " ", py::arg("pCell"))
    ;
}
