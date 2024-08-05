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
#include "DeltaNotchInteriorSrnModel.hpp"

#include "DeltaNotchInteriorSrnModel.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchInteriorSrnModel DeltaNotchInteriorSrnModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractSrnModel * _AbstractSrnModelPtr;

class DeltaNotchInteriorSrnModel_Overrides : public DeltaNotchInteriorSrnModel{
    public:
    using DeltaNotchInteriorSrnModel::DeltaNotchInteriorSrnModel;
    ::AbstractSrnModel * CreateSrnModel() override {
        PYBIND11_OVERRIDE(
            _AbstractSrnModelPtr,
            DeltaNotchInteriorSrnModel,
            CreateSrnModel,
            );
    }
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchInteriorSrnModel,
            ResetForDivision,
            );
    }
    void Initialise() override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchInteriorSrnModel,
            Initialise,
            );
    }
    void SimulateToCurrentTime() override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchInteriorSrnModel,
            SimulateToCurrentTime,
            );
    }
    void OutputSrnModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchInteriorSrnModel,
            OutputSrnModelParameters,
                    rParamsFile);
    }
    void AddShrunkEdgeToInterior(::AbstractSrnModel * pShrunkEdgeSrn) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchInteriorSrnModel,
            AddShrunkEdgeToInterior,
                    pShrunkEdgeSrn);
    }

};
void register_DeltaNotchInteriorSrnModel_class(py::module &m){
py::class_<DeltaNotchInteriorSrnModel , DeltaNotchInteriorSrnModel_Overrides , boost::shared_ptr<DeltaNotchInteriorSrnModel >  , AbstractOdeSrnModel  >(m, "DeltaNotchInteriorSrnModel")
        .def(py::init<::boost::shared_ptr<AbstractCellCycleModelOdeSolver> >(), py::arg("pOdeSolver") = boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
        .def(
            "CreateSrnModel",
            (::AbstractSrnModel *(DeltaNotchInteriorSrnModel::*)()) &DeltaNotchInteriorSrnModel::CreateSrnModel,
            " "  , py::return_value_policy::reference)
        .def(
            "ResetForDivision",
            (void(DeltaNotchInteriorSrnModel::*)()) &DeltaNotchInteriorSrnModel::ResetForDivision,
            " "  )
        .def(
            "Initialise",
            (void(DeltaNotchInteriorSrnModel::*)()) &DeltaNotchInteriorSrnModel::Initialise,
            " "  )
        .def(
            "SimulateToCurrentTime",
            (void(DeltaNotchInteriorSrnModel::*)()) &DeltaNotchInteriorSrnModel::SimulateToCurrentTime,
            " "  )
        .def(
            "UpdateDeltaNotch",
            (void(DeltaNotchInteriorSrnModel::*)()) &DeltaNotchInteriorSrnModel::UpdateDeltaNotch,
            " "  )
        .def(
            "GetNotch",
            (double(DeltaNotchInteriorSrnModel::*)()) &DeltaNotchInteriorSrnModel::GetNotch,
            " "  )
        .def(
            "SetNotch",
            (void(DeltaNotchInteriorSrnModel::*)(double)) &DeltaNotchInteriorSrnModel::SetNotch,
            " " , py::arg("value") )
        .def(
            "GetDelta",
            (double(DeltaNotchInteriorSrnModel::*)()) &DeltaNotchInteriorSrnModel::GetDelta,
            " "  )
        .def(
            "SetDelta",
            (void(DeltaNotchInteriorSrnModel::*)(double)) &DeltaNotchInteriorSrnModel::SetDelta,
            " " , py::arg("value") )
        .def(
            "GetTotalEdgeDelta",
            (double(DeltaNotchInteriorSrnModel::*)()) &DeltaNotchInteriorSrnModel::GetTotalEdgeDelta,
            " "  )
        .def(
            "GetTotalEdgeNotch",
            (double(DeltaNotchInteriorSrnModel::*)()) &DeltaNotchInteriorSrnModel::GetTotalEdgeNotch,
            " "  )
        .def(
            "OutputSrnModelParameters",
            (void(DeltaNotchInteriorSrnModel::*)(::out_stream &)) &DeltaNotchInteriorSrnModel::OutputSrnModelParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "AddShrunkEdgeToInterior",
            (void(DeltaNotchInteriorSrnModel::*)(::AbstractSrnModel *)) &DeltaNotchInteriorSrnModel::AddShrunkEdgeToInterior,
            " " , py::arg("pShrunkEdgeSrn") )
    ;
}
