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
#include "DeltaNotchEdgeSrnModel.hpp"

#include "DeltaNotchEdgeSrnModel.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchEdgeSrnModel DeltaNotchEdgeSrnModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractSrnModel * _AbstractSrnModelPtr;

class DeltaNotchEdgeSrnModel_Overrides : public DeltaNotchEdgeSrnModel{
    public:
    using DeltaNotchEdgeSrnModel::DeltaNotchEdgeSrnModel;
    ::AbstractSrnModel * CreateSrnModel() override {
        PYBIND11_OVERRIDE(
            _AbstractSrnModelPtr,
            DeltaNotchEdgeSrnModel,
            CreateSrnModel,
            );
    }
    void Initialise() override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeSrnModel,
            Initialise,
            );
    }
    void InitialiseDaughterCell() override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeSrnModel,
            InitialiseDaughterCell,
            );
    }
    void SimulateToCurrentTime() override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeSrnModel,
            SimulateToCurrentTime,
            );
    }
    void OutputSrnModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeSrnModel,
            OutputSrnModelParameters,
                    rParamsFile);
    }
    void AddSrnQuantities(::AbstractSrnModel * pOtherSrn, double const scale) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeSrnModel,
            AddSrnQuantities,
                    pOtherSrn,
        scale);
    }
    void AddShrunkEdgeSrn(::AbstractSrnModel * pShrunkEdgeSrn) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeSrnModel,
            AddShrunkEdgeSrn,
                    pShrunkEdgeSrn);
    }
    void AddMergedEdgeSrn(::AbstractSrnModel * pMergedEdgeSrn) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeSrnModel,
            AddMergedEdgeSrn,
                    pMergedEdgeSrn);
    }
    void SplitEdgeSrn(double const relativePosition) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeSrnModel,
            SplitEdgeSrn,
                    relativePosition);
    }

};
void register_DeltaNotchEdgeSrnModel_class(py::module &m){
py::class_<DeltaNotchEdgeSrnModel , DeltaNotchEdgeSrnModel_Overrides , boost::shared_ptr<DeltaNotchEdgeSrnModel >  , AbstractOdeSrnModel  >(m, "DeltaNotchEdgeSrnModel")
        .def(py::init<::boost::shared_ptr<AbstractCellCycleModelOdeSolver> >(), py::arg("pOdeSolver") = boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
        .def(
            "CreateSrnModel",
            (::AbstractSrnModel *(DeltaNotchEdgeSrnModel::*)()) &DeltaNotchEdgeSrnModel::CreateSrnModel,
            " "  , py::return_value_policy::reference)
        .def(
            "Initialise",
            (void(DeltaNotchEdgeSrnModel::*)()) &DeltaNotchEdgeSrnModel::Initialise,
            " "  )
        .def(
            "InitialiseDaughterCell",
            (void(DeltaNotchEdgeSrnModel::*)()) &DeltaNotchEdgeSrnModel::InitialiseDaughterCell,
            " "  )
        .def(
            "SimulateToCurrentTime",
            (void(DeltaNotchEdgeSrnModel::*)()) &DeltaNotchEdgeSrnModel::SimulateToCurrentTime,
            " "  )
        .def(
            "UpdateDeltaNotch",
            (void(DeltaNotchEdgeSrnModel::*)()) &DeltaNotchEdgeSrnModel::UpdateDeltaNotch,
            " "  )
        .def(
            "GetNotch",
            (double(DeltaNotchEdgeSrnModel::*)()) &DeltaNotchEdgeSrnModel::GetNotch,
            " "  )
        .def(
            "SetNotch",
            (void(DeltaNotchEdgeSrnModel::*)(double)) &DeltaNotchEdgeSrnModel::SetNotch,
            " " , py::arg("value") )
        .def(
            "GetDelta",
            (double(DeltaNotchEdgeSrnModel::*)()) &DeltaNotchEdgeSrnModel::GetDelta,
            " "  )
        .def(
            "SetDelta",
            (void(DeltaNotchEdgeSrnModel::*)(double)) &DeltaNotchEdgeSrnModel::SetDelta,
            " " , py::arg("value") )
        .def(
            "GetNeighbouringDelta",
            (double(DeltaNotchEdgeSrnModel::*)() const ) &DeltaNotchEdgeSrnModel::GetNeighbouringDelta,
            " "  )
        .def(
            "GetInteriorDelta",
            (double(DeltaNotchEdgeSrnModel::*)() const ) &DeltaNotchEdgeSrnModel::GetInteriorDelta,
            " "  )
        .def(
            "GetInteriorNotch",
            (double(DeltaNotchEdgeSrnModel::*)() const ) &DeltaNotchEdgeSrnModel::GetInteriorNotch,
            " "  )
        .def(
            "OutputSrnModelParameters",
            (void(DeltaNotchEdgeSrnModel::*)(::out_stream &)) &DeltaNotchEdgeSrnModel::OutputSrnModelParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "AddSrnQuantities",
            (void(DeltaNotchEdgeSrnModel::*)(::AbstractSrnModel *, double const)) &DeltaNotchEdgeSrnModel::AddSrnQuantities,
            " " , py::arg("pOtherSrn"), py::arg("scale") = 1. )
        .def(
            "AddShrunkEdgeSrn",
            (void(DeltaNotchEdgeSrnModel::*)(::AbstractSrnModel *)) &DeltaNotchEdgeSrnModel::AddShrunkEdgeSrn,
            " " , py::arg("pShrunkEdgeSrn") )
        .def(
            "AddMergedEdgeSrn",
            (void(DeltaNotchEdgeSrnModel::*)(::AbstractSrnModel *)) &DeltaNotchEdgeSrnModel::AddMergedEdgeSrn,
            " " , py::arg("pMergedEdgeSrn") )
        .def(
            "SplitEdgeSrn",
            (void(DeltaNotchEdgeSrnModel::*)(double const)) &DeltaNotchEdgeSrnModel::SplitEdgeSrn,
            " " , py::arg("relativePosition") )
    ;
}
