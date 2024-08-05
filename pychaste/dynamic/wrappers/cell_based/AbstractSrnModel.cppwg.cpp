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
#include "AbstractSrnModel.hpp"

#include "AbstractSrnModel.cppwg.hpp"

namespace py = pybind11;
typedef AbstractSrnModel AbstractSrnModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractSrnModel * _AbstractSrnModelPtr;

class AbstractSrnModel_Overrides : public AbstractSrnModel{
    public:
    using AbstractSrnModel::AbstractSrnModel;
    void SetCell(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            SetCell,
                    pCell);
    }
    void Initialise() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            Initialise,
            );
    }
    void InitialiseDaughterCell() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            InitialiseDaughterCell,
            );
    }
    void SimulateToCurrentTime() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractSrnModel,
            SimulateToCurrentTime,
            );
    }
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            ResetForDivision,
            );
    }
    ::AbstractSrnModel * CreateSrnModel() override {
        PYBIND11_OVERRIDE_PURE(
            _AbstractSrnModelPtr,
            AbstractSrnModel,
            CreateSrnModel,
            );
    }
    void OutputSrnModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            OutputSrnModelParameters,
                    rParamsFile);
    }
    void ScaleSrnVariables(double const theta) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            ScaleSrnVariables,
                    theta);
    }
    void AddSrnQuantities(::AbstractSrnModel * pOtherSrn, double const scale) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            AddSrnQuantities,
                    pOtherSrn,
        scale);
    }
    void AddShrunkEdgeSrn(::AbstractSrnModel * pShrunkEdgeSrn) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            AddShrunkEdgeSrn,
                    pShrunkEdgeSrn);
    }
    void AddMergedEdgeSrn(::AbstractSrnModel * pMergedEdgeSrn) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            AddMergedEdgeSrn,
                    pMergedEdgeSrn);
    }
    void AddShrunkEdgeToInterior(::AbstractSrnModel * pShrunkEdgeSrn) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            AddShrunkEdgeToInterior,
                    pShrunkEdgeSrn);
    }
    void SplitEdgeSrn(double const relativePosition) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSrnModel,
            SplitEdgeSrn,
                    relativePosition);
    }

};
void register_AbstractSrnModel_class(py::module &m){
py::class_<AbstractSrnModel , AbstractSrnModel_Overrides , boost::shared_ptr<AbstractSrnModel >   >(m, "AbstractSrnModel")
        .def(py::init< >())
        .def(
            "SetCell",
            (void(AbstractSrnModel::*)(::CellPtr)) &AbstractSrnModel::SetCell,
            " " , py::arg("pCell") )
        .def(
            "Initialise",
            (void(AbstractSrnModel::*)()) &AbstractSrnModel::Initialise,
            " "  )
        .def(
            "InitialiseDaughterCell",
            (void(AbstractSrnModel::*)()) &AbstractSrnModel::InitialiseDaughterCell,
            " "  )
        .def(
            "GetCell",
            (::CellPtr(AbstractSrnModel::*)()) &AbstractSrnModel::GetCell,
            " "  )
        .def(
            "SetSimulatedToTime",
            (void(AbstractSrnModel::*)(double)) &AbstractSrnModel::SetSimulatedToTime,
            " " , py::arg("simulatedToTime") )
        .def(
            "GetSimulatedToTime",
            (double(AbstractSrnModel::*)() const ) &AbstractSrnModel::GetSimulatedToTime,
            " "  )
        .def(
            "SimulateToCurrentTime",
            (void(AbstractSrnModel::*)()) &AbstractSrnModel::SimulateToCurrentTime,
            " "  )
        .def(
            "ResetForDivision",
            (void(AbstractSrnModel::*)()) &AbstractSrnModel::ResetForDivision,
            " "  )
        .def(
            "CreateSrnModel",
            (::AbstractSrnModel *(AbstractSrnModel::*)()) &AbstractSrnModel::CreateSrnModel,
            " "  , py::return_value_policy::reference)
        .def(
            "OutputSrnModelInfo",
            (void(AbstractSrnModel::*)(::out_stream &)) &AbstractSrnModel::OutputSrnModelInfo,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputSrnModelParameters",
            (void(AbstractSrnModel::*)(::out_stream &)) &AbstractSrnModel::OutputSrnModelParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "SetEdgeLocalIndex",
            (void(AbstractSrnModel::*)(unsigned int)) &AbstractSrnModel::SetEdgeLocalIndex,
            " " , py::arg("index") )
        .def(
            "GetEdgeLocalIndex",
            (unsigned int(AbstractSrnModel::*)()) &AbstractSrnModel::GetEdgeLocalIndex,
            " "  )
        .def(
            "HasEdgeModel",
            (bool(AbstractSrnModel::*)() const ) &AbstractSrnModel::HasEdgeModel,
            " "  )
        .def(
            "SetEdgeModelIndicator",
            (void(AbstractSrnModel::*)(bool const)) &AbstractSrnModel::SetEdgeModelIndicator,
            " " , py::arg("isEdgeModel") )
        .def(
            "ScaleSrnVariables",
            (void(AbstractSrnModel::*)(double const)) &AbstractSrnModel::ScaleSrnVariables,
            " " , py::arg("theta") )
        .def(
            "AddSrnQuantities",
            (void(AbstractSrnModel::*)(::AbstractSrnModel *, double const)) &AbstractSrnModel::AddSrnQuantities,
            " " , py::arg("pOtherSrn"), py::arg("scale") = 1. )
        .def(
            "AddShrunkEdgeSrn",
            (void(AbstractSrnModel::*)(::AbstractSrnModel *)) &AbstractSrnModel::AddShrunkEdgeSrn,
            " " , py::arg("pShrunkEdgeSrn") )
        .def(
            "AddMergedEdgeSrn",
            (void(AbstractSrnModel::*)(::AbstractSrnModel *)) &AbstractSrnModel::AddMergedEdgeSrn,
            " " , py::arg("pMergedEdgeSrn") )
        .def(
            "AddShrunkEdgeToInterior",
            (void(AbstractSrnModel::*)(::AbstractSrnModel *)) &AbstractSrnModel::AddShrunkEdgeToInterior,
            " " , py::arg("pShrunkEdgeSrn") )
        .def(
            "SplitEdgeSrn",
            (void(AbstractSrnModel::*)(double const)) &AbstractSrnModel::SplitEdgeSrn,
            " " , py::arg("relativePosition") )
    ;
}
