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
#include "OnLatticeSimulation.hpp"

#include "OnLatticeSimulation2.cppwg.hpp"

namespace py = pybind11;
typedef OnLatticeSimulation<2 > OnLatticeSimulation2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class OnLatticeSimulation2_Overrides : public OnLatticeSimulation2{
    public:
    using OnLatticeSimulation2::OnLatticeSimulation;
    void OutputAdditionalSimulationSetup(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation2,
            OutputAdditionalSimulationSetup,
                    rParamsFile);
    }
    void OutputSimulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation2,
            OutputSimulationParameters,
                    rParamsFile);
    }
    void UpdateCellPopulation() override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation2,
            UpdateCellPopulation,
            );
    }
    void UpdateCellLocationsAndTopology() override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation2,
            UpdateCellLocationsAndTopology,
            );
    }

};
void register_OnLatticeSimulation2_class(py::module &m){
py::class_<OnLatticeSimulation2 , OnLatticeSimulation2_Overrides , boost::shared_ptr<OnLatticeSimulation2 >  , AbstractCellBasedSimulation<2, 2>  >(m, "OnLatticeSimulation2")
        .def(py::init<::AbstractCellPopulation<2> &, bool, bool >(), py::arg("rCellPopulation"), py::arg("deleteCellPopulationInDestructor") = false, py::arg("initialiseCells") = true)
        .def(
            "AddUpdateRule",
            (void(OnLatticeSimulation2::*)(::boost::shared_ptr<AbstractUpdateRule<2>>)) &OnLatticeSimulation2::AddUpdateRule,
            " " , py::arg("pUpdateRule") )
        .def(
            "RemoveAllUpdateRules",
            (void(OnLatticeSimulation2::*)()) &OnLatticeSimulation2::RemoveAllUpdateRules,
            " "  )
        .def(
            "OutputAdditionalSimulationSetup",
            (void(OnLatticeSimulation2::*)(::out_stream &)) &OnLatticeSimulation2::OutputAdditionalSimulationSetup,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputSimulationParameters",
            (void(OnLatticeSimulation2::*)(::out_stream &)) &OnLatticeSimulation2::OutputSimulationParameters,
            " " , py::arg("rParamsFile") )
    ;
}
