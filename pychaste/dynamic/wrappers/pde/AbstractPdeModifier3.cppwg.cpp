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
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include "PythonPetscObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractPdeModifier.hpp"

#include "AbstractPdeModifier3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractPdeModifier<3 > AbstractPdeModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractPdeModifier3_Overrides : public AbstractPdeModifier3{
    public:
    using AbstractPdeModifier3::AbstractPdeModifier;
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractPdeModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void UpdateAtEndOfOutputTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier3,
            UpdateAtEndOfOutputTimeStep,
                    rCellPopulation);
    }
    void UpdateAtEndOfSolve(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier3,
            UpdateAtEndOfSolve,
                    rCellPopulation);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractPdeModifier3_class(py::module &m){
py::class_<AbstractPdeModifier3 , AbstractPdeModifier3_Overrides , boost::shared_ptr<AbstractPdeModifier3 >   >(m, "AbstractPdeModifier3")
        .def(
            "GetPde",
            (::boost::shared_ptr<AbstractLinearPde<3, 3>>(AbstractPdeModifier3::*)()) &AbstractPdeModifier3::GetPde,
            " "  )
        .def(
            "GetBoundaryCondition",
            (::boost::shared_ptr<AbstractBoundaryCondition<3>>(AbstractPdeModifier3::*)()) &AbstractPdeModifier3::GetBoundaryCondition,
            " "  )
        .def(
            "IsNeumannBoundaryCondition",
            (bool(AbstractPdeModifier3::*)()) &AbstractPdeModifier3::IsNeumannBoundaryCondition,
            " "  )
        .def(
            "SetDependentVariableName",
            (void(AbstractPdeModifier3::*)(::std::string const &)) &AbstractPdeModifier3::SetDependentVariableName,
            " " , py::arg("rName") )
        .def(
            "rGetDependentVariableName",
            (::std::string &(AbstractPdeModifier3::*)()) &AbstractPdeModifier3::rGetDependentVariableName,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "HasAveragedSourcePde",
            (bool(AbstractPdeModifier3::*)()) &AbstractPdeModifier3::HasAveragedSourcePde,
            " "  )
        .def(
            "SetUpSourceTermsForAveragedSourcePde",
            (void(AbstractPdeModifier3::*)(::TetrahedralMesh<3, 3> *, ::std::map<boost::shared_ptr<Cell>, unsigned int> *)) &AbstractPdeModifier3::SetUpSourceTermsForAveragedSourcePde,
            " " , py::arg("pMesh"), py::arg("pCellPdeElementMap") = nullptr )
        .def(
            "SetupSolve",
            (void(AbstractPdeModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &AbstractPdeModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateAtEndOfTimeStep",
            (void(AbstractPdeModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractPdeModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateAtEndOfOutputTimeStep",
            (void(AbstractPdeModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractPdeModifier3::UpdateAtEndOfOutputTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateAtEndOfSolve",
            (void(AbstractPdeModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractPdeModifier3::UpdateAtEndOfSolve,
            " " , py::arg("rCellPopulation") )
        .def(
            "GetOutputGradient",
            (bool(AbstractPdeModifier3::*)()) &AbstractPdeModifier3::GetOutputGradient,
            " "  )
        .def(
            "SetOutputGradient",
            (void(AbstractPdeModifier3::*)(bool)) &AbstractPdeModifier3::SetOutputGradient,
            " " , py::arg("outputGradient") )
        .def(
            "SetOutputSolutionAtPdeNodes",
            (void(AbstractPdeModifier3::*)(bool)) &AbstractPdeModifier3::SetOutputSolutionAtPdeNodes,
            " " , py::arg("outputSolutionAtPdeNodes") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractPdeModifier3::*)(::out_stream &)) &AbstractPdeModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
