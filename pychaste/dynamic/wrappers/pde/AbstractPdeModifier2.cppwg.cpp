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

#include "AbstractPdeModifier2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractPdeModifier<2 > AbstractPdeModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractPdeModifier2_Overrides : public AbstractPdeModifier2{
    public:
    using AbstractPdeModifier2::AbstractPdeModifier;
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractPdeModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void UpdateAtEndOfOutputTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier2,
            UpdateAtEndOfOutputTimeStep,
                    rCellPopulation);
    }
    void UpdateAtEndOfSolve(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier2,
            UpdateAtEndOfSolve,
                    rCellPopulation);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractPdeModifier2_class(py::module &m){
py::class_<AbstractPdeModifier2 , AbstractPdeModifier2_Overrides , boost::shared_ptr<AbstractPdeModifier2 >   >(m, "AbstractPdeModifier2")
        .def(
            "GetPde",
            (::boost::shared_ptr<AbstractLinearPde<2, 2>>(AbstractPdeModifier2::*)()) &AbstractPdeModifier2::GetPde,
            " "  )
        .def(
            "GetBoundaryCondition",
            (::boost::shared_ptr<AbstractBoundaryCondition<2>>(AbstractPdeModifier2::*)()) &AbstractPdeModifier2::GetBoundaryCondition,
            " "  )
        .def(
            "IsNeumannBoundaryCondition",
            (bool(AbstractPdeModifier2::*)()) &AbstractPdeModifier2::IsNeumannBoundaryCondition,
            " "  )
        .def(
            "SetDependentVariableName",
            (void(AbstractPdeModifier2::*)(::std::string const &)) &AbstractPdeModifier2::SetDependentVariableName,
            " " , py::arg("rName") )
        .def(
            "rGetDependentVariableName",
            (::std::string &(AbstractPdeModifier2::*)()) &AbstractPdeModifier2::rGetDependentVariableName,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "HasAveragedSourcePde",
            (bool(AbstractPdeModifier2::*)()) &AbstractPdeModifier2::HasAveragedSourcePde,
            " "  )
        .def(
            "SetUpSourceTermsForAveragedSourcePde",
            (void(AbstractPdeModifier2::*)(::TetrahedralMesh<2, 2> *, ::std::map<boost::shared_ptr<Cell>, unsigned int> *)) &AbstractPdeModifier2::SetUpSourceTermsForAveragedSourcePde,
            " " , py::arg("pMesh"), py::arg("pCellPdeElementMap") = nullptr )
        .def(
            "SetupSolve",
            (void(AbstractPdeModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &AbstractPdeModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateAtEndOfTimeStep",
            (void(AbstractPdeModifier2::*)(::AbstractCellPopulation<2> &)) &AbstractPdeModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateAtEndOfOutputTimeStep",
            (void(AbstractPdeModifier2::*)(::AbstractCellPopulation<2> &)) &AbstractPdeModifier2::UpdateAtEndOfOutputTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateAtEndOfSolve",
            (void(AbstractPdeModifier2::*)(::AbstractCellPopulation<2> &)) &AbstractPdeModifier2::UpdateAtEndOfSolve,
            " " , py::arg("rCellPopulation") )
        .def(
            "GetOutputGradient",
            (bool(AbstractPdeModifier2::*)()) &AbstractPdeModifier2::GetOutputGradient,
            " "  )
        .def(
            "SetOutputGradient",
            (void(AbstractPdeModifier2::*)(bool)) &AbstractPdeModifier2::SetOutputGradient,
            " " , py::arg("outputGradient") )
        .def(
            "SetOutputSolutionAtPdeNodes",
            (void(AbstractPdeModifier2::*)(bool)) &AbstractPdeModifier2::SetOutputSolutionAtPdeNodes,
            " " , py::arg("outputSolutionAtPdeNodes") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractPdeModifier2::*)(::out_stream &)) &AbstractPdeModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
