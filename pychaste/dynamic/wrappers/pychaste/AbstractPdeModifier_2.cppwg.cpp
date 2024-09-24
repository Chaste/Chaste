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

#include "AbstractPdeModifier_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractPdeModifier<2> AbstractPdeModifier_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractPdeModifier_2_Overrides : public AbstractPdeModifier_2
{
public:
    using AbstractPdeModifier_2::AbstractPdeModifier;
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier_2,
            SetupSolve,
            rCellPopulation,
            outputDirectory);
    }
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractPdeModifier_2,
            UpdateAtEndOfTimeStep,
            rCellPopulation);
    }
    void UpdateAtEndOfOutputTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier_2,
            UpdateAtEndOfOutputTimeStep,
            rCellPopulation);
    }
    void UpdateAtEndOfSolve(::AbstractCellPopulation<2> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier_2,
            UpdateAtEndOfSolve,
            rCellPopulation);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractPdeModifier_2,
            OutputSimulationModifierParameters,
            rParamsFile);
    }
};

void register_AbstractPdeModifier_2_class(py::module &m)
{
    py::class_<AbstractPdeModifier_2, AbstractPdeModifier_2_Overrides, boost::shared_ptr<AbstractPdeModifier_2>, AbstractCellBasedSimulationModifier<2>>(m, "AbstractPdeModifier_2")
        .def("GetPde",
            (::boost::shared_ptr<AbstractLinearPde<2, 2>>(AbstractPdeModifier_2::*)()) &AbstractPdeModifier_2::GetPde,
            " ")
        .def("GetBoundaryCondition",
            (::boost::shared_ptr<AbstractBoundaryCondition<2>>(AbstractPdeModifier_2::*)()) &AbstractPdeModifier_2::GetBoundaryCondition,
            " ")
        .def("IsNeumannBoundaryCondition",
            (bool(AbstractPdeModifier_2::*)()) &AbstractPdeModifier_2::IsNeumannBoundaryCondition,
            " ")
        .def("SetDependentVariableName",
            (void(AbstractPdeModifier_2::*)(::std::string const &)) &AbstractPdeModifier_2::SetDependentVariableName,
            " ", py::arg("rName"))
        .def("rGetDependentVariableName",
            (::std::string &(AbstractPdeModifier_2::*)()) &AbstractPdeModifier_2::rGetDependentVariableName,
            " ", py::return_value_policy::reference_internal)
        .def("HasAveragedSourcePde",
            (bool(AbstractPdeModifier_2::*)()) &AbstractPdeModifier_2::HasAveragedSourcePde,
            " ")
        .def("SetUpSourceTermsForAveragedSourcePde",
            (void(AbstractPdeModifier_2::*)(::TetrahedralMesh<2, 2> *, ::std::map<boost::shared_ptr<Cell>, unsigned int> *)) &AbstractPdeModifier_2::SetUpSourceTermsForAveragedSourcePde,
            " ", py::arg("pMesh"), py::arg("pCellPdeElementMap") = nullptr)
        .def("SetupSolve",
            (void(AbstractPdeModifier_2::*)(::AbstractCellPopulation<2> &, ::std::string)) &AbstractPdeModifier_2::SetupSolve,
            " ", py::arg("rCellPopulation"), py::arg("outputDirectory"))
        .def("UpdateAtEndOfTimeStep",
            (void(AbstractPdeModifier_2::*)(::AbstractCellPopulation<2> &)) &AbstractPdeModifier_2::UpdateAtEndOfTimeStep,
            " ", py::arg("rCellPopulation"))
        .def("UpdateAtEndOfOutputTimeStep",
            (void(AbstractPdeModifier_2::*)(::AbstractCellPopulation<2> &)) &AbstractPdeModifier_2::UpdateAtEndOfOutputTimeStep,
            " ", py::arg("rCellPopulation"))
        .def("UpdateAtEndOfSolve",
            (void(AbstractPdeModifier_2::*)(::AbstractCellPopulation<2> &)) &AbstractPdeModifier_2::UpdateAtEndOfSolve,
            " ", py::arg("rCellPopulation"))
        .def("GetOutputGradient",
            (bool(AbstractPdeModifier_2::*)()) &AbstractPdeModifier_2::GetOutputGradient,
            " ")
        .def("SetOutputGradient",
            (void(AbstractPdeModifier_2::*)(bool)) &AbstractPdeModifier_2::SetOutputGradient,
            " ", py::arg("outputGradient"))
        .def("SetOutputSolutionAtPdeNodes",
            (void(AbstractPdeModifier_2::*)(bool)) &AbstractPdeModifier_2::SetOutputSolutionAtPdeNodes,
            " ", py::arg("outputSolutionAtPdeNodes"))
        .def("OutputSimulationModifierParameters",
            (void(AbstractPdeModifier_2::*)(::out_stream &)) &AbstractPdeModifier_2::OutputSimulationModifierParameters,
            " ", py::arg("rParamsFile"))
    ;
}
