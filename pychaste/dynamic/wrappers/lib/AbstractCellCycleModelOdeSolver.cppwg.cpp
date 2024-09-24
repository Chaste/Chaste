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
#include "AbstractCellCycleModelOdeSolver.hpp"

#include "AbstractCellCycleModelOdeSolver.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellCycleModelOdeSolver AbstractCellCycleModelOdeSolver;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCellCycleModelOdeSolver_Overrides : public AbstractCellCycleModelOdeSolver
{
public:
    using AbstractCellCycleModelOdeSolver::AbstractCellCycleModelOdeSolver;
    bool IsSetUp() override
    {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractCellCycleModelOdeSolver,
            IsSetUp,
            );
    }
    void Reset() override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellCycleModelOdeSolver,
            Reset,
            );
    }
    void Initialise() override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellCycleModelOdeSolver,
            Initialise,
            );
    }
    bool IsAdaptive() override
    {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCellCycleModelOdeSolver,
            IsAdaptive,
            );
    }
};

void register_AbstractCellCycleModelOdeSolver_class(py::module &m)
{
    py::class_<AbstractCellCycleModelOdeSolver, AbstractCellCycleModelOdeSolver_Overrides, boost::shared_ptr<AbstractCellCycleModelOdeSolver>>(m, "AbstractCellCycleModelOdeSolver")
        .def(py::init<>())
        .def("IsSetUp",
            (bool(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::IsSetUp,
            " ")
        .def("Reset",
            (void(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::Reset,
            " ")
        .def("SolveAndUpdateStateVariable",
            (void(AbstractCellCycleModelOdeSolver::*)(::AbstractOdeSystem *, double, double, double)) &AbstractCellCycleModelOdeSolver::SolveAndUpdateStateVariable,
            " ", py::arg("pAbstractOdeSystem"), py::arg("startTime"), py::arg("endTime"), py::arg("timeStep"))
        .def("Initialise",
            (void(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::Initialise,
            " ")
        .def("StoppingEventOccurred",
            (bool(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::StoppingEventOccurred,
            " ")
        .def("GetStoppingTime",
            (double(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::GetStoppingTime,
            " ")
        .def("SetSizeOfOdeSystem",
            (void(AbstractCellCycleModelOdeSolver::*)(unsigned int)) &AbstractCellCycleModelOdeSolver::SetSizeOfOdeSystem,
            " ", py::arg("sizeOfOdeSystem"))
        .def("GetSizeOfOdeSystem",
            (unsigned int(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::GetSizeOfOdeSystem,
            " ")
        .def("CheckForStoppingEvents",
            (void(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::CheckForStoppingEvents,
            " ")
        .def("SetMaxSteps",
            (void(AbstractCellCycleModelOdeSolver::*)(long int)) &AbstractCellCycleModelOdeSolver::SetMaxSteps,
            " ", py::arg("numSteps"))
        .def("SetTolerances",
            (void(AbstractCellCycleModelOdeSolver::*)(double, double)) &AbstractCellCycleModelOdeSolver::SetTolerances,
            " ", py::arg("relTol") = 1.0E-4, py::arg("absTol") = 9.9999999999999995E-7)
        .def("IsAdaptive",
            (bool(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::IsAdaptive,
            " ")
    ;
}
