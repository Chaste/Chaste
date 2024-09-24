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
#include "AbstractOdeSystem.hpp"

#include "AbstractOdeSystem.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOdeSystem AbstractOdeSystem;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractOdeSystem_Overrides : public AbstractOdeSystem
{
public:
    using AbstractOdeSystem::AbstractOdeSystem;
    void EvaluateYDerivatives(double time, ::std::vector<double> const & rY, ::std::vector<double> & rDY) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOdeSystem,
            EvaluateYDerivatives,
            time,
            rY,
            rDY);
    }
    bool CalculateStoppingEvent(double time, ::std::vector<double> const & rY) override
    {
        PYBIND11_OVERRIDE(
            bool,
            AbstractOdeSystem,
            CalculateStoppingEvent,
            time,
            rY);
    }
    double CalculateRootFunction(double time, ::std::vector<double> const & rY) override
    {
        PYBIND11_OVERRIDE(
            double,
            AbstractOdeSystem,
            CalculateRootFunction,
            time,
            rY);
    }
};

void register_AbstractOdeSystem_class(py::module &m)
{
    py::class_<AbstractOdeSystem, AbstractOdeSystem_Overrides, boost::shared_ptr<AbstractOdeSystem>>(m, "AbstractOdeSystem")
        .def("EvaluateYDerivatives",
            (void(AbstractOdeSystem::*)(double, ::std::vector<double> const &, ::std::vector<double> &)) &AbstractOdeSystem::EvaluateYDerivatives,
            " ", py::arg("time"), py::arg("rY"), py::arg("rDY"))
        .def("CalculateStoppingEvent",
            (bool(AbstractOdeSystem::*)(double, ::std::vector<double> const &)) &AbstractOdeSystem::CalculateStoppingEvent,
            " ", py::arg("time"), py::arg("rY"))
        .def("CalculateRootFunction",
            (double(AbstractOdeSystem::*)(double, ::std::vector<double> const &)) &AbstractOdeSystem::CalculateRootFunction,
            " ", py::arg("time"), py::arg("rY"))
        .def("GetUseAnalyticJacobian",
            (bool(AbstractOdeSystem::*)()) &AbstractOdeSystem::GetUseAnalyticJacobian,
            " ")
        .def("rGetConstStateVariables",
            (::std::vector<double> const &(AbstractOdeSystem::*)() const) &AbstractOdeSystem::rGetConstStateVariables,
            " ", py::return_value_policy::reference_internal)
    ;
}
