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
#include "ImmersedBoundaryLinearMembraneForce.hpp"

#include "ImmersedBoundaryLinearMembraneForce3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryLinearMembraneForce<3 > ImmersedBoundaryLinearMembraneForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryLinearMembraneForce3_Overrides : public ImmersedBoundaryLinearMembraneForce3{
    public:
    using ImmersedBoundaryLinearMembraneForce3::ImmersedBoundaryLinearMembraneForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<3> *, Node<3> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearMembraneForce3,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearMembraneForce3,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryLinearMembraneForce3_class(py::module &m){
py::class_<ImmersedBoundaryLinearMembraneForce3 , ImmersedBoundaryLinearMembraneForce3_Overrides , boost::shared_ptr<ImmersedBoundaryLinearMembraneForce3 >  , AbstractImmersedBoundaryForce<3>  >(m, "ImmersedBoundaryLinearMembraneForce3")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryLinearMembraneForce3::*)(::std::vector<std::pair<Node<3> *, Node<3> *>> &, ::ImmersedBoundaryCellPopulation<3> &)) &ImmersedBoundaryLinearMembraneForce3::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryLinearMembraneForce3::*)(::out_stream &)) &ImmersedBoundaryLinearMembraneForce3::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetElementSpringConst",
            (double(ImmersedBoundaryLinearMembraneForce3::*)() const ) &ImmersedBoundaryLinearMembraneForce3::GetElementSpringConst,
            " "  )
        .def(
            "SetElementSpringConst",
            (void(ImmersedBoundaryLinearMembraneForce3::*)(double)) &ImmersedBoundaryLinearMembraneForce3::SetElementSpringConst,
            " " , py::arg("elementSpringConst") )
        .def(
            "GetElementRestLength",
            (double(ImmersedBoundaryLinearMembraneForce3::*)() const ) &ImmersedBoundaryLinearMembraneForce3::GetElementRestLength,
            " "  )
        .def(
            "SetElementRestLength",
            (void(ImmersedBoundaryLinearMembraneForce3::*)(double)) &ImmersedBoundaryLinearMembraneForce3::SetElementRestLength,
            " " , py::arg("elementRestLength") )
        .def(
            "GetLaminaSpringConst",
            (double(ImmersedBoundaryLinearMembraneForce3::*)() const ) &ImmersedBoundaryLinearMembraneForce3::GetLaminaSpringConst,
            " "  )
        .def(
            "SetLaminaSpringConst",
            (void(ImmersedBoundaryLinearMembraneForce3::*)(double)) &ImmersedBoundaryLinearMembraneForce3::SetLaminaSpringConst,
            " " , py::arg("laminaSpringConst") )
        .def(
            "GetLaminaRestLength",
            (double(ImmersedBoundaryLinearMembraneForce3::*)() const ) &ImmersedBoundaryLinearMembraneForce3::GetLaminaRestLength,
            " "  )
        .def(
            "SetLaminaRestLength",
            (void(ImmersedBoundaryLinearMembraneForce3::*)(double)) &ImmersedBoundaryLinearMembraneForce3::SetLaminaRestLength,
            " " , py::arg("laminaRestLength") )
    ;
}
