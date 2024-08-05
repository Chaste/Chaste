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
#include "ImmersedBoundaryLinearDifferentialAdhesionForce.hpp"

#include "ImmersedBoundaryLinearDifferentialAdhesionForce2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryLinearDifferentialAdhesionForce<2 > ImmersedBoundaryLinearDifferentialAdhesionForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryLinearDifferentialAdhesionForce2_Overrides : public ImmersedBoundaryLinearDifferentialAdhesionForce2{
    public:
    using ImmersedBoundaryLinearDifferentialAdhesionForce2::ImmersedBoundaryLinearDifferentialAdhesionForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<2> *, Node<2> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearDifferentialAdhesionForce2,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearDifferentialAdhesionForce2,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryLinearDifferentialAdhesionForce2_class(py::module &m){
py::class_<ImmersedBoundaryLinearDifferentialAdhesionForce2 , ImmersedBoundaryLinearDifferentialAdhesionForce2_Overrides , boost::shared_ptr<ImmersedBoundaryLinearDifferentialAdhesionForce2 >  , AbstractImmersedBoundaryForce<2>  >(m, "ImmersedBoundaryLinearDifferentialAdhesionForce2")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &, ::ImmersedBoundaryCellPopulation<2> &)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(::out_stream &)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetLabelledCellToLabelledCellSpringConst",
            (double(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)() const ) &ImmersedBoundaryLinearDifferentialAdhesionForce2::GetLabelledCellToLabelledCellSpringConst,
            " "  )
        .def(
            "SetLabelledCellToLabelledCellSpringConst",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(double)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::SetLabelledCellToLabelledCellSpringConst,
            " " , py::arg("labelledCellToLabelledCellSpringConst") )
        .def(
            "GetLabelledCellToCellSpringConst",
            (double(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)() const ) &ImmersedBoundaryLinearDifferentialAdhesionForce2::GetLabelledCellToCellSpringConst,
            " "  )
        .def(
            "SetLabelledCellToCellSpringConst",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(double)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::SetLabelledCellToCellSpringConst,
            " " , py::arg("labelledCellToCellSpringConst") )
        .def(
            "GetCellToCellSpringConst",
            (double(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)() const ) &ImmersedBoundaryLinearDifferentialAdhesionForce2::GetCellToCellSpringConst,
            " "  )
        .def(
            "SetCellToCellSpringConst",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(double)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::SetCellToCellSpringConst,
            " " , py::arg("cellToCellSpringConst") )
        .def(
            "GetRestLength",
            (double(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)() const ) &ImmersedBoundaryLinearDifferentialAdhesionForce2::GetRestLength,
            " "  )
        .def(
            "SetRestLength",
            (void(ImmersedBoundaryLinearDifferentialAdhesionForce2::*)(double)) &ImmersedBoundaryLinearDifferentialAdhesionForce2::SetRestLength,
            " " , py::arg("restLength") )
    ;
}
