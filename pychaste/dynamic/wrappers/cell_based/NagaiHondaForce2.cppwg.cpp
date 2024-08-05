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
#include "NagaiHondaForce.hpp"

#include "NagaiHondaForce2.cppwg.hpp"

namespace py = pybind11;
typedef NagaiHondaForce<2 > NagaiHondaForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class NagaiHondaForce2_Overrides : public NagaiHondaForce2{
    public:
    using NagaiHondaForce2::NagaiHondaForce;
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NagaiHondaForce2,
            AddForceContribution,
                    rCellPopulation);
    }
    double GetAdhesionParameter(::Node<2> * pNodeA, ::Node<2> * pNodeB, ::VertexBasedCellPopulation<2> & rVertexCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            NagaiHondaForce2,
            GetAdhesionParameter,
                    pNodeA,
        pNodeB,
        rVertexCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            NagaiHondaForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_NagaiHondaForce2_class(py::module &m){
py::class_<NagaiHondaForce2 , NagaiHondaForce2_Overrides , boost::shared_ptr<NagaiHondaForce2 >  , AbstractForce<2>  >(m, "NagaiHondaForce2")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(NagaiHondaForce2::*)(::AbstractCellPopulation<2> &)) &NagaiHondaForce2::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "GetAdhesionParameter",
            (double(NagaiHondaForce2::*)(::Node<2> *, ::Node<2> *, ::VertexBasedCellPopulation<2> &)) &NagaiHondaForce2::GetAdhesionParameter,
            " " , py::arg("pNodeA"), py::arg("pNodeB"), py::arg("rVertexCellPopulation") )
        .def(
            "GetNagaiHondaDeformationEnergyParameter",
            (double(NagaiHondaForce2::*)()) &NagaiHondaForce2::GetNagaiHondaDeformationEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaMembraneSurfaceEnergyParameter",
            (double(NagaiHondaForce2::*)()) &NagaiHondaForce2::GetNagaiHondaMembraneSurfaceEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaCellCellAdhesionEnergyParameter",
            (double(NagaiHondaForce2::*)()) &NagaiHondaForce2::GetNagaiHondaCellCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaCellBoundaryAdhesionEnergyParameter",
            (double(NagaiHondaForce2::*)()) &NagaiHondaForce2::GetNagaiHondaCellBoundaryAdhesionEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaTargetAreaParameter",
            (double(NagaiHondaForce2::*)()) &NagaiHondaForce2::GetNagaiHondaTargetAreaParameter,
            " "  )
        .def(
            "SetNagaiHondaDeformationEnergyParameter",
            (void(NagaiHondaForce2::*)(double)) &NagaiHondaForce2::SetNagaiHondaDeformationEnergyParameter,
            " " , py::arg("nagaiHondaDeformationEnergyParameter") )
        .def(
            "SetNagaiHondaMembraneSurfaceEnergyParameter",
            (void(NagaiHondaForce2::*)(double)) &NagaiHondaForce2::SetNagaiHondaMembraneSurfaceEnergyParameter,
            " " , py::arg("nagaiHondaMembraneSurfaceEnergyParameter") )
        .def(
            "SetNagaiHondaCellCellAdhesionEnergyParameter",
            (void(NagaiHondaForce2::*)(double)) &NagaiHondaForce2::SetNagaiHondaCellCellAdhesionEnergyParameter,
            " " , py::arg("nagaiHondaCellCellAdhesionEnergyEnergyParameter") )
        .def(
            "SetNagaiHondaCellBoundaryAdhesionEnergyParameter",
            (void(NagaiHondaForce2::*)(double)) &NagaiHondaForce2::SetNagaiHondaCellBoundaryAdhesionEnergyParameter,
            " " , py::arg("nagaiHondaCellBoundaryAdhesionEnergyParameter") )
        .def(
            "SetNagaiHondaTargetAreaParameter",
            (void(NagaiHondaForce2::*)(double)) &NagaiHondaForce2::SetNagaiHondaTargetAreaParameter,
            " " , py::arg("nagaiHondaTargetAreaParameter") )
        .def(
            "OutputForceParameters",
            (void(NagaiHondaForce2::*)(::out_stream &)) &NagaiHondaForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
