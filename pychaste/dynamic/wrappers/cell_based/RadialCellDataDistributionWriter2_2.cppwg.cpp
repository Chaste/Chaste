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
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RadialCellDataDistributionWriter.hpp"

#include "RadialCellDataDistributionWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef RadialCellDataDistributionWriter<2,2 > RadialCellDataDistributionWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class RadialCellDataDistributionWriter2_2_Overrides : public RadialCellDataDistributionWriter2_2{
    public:
    using RadialCellDataDistributionWriter2_2::RadialCellDataDistributionWriter;
    void Visit(::MeshBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }

};
void register_RadialCellDataDistributionWriter2_2_class(py::module &m){
py::class_<RadialCellDataDistributionWriter2_2 , RadialCellDataDistributionWriter2_2_Overrides , boost::shared_ptr<RadialCellDataDistributionWriter2_2 >  , AbstractCellPopulationWriter<2, 2>  >(m, "RadialCellDataDistributionWriter2_2")
        .def(py::init< >())
        .def(
            "VisitAnyPopulation",
            (void(RadialCellDataDistributionWriter2_2::*)(::AbstractCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::VisitAnyPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::MeshBasedCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::CaBasedCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::NodeBasedCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::PottsBasedCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::VertexBasedCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::ImmersedBoundaryCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "SetVariableName",
            (void(RadialCellDataDistributionWriter2_2::*)(::std::string)) &RadialCellDataDistributionWriter2_2::SetVariableName,
            " " , py::arg("variableName") )
        .def(
            "GetVariableName",
            (::std::string(RadialCellDataDistributionWriter2_2::*)() const ) &RadialCellDataDistributionWriter2_2::GetVariableName,
            " "  )
        .def(
            "SetNumRadialBins",
            (void(RadialCellDataDistributionWriter2_2::*)(unsigned int)) &RadialCellDataDistributionWriter2_2::SetNumRadialBins,
            " " , py::arg("numRadialBins") )
        .def(
            "GetNumRadialBins",
            (unsigned int(RadialCellDataDistributionWriter2_2::*)() const ) &RadialCellDataDistributionWriter2_2::GetNumRadialBins,
            " "  )
    ;
}
