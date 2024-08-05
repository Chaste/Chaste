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
#include "CellPopulationPyChasteActorGenerator.hpp"

#include "CellPopulationPyChasteActorGenerator2.cppwg.hpp"

namespace py = pybind11;
typedef CellPopulationPyChasteActorGenerator<2 > CellPopulationPyChasteActorGenerator2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellPopulationPyChasteActorGenerator2_Overrides : public CellPopulationPyChasteActorGenerator2{
    public:
    using CellPopulationPyChasteActorGenerator2::CellPopulationPyChasteActorGenerator;
    void AddActor(::vtkSmartPointer<vtkRenderer> pRenderer) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationPyChasteActorGenerator2,
            AddActor,
                    pRenderer);
    }

};
void register_CellPopulationPyChasteActorGenerator2_class(py::module &m){
py::class_<CellPopulationPyChasteActorGenerator2 , CellPopulationPyChasteActorGenerator2_Overrides , boost::shared_ptr<CellPopulationPyChasteActorGenerator2 >  , AbstractPyChasteActorGenerator<2>  >(m, "CellPopulationPyChasteActorGenerator2")
        .def(py::init< >())
        .def(
            "AddActor",
            (void(CellPopulationPyChasteActorGenerator2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator2::AddActor,
            " " , py::arg("pRenderer") )
        .def(
            "AddMeshBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator2::AddMeshBasedCellPopulationActor,
            " " , py::arg("pRenderer") )
        .def(
            "AddVertexBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator2::AddVertexBasedCellPopulationActor,
            " " , py::arg("pRenderer") )
        .def(
            "AddImmersedBoundaryCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator2::AddImmersedBoundaryCellPopulationActor,
            " " , py::arg("pRenderer") )
        .def(
            "AddCaBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator2::AddCaBasedCellPopulationActor,
            " " , py::arg("pRenderer") )
        .def(
            "AddPottsBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator2::AddPottsBasedCellPopulationActor,
            " " , py::arg("pRenderer") )
        .def(
            "SetCellPopulation",
            (void(CellPopulationPyChasteActorGenerator2::*)(::boost::shared_ptr<AbstractCellPopulation<2, 2>>)) &CellPopulationPyChasteActorGenerator2::SetCellPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "SetShowVoronoiMeshEdges",
            (void(CellPopulationPyChasteActorGenerator2::*)(bool)) &CellPopulationPyChasteActorGenerator2::SetShowVoronoiMeshEdges,
            " " , py::arg("showEdges") )
        .def(
            "SetShowMutableMeshEdges",
            (void(CellPopulationPyChasteActorGenerator2::*)(bool)) &CellPopulationPyChasteActorGenerator2::SetShowMutableMeshEdges,
            " " , py::arg("showEdges") )
        .def(
            "SetShowPottsMeshEdges",
            (void(CellPopulationPyChasteActorGenerator2::*)(bool)) &CellPopulationPyChasteActorGenerator2::SetShowPottsMeshEdges,
            " " , py::arg("showEdges") )
        .def(
            "SetShowPottsMeshOutlines",
            (void(CellPopulationPyChasteActorGenerator2::*)(bool)) &CellPopulationPyChasteActorGenerator2::SetShowPottsMeshOutlines,
            " " , py::arg("showOutlines") )
        .def(
            "SetColorByCellType",
            (void(CellPopulationPyChasteActorGenerator2::*)(bool)) &CellPopulationPyChasteActorGenerator2::SetColorByCellType,
            " " , py::arg("colorByCellType") )
        .def(
            "SetColorByCellMutationState",
            (void(CellPopulationPyChasteActorGenerator2::*)(bool)) &CellPopulationPyChasteActorGenerator2::SetColorByCellMutationState,
            " " , py::arg("colorByCellMutationState") )
        .def(
            "SetColorByCellLabel",
            (void(CellPopulationPyChasteActorGenerator2::*)(bool)) &CellPopulationPyChasteActorGenerator2::SetColorByCellLabel,
            " " , py::arg("colorByCellLabel") )
        .def(
            "SetColorByUserDefined",
            (void(CellPopulationPyChasteActorGenerator2::*)(bool)) &CellPopulationPyChasteActorGenerator2::SetColorByUserDefined,
            " " , py::arg("colorByCellUserDefined") )
        .def(
            "SetColorByCellData",
            (void(CellPopulationPyChasteActorGenerator2::*)(bool)) &CellPopulationPyChasteActorGenerator2::SetColorByCellData,
            " " , py::arg("colorByCellData") )
        .def(
            "SetShowCellCentres",
            (void(CellPopulationPyChasteActorGenerator2::*)(bool)) &CellPopulationPyChasteActorGenerator2::SetShowCellCentres,
            " " , py::arg("showCentres") )
    ;
}
