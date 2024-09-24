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
#include "CellPopulationPyChasteActorGenerator.hpp"

#include "CellPopulationPyChasteActorGenerator_2.cppwg.hpp"

namespace py = pybind11;
typedef CellPopulationPyChasteActorGenerator<2> CellPopulationPyChasteActorGenerator_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellPopulationPyChasteActorGenerator_2_Overrides : public CellPopulationPyChasteActorGenerator_2
{
public:
    using CellPopulationPyChasteActorGenerator_2::CellPopulationPyChasteActorGenerator;
    void AddActor(::vtkSmartPointer<vtkRenderer> pRenderer) override
    {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationPyChasteActorGenerator_2,
            AddActor,
            pRenderer);
    }
};

void register_CellPopulationPyChasteActorGenerator_2_class(py::module &m)
{
    py::class_<CellPopulationPyChasteActorGenerator_2, CellPopulationPyChasteActorGenerator_2_Overrides, boost::shared_ptr<CellPopulationPyChasteActorGenerator_2>, AbstractPyChasteActorGenerator<2>>(m, "CellPopulationPyChasteActorGenerator_2")
        .def(py::init<>())
        .def("AddActor",
            (void(CellPopulationPyChasteActorGenerator_2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator_2::AddActor,
            " ", py::arg("pRenderer"))
        .def("AddMeshBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator_2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator_2::AddMeshBasedCellPopulationActor,
            " ", py::arg("pRenderer"))
        .def("AddVertexBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator_2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator_2::AddVertexBasedCellPopulationActor,
            " ", py::arg("pRenderer"))
        .def("AddImmersedBoundaryCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator_2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator_2::AddImmersedBoundaryCellPopulationActor,
            " ", py::arg("pRenderer"))
        .def("AddCaBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator_2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator_2::AddCaBasedCellPopulationActor,
            " ", py::arg("pRenderer"))
        .def("AddPottsBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator_2::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator_2::AddPottsBasedCellPopulationActor,
            " ", py::arg("pRenderer"))
        .def("SetCellPopulation",
            (void(CellPopulationPyChasteActorGenerator_2::*)(::boost::shared_ptr<AbstractCellPopulation<2, 2>>)) &CellPopulationPyChasteActorGenerator_2::SetCellPopulation,
            " ", py::arg("pCellPopulation"))
        .def("SetShowVoronoiMeshEdges",
            (void(CellPopulationPyChasteActorGenerator_2::*)(bool)) &CellPopulationPyChasteActorGenerator_2::SetShowVoronoiMeshEdges,
            " ", py::arg("showEdges"))
        .def("SetShowMutableMeshEdges",
            (void(CellPopulationPyChasteActorGenerator_2::*)(bool)) &CellPopulationPyChasteActorGenerator_2::SetShowMutableMeshEdges,
            " ", py::arg("showEdges"))
        .def("SetShowPottsMeshEdges",
            (void(CellPopulationPyChasteActorGenerator_2::*)(bool)) &CellPopulationPyChasteActorGenerator_2::SetShowPottsMeshEdges,
            " ", py::arg("showEdges"))
        .def("SetShowPottsMeshOutlines",
            (void(CellPopulationPyChasteActorGenerator_2::*)(bool)) &CellPopulationPyChasteActorGenerator_2::SetShowPottsMeshOutlines,
            " ", py::arg("showOutlines"))
        .def("SetColorByCellType",
            (void(CellPopulationPyChasteActorGenerator_2::*)(bool)) &CellPopulationPyChasteActorGenerator_2::SetColorByCellType,
            " ", py::arg("colorByCellType"))
        .def("SetColorByCellMutationState",
            (void(CellPopulationPyChasteActorGenerator_2::*)(bool)) &CellPopulationPyChasteActorGenerator_2::SetColorByCellMutationState,
            " ", py::arg("colorByCellMutationState"))
        .def("SetColorByCellLabel",
            (void(CellPopulationPyChasteActorGenerator_2::*)(bool)) &CellPopulationPyChasteActorGenerator_2::SetColorByCellLabel,
            " ", py::arg("colorByCellLabel"))
        .def("SetColorByUserDefined",
            (void(CellPopulationPyChasteActorGenerator_2::*)(bool)) &CellPopulationPyChasteActorGenerator_2::SetColorByUserDefined,
            " ", py::arg("colorByCellUserDefined"))
        .def("SetColorByCellData",
            (void(CellPopulationPyChasteActorGenerator_2::*)(bool)) &CellPopulationPyChasteActorGenerator_2::SetColorByCellData,
            " ", py::arg("colorByCellData"))
        .def("SetShowCellCentres",
            (void(CellPopulationPyChasteActorGenerator_2::*)(bool)) &CellPopulationPyChasteActorGenerator_2::SetShowCellCentres,
            " ", py::arg("showCentres"))
    ;
}
