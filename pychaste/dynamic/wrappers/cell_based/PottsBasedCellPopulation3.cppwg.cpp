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
#include "AbstractCellBasedSimulation.hpp"
#include "AbstractImmersedBoundaryDivisionRule.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"
#include "BoundaryNodeWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAppliedForceWriter.hpp"
#include "CellCycleModelProteinConcentrationsWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "CellDeltaNotchWriter.hpp"
#include "CellDivisionLocationsWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"
#include "CellPopulationAreaWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellRadiusWriter.hpp"
#include "CellRemovalLocationsWriter.hpp"
#include "CellRosetteRankWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "HeterotypicBoundaryLengthWriter.hpp"
#include "LegacyCellProliferativeTypesWriter.hpp"
#include "NodeLocationWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "PottsMeshWriter.hpp"
#include "RadialCellDataDistributionWriter.hpp"
#include "VertexIntersectionSwapLocationsWriter.hpp"
#include "VertexT1SwapLocationsWriter.hpp"
#include "VertexT2SwapLocationsWriter.hpp"
#include "VertexT3SwapLocationsWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PottsBasedCellPopulation.hpp"

#include "PottsBasedCellPopulation3.cppwg.hpp"

namespace py = pybind11;
typedef PottsBasedCellPopulation<3 > PottsBasedCellPopulation3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::TetrahedralMesh<3, 3> * _TetrahedralMesh_lt_3_3_gt_Ptr;
typedef ::Node<3> * _Node_lt_3_gt_Ptr;
typedef unsigned int unsignedint;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::CellPtr _CellPtr;
typedef unsigned int unsignedint;

class PottsBasedCellPopulation3_Overrides : public PottsBasedCellPopulation3{
    public:
    using PottsBasedCellPopulation3::PottsBasedCellPopulation;
    ::Node<3> * GetNode(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _Node_lt_3_gt_Ptr,
            PottsBasedCellPopulation3,
            GetNode,
                    index);
    }
    unsigned int GetNumNodes() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            PottsBasedCellPopulation3,
            GetNumNodes,
            );
    }
    ::std::set<unsigned int> GetNeighbouringLocationIndices(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            PottsBasedCellPopulation3,
            GetNeighbouringLocationIndices,
                    pCell);
    }
    ::boost::numeric::ublas::c_vector<double, 3> GetLocationOfCellCentre(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            PottsBasedCellPopulation3,
            GetLocationOfCellCentre,
                    pCell);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override {
        PYBIND11_OVERRIDE(
            _CellPtr,
            PottsBasedCellPopulation3,
            AddCell,
                    pNewCell,
        pParentCell);
    }
    unsigned int RemoveDeadCells() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            PottsBasedCellPopulation3,
            RemoveDeadCells,
            );
    }
    void UpdateCellLocations(double dt) override {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation3,
            UpdateCellLocations,
                    dt);
    }
    bool IsCellAssociatedWithADeletedLocation(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            bool,
            PottsBasedCellPopulation3,
            IsCellAssociatedWithADeletedLocation,
                    pCell);
    }
    void Update(bool hasHadBirthsOrDeaths) override {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation3,
            Update,
                    hasHadBirthsOrDeaths);
    }
    void OpenWritersFiles(::OutputFileHandler & rOutputFileHandler) override {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation3,
            OpenWritersFiles,
                    rOutputFileHandler);
    }
    void WriteResultsToFiles(::std::string const & rDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation3,
            WriteResultsToFiles,
                    rDirectory);
    }
    void AcceptPopulationWriter(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>> pPopulationWriter) override {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation3,
            AcceptPopulationWriter,
                    pPopulationWriter);
    }
    void AcceptPopulationCountWriter(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>> pPopulationCountWriter) override {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation3,
            AcceptPopulationCountWriter,
                    pPopulationCountWriter);
    }
    void AcceptPopulationEventWriter(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>> pPopulationEventWriter) override {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation3,
            AcceptPopulationEventWriter,
                    pPopulationEventWriter);
    }
    void AcceptCellWriter(::boost::shared_ptr<AbstractCellWriter<3, 3>> pCellWriter, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation3,
            AcceptCellWriter,
                    pCellWriter,
        pCell);
    }
    double GetVolumeOfCell(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            double,
            PottsBasedCellPopulation3,
            GetVolumeOfCell,
                    pCell);
    }
    double GetWidth(unsigned int const & rDimension) override {
        PYBIND11_OVERRIDE(
            double,
            PottsBasedCellPopulation3,
            GetWidth,
                    rDimension);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation3,
            OutputCellPopulationParameters,
                    rParamsFile);
    }
    void AddUpdateRule(::boost::shared_ptr<AbstractUpdateRule<3>> pUpdateRule) override {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation3,
            AddUpdateRule,
                    pUpdateRule);
    }
    double GetCellDataItemAtPdeNode(unsigned int pdeNodeIndex, ::std::string & rVariableName, bool dirichletBoundaryConditionApplies, double dirichletBoundaryValue) override {
        PYBIND11_OVERRIDE(
            double,
            PottsBasedCellPopulation3,
            GetCellDataItemAtPdeNode,
                    pdeNodeIndex,
        rVariableName,
        dirichletBoundaryConditionApplies,
        dirichletBoundaryValue);
    }

};
void register_PottsBasedCellPopulation3_class(py::module &m){
py::class_<PottsBasedCellPopulation3 , PottsBasedCellPopulation3_Overrides , boost::shared_ptr<PottsBasedCellPopulation3 >  , AbstractOnLatticeCellPopulation<3>  >(m, "PottsBasedCellPopulation3")
        .def(py::init<::PottsMesh<3> &, ::std::vector<boost::shared_ptr<Cell>> &, bool, bool, ::std::vector<unsigned int> const >(), py::arg("rMesh"), py::arg("rCells"), py::arg("deleteMesh") = false, py::arg("validate") = true, py::arg("locationIndices") = std::vector<unsigned int>())
        .def(py::init<::PottsMesh<3> & >(), py::arg("rMesh"))
        .def(
            "GetNumElements",
            (unsigned int(PottsBasedCellPopulation3::*)()) &PottsBasedCellPopulation3::GetNumElements,
            " "  )
        .def(
            "GetNode",
            (::Node<3> *(PottsBasedCellPopulation3::*)(unsigned int)) &PottsBasedCellPopulation3::GetNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNumNodes",
            (unsigned int(PottsBasedCellPopulation3::*)()) &PottsBasedCellPopulation3::GetNumNodes,
            " "  )
        .def(
            "GetNeighbouringLocationIndices",
            (::std::set<unsigned int>(PottsBasedCellPopulation3::*)(::CellPtr)) &PottsBasedCellPopulation3::GetNeighbouringLocationIndices,
            " " , py::arg("pCell") )
        .def(
            "GetLocationOfCellCentre",
            (::boost::numeric::ublas::c_vector<double, 3>(PottsBasedCellPopulation3::*)(::CellPtr)) &PottsBasedCellPopulation3::GetLocationOfCellCentre,
            " " , py::arg("pCell") )
        .def(
            "AddCell",
            (::CellPtr(PottsBasedCellPopulation3::*)(::CellPtr, ::CellPtr)) &PottsBasedCellPopulation3::AddCell,
            " " , py::arg("pNewCell"), py::arg("pParentCell") = ::CellPtr( ) )
        .def(
            "RemoveDeadCells",
            (unsigned int(PottsBasedCellPopulation3::*)()) &PottsBasedCellPopulation3::RemoveDeadCells,
            " "  )
        .def(
            "UpdateCellLocations",
            (void(PottsBasedCellPopulation3::*)(double)) &PottsBasedCellPopulation3::UpdateCellLocations,
            " " , py::arg("dt") )
        .def(
            "IsCellAssociatedWithADeletedLocation",
            (bool(PottsBasedCellPopulation3::*)(::CellPtr)) &PottsBasedCellPopulation3::IsCellAssociatedWithADeletedLocation,
            " " , py::arg("pCell") )
        .def(
            "Update",
            (void(PottsBasedCellPopulation3::*)(bool)) &PottsBasedCellPopulation3::Update,
            " " , py::arg("hasHadBirthsOrDeaths") = true )
        .def(
            "OpenWritersFiles",
            (void(PottsBasedCellPopulation3::*)(::OutputFileHandler &)) &PottsBasedCellPopulation3::OpenWritersFiles,
            " " , py::arg("rOutputFileHandler") )
        .def(
            "WriteResultsToFiles",
            (void(PottsBasedCellPopulation3::*)(::std::string const &)) &PottsBasedCellPopulation3::WriteResultsToFiles,
            " " , py::arg("rDirectory") )
        .def(
            "AcceptPopulationWriter",
            (void(PottsBasedCellPopulation3::*)(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>>)) &PottsBasedCellPopulation3::AcceptPopulationWriter,
            " " , py::arg("pPopulationWriter") )
        .def(
            "AcceptPopulationCountWriter",
            (void(PottsBasedCellPopulation3::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>>)) &PottsBasedCellPopulation3::AcceptPopulationCountWriter,
            " " , py::arg("pPopulationCountWriter") )
        .def(
            "AcceptPopulationEventWriter",
            (void(PottsBasedCellPopulation3::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>>)) &PottsBasedCellPopulation3::AcceptPopulationEventWriter,
            " " , py::arg("pPopulationEventWriter") )
        .def(
            "AcceptCellWriter",
            (void(PottsBasedCellPopulation3::*)(::boost::shared_ptr<AbstractCellWriter<3, 3>>, ::CellPtr)) &PottsBasedCellPopulation3::AcceptCellWriter,
            " " , py::arg("pCellWriter"), py::arg("pCell") )
        .def(
            "GetVolumeOfCell",
            (double(PottsBasedCellPopulation3::*)(::CellPtr)) &PottsBasedCellPopulation3::GetVolumeOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetWidth",
            (double(PottsBasedCellPopulation3::*)(unsigned int const &)) &PottsBasedCellPopulation3::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "OutputCellPopulationParameters",
            (void(PottsBasedCellPopulation3::*)(::out_stream &)) &PottsBasedCellPopulation3::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "SetTemperature",
            (void(PottsBasedCellPopulation3::*)(double)) &PottsBasedCellPopulation3::SetTemperature,
            " " , py::arg("temperature") )
        .def(
            "GetTemperature",
            (double(PottsBasedCellPopulation3::*)()) &PottsBasedCellPopulation3::GetTemperature,
            " "  )
        .def(
            "SetNumSweepsPerTimestep",
            (void(PottsBasedCellPopulation3::*)(unsigned int)) &PottsBasedCellPopulation3::SetNumSweepsPerTimestep,
            " " , py::arg("numSweepsPerTimestep") )
        .def(
            "GetNumSweepsPerTimestep",
            (unsigned int(PottsBasedCellPopulation3::*)()) &PottsBasedCellPopulation3::GetNumSweepsPerTimestep,
            " "  )
        .def(
            "CreateElementTessellation",
            (void(PottsBasedCellPopulation3::*)()) &PottsBasedCellPopulation3::CreateElementTessellation,
            " "  )
        .def(
            "CreateMutableMesh",
            (void(PottsBasedCellPopulation3::*)()) &PottsBasedCellPopulation3::CreateMutableMesh,
            " "  )
        .def(
            "AddUpdateRule",
            (void(PottsBasedCellPopulation3::*)(::boost::shared_ptr<AbstractUpdateRule<3>>)) &PottsBasedCellPopulation3::AddUpdateRule,
            " " , py::arg("pUpdateRule") )
        .def(
            "GetCellDataItemAtPdeNode",
            (double(PottsBasedCellPopulation3::*)(unsigned int, ::std::string &, bool, double)) &PottsBasedCellPopulation3::GetCellDataItemAtPdeNode,
            " " , py::arg("pdeNodeIndex"), py::arg("rVariableName"), py::arg("dirichletBoundaryConditionApplies") = false, py::arg("dirichletBoundaryValue") = 0. )
        .def("AddPopulationWriterVoronoiDataWriter", &PottsBasedCellPopulation3::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &PottsBasedCellPopulation3::AddCellWriter<CellLabelWriter>)
    ;
}
