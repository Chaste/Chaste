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
#include "CaBasedCellPopulation.hpp"

#include "CaBasedCellPopulation3.cppwg.hpp"

namespace py = pybind11;
typedef CaBasedCellPopulation<3 > CaBasedCellPopulation3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::TetrahedralMesh<3, 3> * _TetrahedralMesh_lt_3_3_gt_Ptr;
typedef ::Node<3> * _Node_lt_3_gt_Ptr;
typedef unsigned int unsignedint;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::CellPtr _CellPtr;
typedef unsigned int unsignedint;
typedef ::std::vector<boost::shared_ptr<AbstractUpdateRule<3>>> const _std_vector_lt_boost_shared_ptr_lt_AbstractUpdateRule_lt_3_gt__gt__gt_const;

class CaBasedCellPopulation3_Overrides : public CaBasedCellPopulation3{
    public:
    using CaBasedCellPopulation3::CaBasedCellPopulation;
    bool IsSiteAvailable(unsigned int index, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            bool,
            CaBasedCellPopulation3,
            IsSiteAvailable,
                    index,
        pCell);
    }
    ::Node<3> * GetNode(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _Node_lt_3_gt_Ptr,
            CaBasedCellPopulation3,
            GetNode,
                    index);
    }
    unsigned int GetNumNodes() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            CaBasedCellPopulation3,
            GetNumNodes,
            );
    }
    ::std::set<unsigned int> GetNeighbouringLocationIndices(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            CaBasedCellPopulation3,
            GetNeighbouringLocationIndices,
                    pCell);
    }
    ::boost::numeric::ublas::c_vector<double, 3> GetLocationOfCellCentre(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            CaBasedCellPopulation3,
            GetLocationOfCellCentre,
                    pCell);
    }
    void AddCellUsingLocationIndex(unsigned int index, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            AddCellUsingLocationIndex,
                    index,
        pCell);
    }
    void RemoveCellUsingLocationIndex(unsigned int index, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            RemoveCellUsingLocationIndex,
                    index,
        pCell);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override {
        PYBIND11_OVERRIDE(
            _CellPtr,
            CaBasedCellPopulation3,
            AddCell,
                    pNewCell,
        pParentCell);
    }
    double EvaluateDivisionPropensity(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            double,
            CaBasedCellPopulation3,
            EvaluateDivisionPropensity,
                    currentNodeIndex,
        targetNodeIndex,
        pCell);
    }
    unsigned int RemoveDeadCells() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            CaBasedCellPopulation3,
            RemoveDeadCells,
            );
    }
    void OpenWritersFiles(::OutputFileHandler & rOutputFileHandler) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            OpenWritersFiles,
                    rOutputFileHandler);
    }
    void UpdateCellLocations(double dt) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            UpdateCellLocations,
                    dt);
    }
    bool IsCellAssociatedWithADeletedLocation(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            bool,
            CaBasedCellPopulation3,
            IsCellAssociatedWithADeletedLocation,
                    pCell);
    }
    void Update(bool hasHadBirthsOrDeaths) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            Update,
                    hasHadBirthsOrDeaths);
    }
    void AcceptPopulationWriter(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>> pPopulationWriter) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            AcceptPopulationWriter,
                    pPopulationWriter);
    }
    void AcceptPopulationCountWriter(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>> pPopulationCountWriter) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            AcceptPopulationCountWriter,
                    pPopulationCountWriter);
    }
    void AcceptPopulationEventWriter(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>> pPopulationEventWriter) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            AcceptPopulationEventWriter,
                    pPopulationEventWriter);
    }
    void AcceptCellWriter(::boost::shared_ptr<AbstractCellWriter<3, 3>> pCellWriter, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            AcceptCellWriter,
                    pCellWriter,
        pCell);
    }
    double GetVolumeOfCell(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            double,
            CaBasedCellPopulation3,
            GetVolumeOfCell,
                    pCell);
    }
    double GetWidth(unsigned int const & rDimension) override {
        PYBIND11_OVERRIDE(
            double,
            CaBasedCellPopulation3,
            GetWidth,
                    rDimension);
    }
    void RemoveAllUpdateRules() override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            RemoveAllUpdateRules,
            );
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            OutputCellPopulationParameters,
                    rParamsFile);
    }
    bool IsRoomToDivide(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            bool,
            CaBasedCellPopulation3,
            IsRoomToDivide,
                    pCell);
    }
    void AddUpdateRule(::boost::shared_ptr<AbstractUpdateRule<3>> pUpdateRule) override {
        PYBIND11_OVERRIDE(
            void,
            CaBasedCellPopulation3,
            AddUpdateRule,
                    pUpdateRule);
    }
    ::std::vector<boost::shared_ptr<AbstractUpdateRule<3>>> const GetUpdateRuleCollection() const  override {
        PYBIND11_OVERRIDE(
            _std_vector_lt_boost_shared_ptr_lt_AbstractUpdateRule_lt_3_gt__gt__gt_const,
            CaBasedCellPopulation3,
            GetUpdateRuleCollection,
            );
    }
    double GetCellDataItemAtPdeNode(unsigned int pdeNodeIndex, ::std::string & rVariableName, bool dirichletBoundaryConditionApplies, double dirichletBoundaryValue) override {
        PYBIND11_OVERRIDE(
            double,
            CaBasedCellPopulation3,
            GetCellDataItemAtPdeNode,
                    pdeNodeIndex,
        rVariableName,
        dirichletBoundaryConditionApplies,
        dirichletBoundaryValue);
    }
    bool IsPdeNodeAssociatedWithNonApoptoticCell(unsigned int pdeNodeIndex) override {
        PYBIND11_OVERRIDE(
            bool,
            CaBasedCellPopulation3,
            IsPdeNodeAssociatedWithNonApoptoticCell,
                    pdeNodeIndex);
    }

};
void register_CaBasedCellPopulation3_class(py::module &m){
py::class_<CaBasedCellPopulation3 , CaBasedCellPopulation3_Overrides , boost::shared_ptr<CaBasedCellPopulation3 >  , AbstractOnLatticeCellPopulation<3>  >(m, "CaBasedCellPopulation3")
        .def(py::init<::PottsMesh<3> &, ::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const, unsigned int, bool, bool >(), py::arg("rMesh"), py::arg("rCells"), py::arg("locationIndices"), py::arg("latticeCarryingCapacity") = 1U, py::arg("deleteMesh") = false, py::arg("validate") = false)
        .def(py::init<::PottsMesh<3> & >(), py::arg("rMesh"))
        .def(
            "IsSiteAvailable",
            (bool(CaBasedCellPopulation3::*)(unsigned int, ::CellPtr)) &CaBasedCellPopulation3::IsSiteAvailable,
            " " , py::arg("index"), py::arg("pCell") )
        .def(
            "GetNode",
            (::Node<3> *(CaBasedCellPopulation3::*)(unsigned int)) &CaBasedCellPopulation3::GetNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNumNodes",
            (unsigned int(CaBasedCellPopulation3::*)()) &CaBasedCellPopulation3::GetNumNodes,
            " "  )
        .def(
            "GetNeighbouringLocationIndices",
            (::std::set<unsigned int>(CaBasedCellPopulation3::*)(::CellPtr)) &CaBasedCellPopulation3::GetNeighbouringLocationIndices,
            " " , py::arg("pCell") )
        .def(
            "GetLocationOfCellCentre",
            (::boost::numeric::ublas::c_vector<double, 3>(CaBasedCellPopulation3::*)(::CellPtr)) &CaBasedCellPopulation3::GetLocationOfCellCentre,
            " " , py::arg("pCell") )
        .def(
            "AddCellUsingLocationIndex",
            (void(CaBasedCellPopulation3::*)(unsigned int, ::CellPtr)) &CaBasedCellPopulation3::AddCellUsingLocationIndex,
            " " , py::arg("index"), py::arg("pCell") )
        .def(
            "RemoveCellUsingLocationIndex",
            (void(CaBasedCellPopulation3::*)(unsigned int, ::CellPtr)) &CaBasedCellPopulation3::RemoveCellUsingLocationIndex,
            " " , py::arg("index"), py::arg("pCell") )
        .def(
            "AddCell",
            (::CellPtr(CaBasedCellPopulation3::*)(::CellPtr, ::CellPtr)) &CaBasedCellPopulation3::AddCell,
            " " , py::arg("pNewCell"), py::arg("pParentCell") = ::CellPtr( ) )
        .def(
            "EvaluateDivisionPropensity",
            (double(CaBasedCellPopulation3::*)(unsigned int, unsigned int, ::CellPtr)) &CaBasedCellPopulation3::EvaluateDivisionPropensity,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("pCell") )
        .def(
            "RemoveDeadCells",
            (unsigned int(CaBasedCellPopulation3::*)()) &CaBasedCellPopulation3::RemoveDeadCells,
            " "  )
        .def(
            "OpenWritersFiles",
            (void(CaBasedCellPopulation3::*)(::OutputFileHandler &)) &CaBasedCellPopulation3::OpenWritersFiles,
            " " , py::arg("rOutputFileHandler") )
        .def(
            "UpdateCellLocations",
            (void(CaBasedCellPopulation3::*)(double)) &CaBasedCellPopulation3::UpdateCellLocations,
            " " , py::arg("dt") )
        .def(
            "IsCellAssociatedWithADeletedLocation",
            (bool(CaBasedCellPopulation3::*)(::CellPtr)) &CaBasedCellPopulation3::IsCellAssociatedWithADeletedLocation,
            " " , py::arg("pCell") )
        .def(
            "Update",
            (void(CaBasedCellPopulation3::*)(bool)) &CaBasedCellPopulation3::Update,
            " " , py::arg("hasHadBirthsOrDeaths") = true )
        .def(
            "AcceptPopulationWriter",
            (void(CaBasedCellPopulation3::*)(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>>)) &CaBasedCellPopulation3::AcceptPopulationWriter,
            " " , py::arg("pPopulationWriter") )
        .def(
            "AcceptPopulationCountWriter",
            (void(CaBasedCellPopulation3::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>>)) &CaBasedCellPopulation3::AcceptPopulationCountWriter,
            " " , py::arg("pPopulationCountWriter") )
        .def(
            "AcceptPopulationEventWriter",
            (void(CaBasedCellPopulation3::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>>)) &CaBasedCellPopulation3::AcceptPopulationEventWriter,
            " " , py::arg("pPopulationEventWriter") )
        .def(
            "AcceptCellWriter",
            (void(CaBasedCellPopulation3::*)(::boost::shared_ptr<AbstractCellWriter<3, 3>>, ::CellPtr)) &CaBasedCellPopulation3::AcceptCellWriter,
            " " , py::arg("pCellWriter"), py::arg("pCell") )
        .def(
            "GetVolumeOfCell",
            (double(CaBasedCellPopulation3::*)(::CellPtr)) &CaBasedCellPopulation3::GetVolumeOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetWidth",
            (double(CaBasedCellPopulation3::*)(unsigned int const &)) &CaBasedCellPopulation3::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "RemoveAllUpdateRules",
            (void(CaBasedCellPopulation3::*)()) &CaBasedCellPopulation3::RemoveAllUpdateRules,
            " "  )
        .def(
            "OutputCellPopulationParameters",
            (void(CaBasedCellPopulation3::*)(::out_stream &)) &CaBasedCellPopulation3::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "IsRoomToDivide",
            (bool(CaBasedCellPopulation3::*)(::CellPtr)) &CaBasedCellPopulation3::IsRoomToDivide,
            " " , py::arg("pCell") )
        .def(
            "GetCaBasedDivisionRule",
            (::boost::shared_ptr<AbstractCaBasedDivisionRule<3>>(CaBasedCellPopulation3::*)()) &CaBasedCellPopulation3::GetCaBasedDivisionRule,
            " "  )
        .def(
            "SetCaBasedDivisionRule",
            (void(CaBasedCellPopulation3::*)(::boost::shared_ptr<AbstractCaBasedDivisionRule<3>>)) &CaBasedCellPopulation3::SetCaBasedDivisionRule,
            " " , py::arg("pCaBasedDivisionRule") )
        .def(
            "AddUpdateRule",
            (void(CaBasedCellPopulation3::*)(::boost::shared_ptr<AbstractUpdateRule<3>>)) &CaBasedCellPopulation3::AddUpdateRule,
            " " , py::arg("pUpdateRule") )
        .def(
            "GetUpdateRuleCollection",
            (::std::vector<boost::shared_ptr<AbstractUpdateRule<3>>> const(CaBasedCellPopulation3::*)() const ) &CaBasedCellPopulation3::GetUpdateRuleCollection,
            " "  )
        .def(
            "GetCellDataItemAtPdeNode",
            (double(CaBasedCellPopulation3::*)(unsigned int, ::std::string &, bool, double)) &CaBasedCellPopulation3::GetCellDataItemAtPdeNode,
            " " , py::arg("pdeNodeIndex"), py::arg("rVariableName"), py::arg("dirichletBoundaryConditionApplies") = false, py::arg("dirichletBoundaryValue") = 0. )
        .def(
            "IsPdeNodeAssociatedWithNonApoptoticCell",
            (bool(CaBasedCellPopulation3::*)(unsigned int)) &CaBasedCellPopulation3::IsPdeNodeAssociatedWithNonApoptoticCell,
            " " , py::arg("pdeNodeIndex") )
        .def("AddPopulationWriterVoronoiDataWriter", &CaBasedCellPopulation3::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &CaBasedCellPopulation3::AddCellWriter<CellLabelWriter>)
    ;
}
