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
#include "PythonUblasObjectConverters.hpp"
#include "ImmersedBoundaryMesh.hpp"

#include "ImmersedBoundaryMesh_2_2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryMesh<2, 2> ImmersedBoundaryMesh_2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef unsigned int unsignedint;

class ImmersedBoundaryMesh_2_2_Overrides : public ImmersedBoundaryMesh_2_2
{
public:
    using ImmersedBoundaryMesh_2_2::ImmersedBoundaryMesh;
    unsigned int GetNumNodes() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryMesh_2_2,
            GetNumNodes,
            );
    }
    unsigned int GetNumElements() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryMesh_2_2,
            GetNumElements,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 2> const & rLocation2) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            ImmersedBoundaryMesh_2_2,
            GetVectorFromAtoB,
            rLocation1,
            rLocation2);
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetCentroidOfElement(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            ImmersedBoundaryMesh_2_2,
            GetCentroidOfElement,
            index);
    }
    void Clear() override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMesh_2_2,
            Clear,
            );
    }
    double GetVolumeOfElement(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryMesh_2_2,
            GetVolumeOfElement,
            index);
    }
    double GetSurfaceAreaOfElement(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryMesh_2_2,
            GetSurfaceAreaOfElement,
            index);
    }
    ::boost::numeric::ublas::c_vector<double, 3> CalculateMomentsOfElement(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            ImmersedBoundaryMesh_2_2,
            CalculateMomentsOfElement,
            index);
    }
    unsigned int SolveNodeMapping(unsigned int index) const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryMesh_2_2,
            SolveNodeMapping,
            index);
    }
};

void register_ImmersedBoundaryMesh_2_2_class(py::module &m)
{
    py::class_<ImmersedBoundaryMesh_2_2, ImmersedBoundaryMesh_2_2_Overrides, boost::shared_ptr<ImmersedBoundaryMesh_2_2>, AbstractMesh<2, 2>>(m, "ImmersedBoundaryMesh_2_2")
        .def(py::init<::std::vector<Node<2> *>, ::std::vector<ImmersedBoundaryElement<2, 2> *>, ::std::vector<ImmersedBoundaryElement<1, 2> *>, unsigned int, unsigned int>(), py::arg("nodes"), py::arg("elements"), py::arg("laminas") = ::std::vector<ImmersedBoundaryElement<1, 2> *> {}, py::arg("numGridPtsX") = 128U, py::arg("numGridPtsY") = 128U)
        .def(py::init<>())
        .def("rGetNodes",
            (::std::vector<Node<2> *> const &(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::rGetNodes,
            " ", py::return_value_policy::reference_internal)
        .def("GetElementIteratorBegin",
            (::ImmersedBoundaryMesh<2, 2>::ImmersedBoundaryElementIterator(ImmersedBoundaryMesh_2_2::*)(bool)) &ImmersedBoundaryMesh_2_2::GetElementIteratorBegin,
            " ", py::arg("skipDeletedElements") = true)
        .def("GetElementIteratorEnd",
            (::ImmersedBoundaryMesh<2, 2>::ImmersedBoundaryElementIterator(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::GetElementIteratorEnd,
            " ")
        .def("GetLaminaIteratorBegin",
            (::ImmersedBoundaryMesh<2, 2>::ImmersedBoundaryLaminaIterator(ImmersedBoundaryMesh_2_2::*)(bool)) &ImmersedBoundaryMesh_2_2::GetLaminaIteratorBegin,
            " ", py::arg("skipDeletedLaminas") = true)
        .def("GetLaminaIteratorEnd",
            (::ImmersedBoundaryMesh<2, 2>::ImmersedBoundaryLaminaIterator(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::GetLaminaIteratorEnd,
            " ")
        .def("GetNumNodes",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetNumNodes,
            " ")
        .def("GetNumElements",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetNumElements,
            " ")
        .def("GetNumAllElements",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetNumAllElements,
            " ")
        .def("GetNumLaminas",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetNumLaminas,
            " ")
        .def("GetNumGridPtsX",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetNumGridPtsX,
            " ")
        .def("GetNumGridPtsY",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetNumGridPtsY,
            " ")
        .def("GetCharacteristicNodeSpacing",
            (double(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetCharacteristicNodeSpacing,
            " ")
        .def("GetSpacingRatio",
            (double(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetSpacingRatio,
            " ")
        .def("GetMaxNodeIndex",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetMaxNodeIndex,
            " ")
        .def("GetMaxElementIndex",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetMaxElementIndex,
            " ")
        .def("GetMaxLaminaIndex",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetMaxLaminaIndex,
            " ")
        .def("GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(ImmersedBoundaryMesh_2_2::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &ImmersedBoundaryMesh_2_2::GetVectorFromAtoB,
            " ", py::arg("rLocation1"), py::arg("rLocation2"))
        .def("SetNode",
            (void(ImmersedBoundaryMesh_2_2::*)(unsigned int, ::ChastePoint<2>)) &ImmersedBoundaryMesh_2_2::SetNode,
            " ", py::arg("nodeIndex"), py::arg("point"))
        .def("ConformToGeometry",
            (void(ImmersedBoundaryMesh_2_2::*)(::boost::numeric::ublas::c_vector<double, 2> &)) &ImmersedBoundaryMesh_2_2::ConformToGeometry,
            " ", py::arg("rLocation"))
        .def("rGet2dVelocityGrids",
            (::boost::multi_array<double, 3> const &(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::rGet2dVelocityGrids,
            " ", py::return_value_policy::reference_internal)
        .def("rGetModifiable2dVelocityGrids",
            (::boost::multi_array<double, 3> &(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::rGetModifiable2dVelocityGrids,
            " ", py::return_value_policy::reference_internal)
        .def("SetNumGridPtsX",
            (void(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::SetNumGridPtsX,
            " ", py::arg("meshPointsX"))
        .def("SetNumGridPtsY",
            (void(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::SetNumGridPtsY,
            " ", py::arg("meshPointsY"))
        .def("SetNumGridPtsXAndY",
            (void(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::SetNumGridPtsXAndY,
            " ", py::arg("numGridPts"))
        .def("SetCharacteristicNodeSpacing",
            (void(ImmersedBoundaryMesh_2_2::*)(double)) &ImmersedBoundaryMesh_2_2::SetCharacteristicNodeSpacing,
            " ", py::arg("nodeSpacing"))
        .def("AddNode",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)(::Node<2> *)) &ImmersedBoundaryMesh_2_2::AddNode,
            " ", py::arg("pNewNode"))
        .def("rGetElementFluidSources",
            (::std::vector<std::shared_ptr<FluidSource<2>>> &(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::rGetElementFluidSources,
            " ", py::return_value_policy::reference_internal)
        .def("rGetBalancingFluidSources",
            (::std::vector<std::shared_ptr<FluidSource<2>>> &(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::rGetBalancingFluidSources,
            " ", py::return_value_policy::reference_internal)
        .def("GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::GetNeighbouringNodeIndices,
            " ", py::arg("nodeIndex"))
        .def("GetElement",
            (::ImmersedBoundaryElement<2, 2> *(ImmersedBoundaryMesh_2_2::*)(unsigned int) const) &ImmersedBoundaryMesh_2_2::GetElement,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("GetLamina",
            (::ImmersedBoundaryElement<1, 2> *(ImmersedBoundaryMesh_2_2::*)(unsigned int) const) &ImmersedBoundaryMesh_2_2::GetLamina,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("GetCentroidOfElement",
            (::boost::numeric::ublas::c_vector<double, 2>(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::GetCentroidOfElement,
            " ", py::arg("index"))
        .def("ConstructFromMeshReader",
            (void(ImmersedBoundaryMesh_2_2::*)(::AbstractMeshReader<2, 2> &)) &ImmersedBoundaryMesh_2_2::ConstructFromMeshReader,
            " ", py::arg("rMeshReader"))
        .def("Clear",
            (void(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::Clear,
            " ")
        .def("GetVolumeOfElement",
            (double(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::GetVolumeOfElement,
            " ", py::arg("index"))
        .def("GetSurfaceAreaOfElement",
            (double(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::GetSurfaceAreaOfElement,
            " ", py::arg("index"))
        .def("GetVoronoiSurfaceAreaOfElement",
            (double(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::GetVoronoiSurfaceAreaOfElement,
            " ", py::arg("elemIdx"))
        .def("GetAverageNodeSpacingOfElement",
            (double(ImmersedBoundaryMesh_2_2::*)(unsigned int, bool)) &ImmersedBoundaryMesh_2_2::GetAverageNodeSpacingOfElement,
            " ", py::arg("index"), py::arg("recalculate") = true)
        .def("GetAverageNodeSpacingOfLamina",
            (double(ImmersedBoundaryMesh_2_2::*)(unsigned int, bool)) &ImmersedBoundaryMesh_2_2::GetAverageNodeSpacingOfLamina,
            " ", py::arg("index"), py::arg("recalculate") = true)
        .def("CalculateMomentsOfElement",
            (::boost::numeric::ublas::c_vector<double, 3>(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::CalculateMomentsOfElement,
            " ", py::arg("index"))
        .def("GetElongationShapeFactorOfElement",
            (double(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::GetElongationShapeFactorOfElement,
            " ", py::arg("elementIndex"))
        .def("GetTortuosityOfMesh",
            (double(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::GetTortuosityOfMesh,
            " ")
        .def("GetSkewnessOfElementMassDistributionAboutAxis",
            (double(ImmersedBoundaryMesh_2_2::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 2>)) &ImmersedBoundaryMesh_2_2::GetSkewnessOfElementMassDistributionAboutAxis,
            " ", py::arg("elemIndex"), py::arg("axis"))
        .def("CalculateBoundingBoxOfElement",
            (::ChasteCuboid<2>(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::CalculateBoundingBoxOfElement,
            " ", py::arg("index"))
        .def("GetShortAxisOfElement",
            (::boost::numeric::ublas::c_vector<double, 2>(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::GetShortAxisOfElement,
            " ", py::arg("index"))
        .def("DivideElementAlongGivenAxis",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)(::ImmersedBoundaryElement<2, 2> *, ::boost::numeric::ublas::c_vector<double, 2>, bool)) &ImmersedBoundaryMesh_2_2::DivideElementAlongGivenAxis,
            " ", py::arg("pElement"), py::arg("axisOfDivision"), py::arg("placeOriginalElementBelow") = false)
        .def("DivideElementAlongShortAxis",
            (unsigned int(ImmersedBoundaryMesh_2_2::*)(::ImmersedBoundaryElement<2, 2> *, bool)) &ImmersedBoundaryMesh_2_2::DivideElementAlongShortAxis,
            " ", py::arg("pElement"), py::arg("placeOriginalElementBelow") = false)
        .def("GetElementDivisionSpacing",
            (double(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::GetElementDivisionSpacing,
            " ")
        .def("SetElementDivisionSpacing",
            (void(ImmersedBoundaryMesh_2_2::*)(double)) &ImmersedBoundaryMesh_2_2::SetElementDivisionSpacing,
            " ", py::arg("elementDivisionSpacing"))
        .def("GetNeighbourDist",
            (double(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetNeighbourDist,
            " ")
        .def("SetCellRearrangementThreshold",
            (void(ImmersedBoundaryMesh_2_2::*)(double)) &ImmersedBoundaryMesh_2_2::SetCellRearrangementThreshold,
            " ", py::arg("cellRearrangementThreshold"))
        .def("GetCellRearrangementThreshold",
            (double(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::GetCellRearrangementThreshold,
            " ")
        .def("SetNeighbourDist",
            (void(ImmersedBoundaryMesh_2_2::*)(double)) &ImmersedBoundaryMesh_2_2::SetNeighbourDist,
            " ", py::arg("neighbourDist"))
        .def("UpdateNodeLocationsVoronoiDiagramIfOutOfDate",
            (void(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::UpdateNodeLocationsVoronoiDiagramIfOutOfDate,
            " ")
        .def("ReMesh",
            (void(ImmersedBoundaryMesh_2_2::*)(bool)) &ImmersedBoundaryMesh_2_2::ReMesh,
            " ", py::arg("randomOrder") = false)
        .def("ReMeshElement",
            (void(ImmersedBoundaryMesh_2_2::*)(::ImmersedBoundaryElement<2, 2> *, bool)) &ImmersedBoundaryMesh_2_2::ReMeshElement,
            " ", py::arg("pElement"), py::arg("randomOrder"))
        .def("ReMeshLamina",
            (void(ImmersedBoundaryMesh_2_2::*)(::ImmersedBoundaryElement<1, 2> *, bool)) &ImmersedBoundaryMesh_2_2::ReMeshLamina,
            " ", py::arg("pLamina"), py::arg("randomOrder"))
        .def("NodesInDifferentElementOrLamina",
            (bool(ImmersedBoundaryMesh_2_2::*)(::Node<2> *, ::Node<2> *)) &ImmersedBoundaryMesh_2_2::NodesInDifferentElementOrLamina,
            " ", py::arg("pNodeA"), py::arg("pNodeB"))
        .def("GetNeighbouringElementIndices",
            (::std::set<unsigned int>(ImmersedBoundaryMesh_2_2::*)(unsigned int)) &ImmersedBoundaryMesh_2_2::GetNeighbouringElementIndices,
            " ", py::arg("elemIdx"))
        .def("CalculateLengthOfVoronoiEdge",
            (double(ImmersedBoundaryMesh_2_2::*)(::boost::polygon::voronoi_diagram<double>::edge_type const &)) &ImmersedBoundaryMesh_2_2::CalculateLengthOfVoronoiEdge,
            " ", py::arg("rEdge"))
        .def("GetPolygonDistribution",
            (::std::array<unsigned int, 13>(ImmersedBoundaryMesh_2_2::*)()) &ImmersedBoundaryMesh_2_2::GetPolygonDistribution,
            " ")
        .def("rGetNodeLocationsVoronoiDiagram",
            (::boost::polygon::voronoi_diagram<double> const &(ImmersedBoundaryMesh_2_2::*)(bool)) &ImmersedBoundaryMesh_2_2::rGetNodeLocationsVoronoiDiagram,
            " ", py::arg("update") = true, py::return_value_policy::reference_internal)
        .def("GetVoronoiCellIdsIndexedByNodeIndex",
            (::std::vector<unsigned int> const &(ImmersedBoundaryMesh_2_2::*)() const) &ImmersedBoundaryMesh_2_2::GetVoronoiCellIdsIndexedByNodeIndex,
            " ", py::return_value_policy::reference_internal)
        .def("ScaleUpToVoronoiCoordinate",
            (int(ImmersedBoundaryMesh_2_2::*)(double) const) &ImmersedBoundaryMesh_2_2::ScaleUpToVoronoiCoordinate,
            " ", py::arg("location"))
        .def("ScaleDistanceDownFromVoronoi",
            (double(ImmersedBoundaryMesh_2_2::*)(double const) const) &ImmersedBoundaryMesh_2_2::ScaleDistanceDownFromVoronoi,
            " ", py::arg("distance"))
    ;
}
