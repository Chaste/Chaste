#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryMesh.hpp"

#include "ImmersedBoundaryMesh2_2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryMesh<2,2 > ImmersedBoundaryMesh2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef unsigned int unsignedint;

class ImmersedBoundaryMesh2_2_Overrides : public ImmersedBoundaryMesh2_2{
    public:
    using ImmersedBoundaryMesh2_2::ImmersedBoundaryMesh;
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryMesh2_2,
            GetNumNodes,
            );
    }
    unsigned int GetNumElements() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryMesh2_2,
            GetNumElements,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 2> const & rLocation2) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            ImmersedBoundaryMesh2_2,
            GetVectorFromAtoB,
                    rLocation1,
        rLocation2);
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetCentroidOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            ImmersedBoundaryMesh2_2,
            GetCentroidOfElement,
                    index);
    }
    void Clear() override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMesh2_2,
            Clear,
            );
    }
    double GetVolumeOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryMesh2_2,
            GetVolumeOfElement,
                    index);
    }
    double GetSurfaceAreaOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryMesh2_2,
            GetSurfaceAreaOfElement,
                    index);
    }
    ::boost::numeric::ublas::c_vector<double, 3> CalculateMomentsOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            ImmersedBoundaryMesh2_2,
            CalculateMomentsOfElement,
                    index);
    }
    unsigned int SolveNodeMapping(unsigned int index) const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryMesh2_2,
            SolveNodeMapping,
                    index);
    }

};
void register_ImmersedBoundaryMesh2_2_class(py::module &m){
py::class_<ImmersedBoundaryMesh2_2 , ImmersedBoundaryMesh2_2_Overrides , boost::shared_ptr<ImmersedBoundaryMesh2_2 >  , AbstractMesh<2, 2>  >(m, "ImmersedBoundaryMesh2_2")
        .def(py::init<::std::vector<Node<2> *>, ::std::vector<ImmersedBoundaryElement<2, 2> *>, ::std::vector<ImmersedBoundaryElement<1, 2> *>, unsigned int, unsigned int >(), py::arg("nodes"), py::arg("elements"), py::arg("laminas") = ::std::vector<ImmersedBoundaryElement<1, 2> *> {}, py::arg("numGridPtsX") = 128U, py::arg("numGridPtsY") = 128U)
        .def(py::init< >())
        .def(
            "rGetNodes",
            (::std::vector<Node<2> *> const &(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::rGetNodes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetElementIteratorBegin",
            (::ImmersedBoundaryMesh<2, 2>::ImmersedBoundaryElementIterator(ImmersedBoundaryMesh2_2::*)(bool)) &ImmersedBoundaryMesh2_2::GetElementIteratorBegin,
            " " , py::arg("skipDeletedElements") = true )
        .def(
            "GetElementIteratorEnd",
            (::ImmersedBoundaryMesh<2, 2>::ImmersedBoundaryElementIterator(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::GetElementIteratorEnd,
            " "  )
        .def(
            "GetLaminaIteratorBegin",
            (::ImmersedBoundaryMesh<2, 2>::ImmersedBoundaryLaminaIterator(ImmersedBoundaryMesh2_2::*)(bool)) &ImmersedBoundaryMesh2_2::GetLaminaIteratorBegin,
            " " , py::arg("skipDeletedLaminas") = true )
        .def(
            "GetLaminaIteratorEnd",
            (::ImmersedBoundaryMesh<2, 2>::ImmersedBoundaryLaminaIterator(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::GetLaminaIteratorEnd,
            " "  )
        .def(
            "GetNumNodes",
            (unsigned int(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetNumNodes,
            " "  )
        .def(
            "GetNumElements",
            (unsigned int(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetNumElements,
            " "  )
        .def(
            "GetNumAllElements",
            (unsigned int(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetNumAllElements,
            " "  )
        .def(
            "GetNumLaminas",
            (unsigned int(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetNumLaminas,
            " "  )
        .def(
            "GetNumGridPtsX",
            (unsigned int(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetNumGridPtsX,
            " "  )
        .def(
            "GetNumGridPtsY",
            (unsigned int(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetNumGridPtsY,
            " "  )
        .def(
            "GetCharacteristicNodeSpacing",
            (double(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetCharacteristicNodeSpacing,
            " "  )
        .def(
            "GetSpacingRatio",
            (double(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetSpacingRatio,
            " "  )
        .def(
            "GetMaxNodeIndex",
            (unsigned int(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetMaxNodeIndex,
            " "  )
        .def(
            "GetMaxElementIndex",
            (unsigned int(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetMaxElementIndex,
            " "  )
        .def(
            "GetMaxLaminaIndex",
            (unsigned int(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetMaxLaminaIndex,
            " "  )
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(ImmersedBoundaryMesh2_2::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &ImmersedBoundaryMesh2_2::GetVectorFromAtoB,
            " " , py::arg("rLocation1"), py::arg("rLocation2") )
        .def(
            "SetNode",
            (void(ImmersedBoundaryMesh2_2::*)(unsigned int, ::ChastePoint<2>)) &ImmersedBoundaryMesh2_2::SetNode,
            " " , py::arg("nodeIndex"), py::arg("point") )
        .def(
            "ConformToGeometry",
            (void(ImmersedBoundaryMesh2_2::*)(::boost::numeric::ublas::c_vector<double, 2> &)) &ImmersedBoundaryMesh2_2::ConformToGeometry,
            " " , py::arg("rLocation") )
        .def(
            "rGet2dVelocityGrids",
            (::boost::multi_array<double, 3> const &(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::rGet2dVelocityGrids,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetModifiable2dVelocityGrids",
            (::boost::multi_array<double, 3> &(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::rGetModifiable2dVelocityGrids,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "SetNumGridPtsX",
            (void(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::SetNumGridPtsX,
            " " , py::arg("meshPointsX") )
        .def(
            "SetNumGridPtsY",
            (void(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::SetNumGridPtsY,
            " " , py::arg("meshPointsY") )
        .def(
            "SetNumGridPtsXAndY",
            (void(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::SetNumGridPtsXAndY,
            " " , py::arg("numGridPts") )
        .def(
            "SetCharacteristicNodeSpacing",
            (void(ImmersedBoundaryMesh2_2::*)(double)) &ImmersedBoundaryMesh2_2::SetCharacteristicNodeSpacing,
            " " , py::arg("nodeSpacing") )
        .def(
            "AddNode",
            (unsigned int(ImmersedBoundaryMesh2_2::*)(::Node<2> *)) &ImmersedBoundaryMesh2_2::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "rGetElementFluidSources",
            (::std::vector<std::shared_ptr<FluidSource<2>>> &(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::rGetElementFluidSources,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetBalancingFluidSources",
            (::std::vector<std::shared_ptr<FluidSource<2>>> &(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::rGetBalancingFluidSources,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::GetNeighbouringNodeIndices,
            " " , py::arg("nodeIndex") )
        .def(
            "GetElement",
            (::ImmersedBoundaryElement<2, 2> *(ImmersedBoundaryMesh2_2::*)(unsigned int) const ) &ImmersedBoundaryMesh2_2::GetElement,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetLamina",
            (::ImmersedBoundaryElement<1, 2> *(ImmersedBoundaryMesh2_2::*)(unsigned int) const ) &ImmersedBoundaryMesh2_2::GetLamina,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetCentroidOfElement",
            (::boost::numeric::ublas::c_vector<double, 2>(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::GetCentroidOfElement,
            " " , py::arg("index") )
        .def(
            "ConstructFromMeshReader",
            (void(ImmersedBoundaryMesh2_2::*)(::AbstractMeshReader<2, 2> &)) &ImmersedBoundaryMesh2_2::ConstructFromMeshReader,
            " " , py::arg("rMeshReader") )
        .def(
            "Clear",
            (void(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::Clear,
            " "  )
        .def(
            "GetVolumeOfElement",
            (double(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::GetVolumeOfElement,
            " " , py::arg("index") )
        .def(
            "GetSurfaceAreaOfElement",
            (double(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::GetSurfaceAreaOfElement,
            " " , py::arg("index") )
        .def(
            "GetVoronoiSurfaceAreaOfElement",
            (double(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::GetVoronoiSurfaceAreaOfElement,
            " " , py::arg("elemIdx") )
        .def(
            "GetAverageNodeSpacingOfElement",
            (double(ImmersedBoundaryMesh2_2::*)(unsigned int, bool)) &ImmersedBoundaryMesh2_2::GetAverageNodeSpacingOfElement,
            " " , py::arg("index"), py::arg("recalculate") = true )
        .def(
            "GetAverageNodeSpacingOfLamina",
            (double(ImmersedBoundaryMesh2_2::*)(unsigned int, bool)) &ImmersedBoundaryMesh2_2::GetAverageNodeSpacingOfLamina,
            " " , py::arg("index"), py::arg("recalculate") = true )
        .def(
            "CalculateMomentsOfElement",
            (::boost::numeric::ublas::c_vector<double, 3>(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::CalculateMomentsOfElement,
            " " , py::arg("index") )
        .def(
            "GetElongationShapeFactorOfElement",
            (double(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::GetElongationShapeFactorOfElement,
            " " , py::arg("elementIndex") )
        .def(
            "GetTortuosityOfMesh",
            (double(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::GetTortuosityOfMesh,
            " "  )
        .def(
            "GetSkewnessOfElementMassDistributionAboutAxis",
            (double(ImmersedBoundaryMesh2_2::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 2>)) &ImmersedBoundaryMesh2_2::GetSkewnessOfElementMassDistributionAboutAxis,
            " " , py::arg("elemIndex"), py::arg("axis") )
        .def(
            "CalculateBoundingBoxOfElement",
            (::ChasteCuboid<2>(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::CalculateBoundingBoxOfElement,
            " " , py::arg("index") )
        .def(
            "GetShortAxisOfElement",
            (::boost::numeric::ublas::c_vector<double, 2>(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::GetShortAxisOfElement,
            " " , py::arg("index") )
        .def(
            "DivideElementAlongGivenAxis",
            (unsigned int(ImmersedBoundaryMesh2_2::*)(::ImmersedBoundaryElement<2, 2> *, ::boost::numeric::ublas::c_vector<double, 2>, bool)) &ImmersedBoundaryMesh2_2::DivideElementAlongGivenAxis,
            " " , py::arg("pElement"), py::arg("axisOfDivision"), py::arg("placeOriginalElementBelow") = false )
        .def(
            "DivideElementAlongShortAxis",
            (unsigned int(ImmersedBoundaryMesh2_2::*)(::ImmersedBoundaryElement<2, 2> *, bool)) &ImmersedBoundaryMesh2_2::DivideElementAlongShortAxis,
            " " , py::arg("pElement"), py::arg("placeOriginalElementBelow") = false )
        .def(
            "GetElementDivisionSpacing",
            (double(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::GetElementDivisionSpacing,
            " "  )
        .def(
            "SetElementDivisionSpacing",
            (void(ImmersedBoundaryMesh2_2::*)(double)) &ImmersedBoundaryMesh2_2::SetElementDivisionSpacing,
            " " , py::arg("elementDivisionSpacing") )
        .def(
            "GetNeighbourDist",
            (double(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetNeighbourDist,
            " "  )
        .def(
            "SetCellRearrangementThreshold",
            (void(ImmersedBoundaryMesh2_2::*)(double)) &ImmersedBoundaryMesh2_2::SetCellRearrangementThreshold,
            " " , py::arg("cellRearrangementThreshold") )
        .def(
            "GetCellRearrangementThreshold",
            (double(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::GetCellRearrangementThreshold,
            " "  )
        .def(
            "SetNeighbourDist",
            (void(ImmersedBoundaryMesh2_2::*)(double)) &ImmersedBoundaryMesh2_2::SetNeighbourDist,
            " " , py::arg("neighbourDist") )
        .def(
            "UpdateNodeLocationsVoronoiDiagramIfOutOfDate",
            (void(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::UpdateNodeLocationsVoronoiDiagramIfOutOfDate,
            " "  )
        .def(
            "ReMesh",
            (void(ImmersedBoundaryMesh2_2::*)(bool)) &ImmersedBoundaryMesh2_2::ReMesh,
            " " , py::arg("randomOrder") = false )
        .def(
            "ReMeshElement",
            (void(ImmersedBoundaryMesh2_2::*)(::ImmersedBoundaryElement<2, 2> *, bool)) &ImmersedBoundaryMesh2_2::ReMeshElement,
            " " , py::arg("pElement"), py::arg("randomOrder") )
        .def(
            "ReMeshLamina",
            (void(ImmersedBoundaryMesh2_2::*)(::ImmersedBoundaryElement<1, 2> *, bool)) &ImmersedBoundaryMesh2_2::ReMeshLamina,
            " " , py::arg("pLamina"), py::arg("randomOrder") )
        .def(
            "NodesInDifferentElementOrLamina",
            (bool(ImmersedBoundaryMesh2_2::*)(::Node<2> *, ::Node<2> *)) &ImmersedBoundaryMesh2_2::NodesInDifferentElementOrLamina,
            " " , py::arg("pNodeA"), py::arg("pNodeB") )
        .def(
            "GetNeighbouringElementIndices",
            (::std::set<unsigned int>(ImmersedBoundaryMesh2_2::*)(unsigned int)) &ImmersedBoundaryMesh2_2::GetNeighbouringElementIndices,
            " " , py::arg("elemIdx") )
        .def(
            "CalculateLengthOfVoronoiEdge",
            (double(ImmersedBoundaryMesh2_2::*)(::boost::polygon::voronoi_diagram<double>::edge_type const &)) &ImmersedBoundaryMesh2_2::CalculateLengthOfVoronoiEdge,
            " " , py::arg("rEdge") )
        .def(
            "GetPolygonDistribution",
            (::std::array<unsigned int, 13>(ImmersedBoundaryMesh2_2::*)()) &ImmersedBoundaryMesh2_2::GetPolygonDistribution,
            " "  )
        .def(
            "rGetNodeLocationsVoronoiDiagram",
            (::boost::polygon::voronoi_diagram<double> const &(ImmersedBoundaryMesh2_2::*)(bool)) &ImmersedBoundaryMesh2_2::rGetNodeLocationsVoronoiDiagram,
            " " , py::arg("update") = true , py::return_value_policy::reference_internal)
        .def(
            "GetVoronoiCellIdsIndexedByNodeIndex",
            (::std::vector<unsigned int> const &(ImmersedBoundaryMesh2_2::*)() const ) &ImmersedBoundaryMesh2_2::GetVoronoiCellIdsIndexedByNodeIndex,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "ScaleUpToVoronoiCoordinate",
            (int(ImmersedBoundaryMesh2_2::*)(double) const ) &ImmersedBoundaryMesh2_2::ScaleUpToVoronoiCoordinate,
            " " , py::arg("location") )
        .def(
            "ScaleDistanceDownFromVoronoi",
            (double(ImmersedBoundaryMesh2_2::*)(double const) const ) &ImmersedBoundaryMesh2_2::ScaleDistanceDownFromVoronoi,
            " " , py::arg("distance") )
    ;
}
