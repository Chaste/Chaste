#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VoronoiVertexMeshGenerator.hpp"

#include "VoronoiVertexMeshGenerator.cppwg.hpp"

namespace py = pybind11;
typedef VoronoiVertexMeshGenerator VoronoiVertexMeshGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::shared_ptr<MutableVertexMesh<2, 2>> _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_;
typedef ::boost::shared_ptr<MutableVertexMesh<2, 2>> _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_;
typedef ::boost::shared_ptr<Toroidal2dVertexMesh> _boost_shared_ptr_lt_Toroidal2dVertexMesh_gt_;

class VoronoiVertexMeshGenerator_Overrides : public VoronoiVertexMeshGenerator{
    public:
    using VoronoiVertexMeshGenerator::VoronoiVertexMeshGenerator;
    ::boost::shared_ptr<MutableVertexMesh<2, 2>> GetMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_,
            VoronoiVertexMeshGenerator,
            GetMesh,
            );
    }
    ::boost::shared_ptr<MutableVertexMesh<2, 2>> GetMeshAfterReMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_,
            VoronoiVertexMeshGenerator,
            GetMeshAfterReMesh,
            );
    }
    ::boost::shared_ptr<Toroidal2dVertexMesh> GetToroidalMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_Toroidal2dVertexMesh_gt_,
            VoronoiVertexMeshGenerator,
            GetToroidalMesh,
            );
    }

};
void register_VoronoiVertexMeshGenerator_class(py::module &m){
py::class_<VoronoiVertexMeshGenerator , VoronoiVertexMeshGenerator_Overrides , boost::shared_ptr<VoronoiVertexMeshGenerator >   >(m, "VoronoiVertexMeshGenerator")
        .def(py::init<unsigned int, unsigned int, unsigned int, double >(), py::arg("numElementsX"), py::arg("numElementsY"), py::arg("numRelaxationSteps"), py::arg("elementTargetArea") = 1.)
        .def(py::init< >())
        .def(
            "GenerateVoronoiMesh",
            (void(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GenerateVoronoiMesh,
            " "  )
        .def(
            "GetMesh",
            (::boost::shared_ptr<MutableVertexMesh<2, 2>>(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetMesh,
            " "  )
        .def(
            "GetMeshAfterReMesh",
            (::boost::shared_ptr<MutableVertexMesh<2, 2>>(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetMeshAfterReMesh,
            " "  )
        .def(
            "GetToroidalMesh",
            (::boost::shared_ptr<Toroidal2dVertexMesh>(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetToroidalMesh,
            " "  )
        .def(
            "GetPolygonDistribution",
            (::std::vector<double>(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetPolygonDistribution,
            " "  )
        .def(
            "GetAreaCoefficientOfVariation",
            (double(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetAreaCoefficientOfVariation,
            " "  )
        .def(
            "RefreshSeedsAndRegenerateMesh",
            (void(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::RefreshSeedsAndRegenerateMesh,
            " "  )
        .def(
            "SetMaxExpectedNumSidesPerPolygon",
            (void(VoronoiVertexMeshGenerator::*)(unsigned int)) &VoronoiVertexMeshGenerator::SetMaxExpectedNumSidesPerPolygon,
            " " , py::arg("maxExpectedNumSidesPerPolygon") )
        .def(
            "GetMaxExpectedNumSidesPerPolygon",
            (unsigned int(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetMaxExpectedNumSidesPerPolygon,
            " "  )
    ;
}
