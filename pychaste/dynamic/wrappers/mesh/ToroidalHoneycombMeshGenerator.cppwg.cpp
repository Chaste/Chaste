#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ToroidalHoneycombMeshGenerator.hpp"

#include "ToroidalHoneycombMeshGenerator.cppwg.hpp"

namespace py = pybind11;
typedef ToroidalHoneycombMeshGenerator ToroidalHoneycombMeshGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::shared_ptr<MutableMesh<2, 2>> _boost_shared_ptr_lt_MutableMesh_lt_2_2_gt__gt_;

class ToroidalHoneycombMeshGenerator_Overrides : public ToroidalHoneycombMeshGenerator{
    public:
    using ToroidalHoneycombMeshGenerator::ToroidalHoneycombMeshGenerator;
    ::boost::shared_ptr<MutableMesh<2, 2>> GetMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_MutableMesh_lt_2_2_gt__gt_,
            ToroidalHoneycombMeshGenerator,
            GetMesh,
            );
    }

};
void register_ToroidalHoneycombMeshGenerator_class(py::module &m){
py::class_<ToroidalHoneycombMeshGenerator , ToroidalHoneycombMeshGenerator_Overrides , boost::shared_ptr<ToroidalHoneycombMeshGenerator >  , HoneycombMeshGenerator  >(m, "ToroidalHoneycombMeshGenerator")
        .def(py::init<unsigned int, unsigned int, double, double >(), py::arg("numNodesAlongWidth"), py::arg("numNodesAlongDepth"), py::arg("widthScaleFactor") = 1., py::arg("depthScaleFactor") = 1.)
        .def(
            "GetMesh",
            (::boost::shared_ptr<MutableMesh<2, 2>>(ToroidalHoneycombMeshGenerator::*)()) &ToroidalHoneycombMeshGenerator::GetMesh,
            " "  )
        .def(
            "GetToroidalMesh",
            (::boost::shared_ptr<Toroidal2dMesh>(ToroidalHoneycombMeshGenerator::*)()) &ToroidalHoneycombMeshGenerator::GetToroidalMesh,
            " "  )
    ;
}
