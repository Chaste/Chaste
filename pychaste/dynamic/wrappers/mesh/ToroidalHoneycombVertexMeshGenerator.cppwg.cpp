#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"

#include "ToroidalHoneycombVertexMeshGenerator.cppwg.hpp"

namespace py = pybind11;
typedef ToroidalHoneycombVertexMeshGenerator ToroidalHoneycombVertexMeshGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::shared_ptr<MutableVertexMesh<2, 2>> _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_;

class ToroidalHoneycombVertexMeshGenerator_Overrides : public ToroidalHoneycombVertexMeshGenerator{
    public:
    using ToroidalHoneycombVertexMeshGenerator::ToroidalHoneycombVertexMeshGenerator;
    ::boost::shared_ptr<MutableVertexMesh<2, 2>> GetMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_,
            ToroidalHoneycombVertexMeshGenerator,
            GetMesh,
            );
    }

};
void register_ToroidalHoneycombVertexMeshGenerator_class(py::module &m){
py::class_<ToroidalHoneycombVertexMeshGenerator , ToroidalHoneycombVertexMeshGenerator_Overrides , boost::shared_ptr<ToroidalHoneycombVertexMeshGenerator >   >(m, "ToroidalHoneycombVertexMeshGenerator")
        .def(py::init<unsigned int, unsigned int, double, double >(), py::arg("numElementsAcross"), py::arg("numElementsUp"), py::arg("cellRearrangementThreshold") = 0.01, py::arg("t2Threshold") = 0.001)
        .def(
            "GetMesh",
            (::boost::shared_ptr<MutableVertexMesh<2, 2>>(ToroidalHoneycombVertexMeshGenerator::*)()) &ToroidalHoneycombVertexMeshGenerator::GetMesh,
            " "  )
        .def(
            "GetToroidalMesh",
            (::boost::shared_ptr<Toroidal2dVertexMesh>(ToroidalHoneycombVertexMeshGenerator::*)()) &ToroidalHoneycombVertexMeshGenerator::GetToroidalMesh,
            " "  )
    ;
}
