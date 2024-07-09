#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

#include "HoneycombVertexMeshGenerator.cppwg.hpp"

namespace py = pybind11;
typedef HoneycombVertexMeshGenerator HoneycombVertexMeshGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::shared_ptr<MutableVertexMesh<2, 2>> _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_;

class HoneycombVertexMeshGenerator_Overrides : public HoneycombVertexMeshGenerator{
    public:
    using HoneycombVertexMeshGenerator::HoneycombVertexMeshGenerator;
    ::boost::shared_ptr<MutableVertexMesh<2, 2>> GetMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_,
            HoneycombVertexMeshGenerator,
            GetMesh,
            );
    }

};
void register_HoneycombVertexMeshGenerator_class(py::module &m){
py::class_<HoneycombVertexMeshGenerator , HoneycombVertexMeshGenerator_Overrides , boost::shared_ptr<HoneycombVertexMeshGenerator >   >(m, "HoneycombVertexMeshGenerator")
        .def(py::init<unsigned int, unsigned int, bool, double, double, double >(), py::arg("numElementsAcross"), py::arg("numElementsUp"), py::arg("isFlatBottom") = false, py::arg("cellRearrangementThreshold") = 0.01, py::arg("t2Threshold") = 0.001, py::arg("elementArea") = 0.5 * sqrt(3.))
        .def(py::init< >())
        .def(
            "GetMesh",
            (::boost::shared_ptr<MutableVertexMesh<2, 2>>(HoneycombVertexMeshGenerator::*)()) &HoneycombVertexMeshGenerator::GetMesh,
            " "  )
    ;
}
