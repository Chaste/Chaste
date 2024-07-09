#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"

#include "CylindricalHoneycombMeshGenerator.cppwg.hpp"

namespace py = pybind11;
typedef CylindricalHoneycombMeshGenerator CylindricalHoneycombMeshGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::shared_ptr<MutableMesh<2, 2>> _boost_shared_ptr_lt_MutableMesh_lt_2_2_gt__gt_;

class CylindricalHoneycombMeshGenerator_Overrides : public CylindricalHoneycombMeshGenerator{
    public:
    using CylindricalHoneycombMeshGenerator::CylindricalHoneycombMeshGenerator;
    ::boost::shared_ptr<MutableMesh<2, 2>> GetMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_MutableMesh_lt_2_2_gt__gt_,
            CylindricalHoneycombMeshGenerator,
            GetMesh,
            );
    }

};
void register_CylindricalHoneycombMeshGenerator_class(py::module &m){
py::class_<CylindricalHoneycombMeshGenerator , CylindricalHoneycombMeshGenerator_Overrides , boost::shared_ptr<CylindricalHoneycombMeshGenerator >  , HoneycombMeshGenerator  >(m, "CylindricalHoneycombMeshGenerator")
        .def(py::init<unsigned int, unsigned int, unsigned int, double >(), py::arg("numNodesAlongWidth"), py::arg("numNodesAlongLength"), py::arg("ghosts") = 3, py::arg("scaleFactor") = 1.)
        .def(
            "GetMesh",
            (::boost::shared_ptr<MutableMesh<2, 2>>(CylindricalHoneycombMeshGenerator::*)()) &CylindricalHoneycombMeshGenerator::GetMesh,
            " "  )
        .def(
            "GetCylindricalMesh",
            (::boost::shared_ptr<Cylindrical2dMesh>(CylindricalHoneycombMeshGenerator::*)()) &CylindricalHoneycombMeshGenerator::GetCylindricalMesh,
            " "  )
    ;
}
