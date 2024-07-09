#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "EdgeHelper.hpp"

#include "EdgeHelper3.cppwg.hpp"

namespace py = pybind11;
typedef EdgeHelper<3 > EdgeHelper3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_EdgeHelper3_class(py::module &m){
py::class_<EdgeHelper3  , boost::shared_ptr<EdgeHelper3 >   >(m, "EdgeHelper3")
        .def(py::init< >())
        .def(
            "GetEdgeFromNodes",
            (::Edge<3> *(EdgeHelper3::*)(::Node<3> *, ::Node<3> *)) &EdgeHelper3::GetEdgeFromNodes,
            " " , py::arg("pNodeA"), py::arg("pNodeB") , py::return_value_policy::reference)
        .def(
            "GetEdgeFromNodes",
            (::Edge<3> *(EdgeHelper3::*)(unsigned int, ::Node<3> *, ::Node<3> *)) &EdgeHelper3::GetEdgeFromNodes,
            " " , py::arg("elementIndex"), py::arg("pNodeA"), py::arg("pNodeB") , py::return_value_policy::reference)
        .def(
            "GetEdge",
            (::Edge<3> *(EdgeHelper3::*)(unsigned int) const ) &EdgeHelper3::GetEdge,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "RemoveDeletedEdges",
            (void(EdgeHelper3::*)()) &EdgeHelper3::RemoveDeletedEdges,
            " "  )
        .def(
            "GetNumEdges",
            (unsigned int(EdgeHelper3::*)() const ) &EdgeHelper3::GetNumEdges,
            " "  )
    ;
}
