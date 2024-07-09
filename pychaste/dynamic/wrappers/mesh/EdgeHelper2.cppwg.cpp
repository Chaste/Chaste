#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "EdgeHelper.hpp"

#include "EdgeHelper2.cppwg.hpp"

namespace py = pybind11;
typedef EdgeHelper<2 > EdgeHelper2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_EdgeHelper2_class(py::module &m){
py::class_<EdgeHelper2  , boost::shared_ptr<EdgeHelper2 >   >(m, "EdgeHelper2")
        .def(py::init< >())
        .def(
            "GetEdgeFromNodes",
            (::Edge<2> *(EdgeHelper2::*)(::Node<2> *, ::Node<2> *)) &EdgeHelper2::GetEdgeFromNodes,
            " " , py::arg("pNodeA"), py::arg("pNodeB") , py::return_value_policy::reference)
        .def(
            "GetEdgeFromNodes",
            (::Edge<2> *(EdgeHelper2::*)(unsigned int, ::Node<2> *, ::Node<2> *)) &EdgeHelper2::GetEdgeFromNodes,
            " " , py::arg("elementIndex"), py::arg("pNodeA"), py::arg("pNodeB") , py::return_value_policy::reference)
        .def(
            "GetEdge",
            (::Edge<2> *(EdgeHelper2::*)(unsigned int) const ) &EdgeHelper2::GetEdge,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "RemoveDeletedEdges",
            (void(EdgeHelper2::*)()) &EdgeHelper2::RemoveDeletedEdges,
            " "  )
        .def(
            "GetNumEdges",
            (unsigned int(EdgeHelper2::*)() const ) &EdgeHelper2::GetNumEdges,
            " "  )
    ;
}
