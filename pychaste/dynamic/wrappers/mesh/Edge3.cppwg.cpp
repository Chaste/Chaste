#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Edge.hpp"

#include "Edge3.cppwg.hpp"

namespace py = pybind11;
typedef Edge<3 > Edge3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_Edge3_class(py::module &m){
py::class_<Edge3  , boost::shared_ptr<Edge3 >   >(m, "Edge3")
        .def(py::init<unsigned int >(), py::arg("index"))
        .def(py::init<unsigned int, ::Node<3> *, ::Node<3> * >(), py::arg("index"), py::arg("pNodeA"), py::arg("pNodeB"))
        .def_static(
            "GenerateMapIndex",
            (::std::pair<unsigned int, unsigned int>(*)(unsigned int, unsigned int)) &Edge3::GenerateMapIndex,
            " " , py::arg("index1"), py::arg("index2") )
        .def(
            "MarkAsDeleted",
            (void(Edge3::*)()) &Edge3::MarkAsDeleted,
            " "  )
        .def(
            "IsDeleted",
            (bool(Edge3::*)()) &Edge3::IsDeleted,
            " "  )
        .def(
            "SetIndex",
            (void(Edge3::*)(unsigned int)) &Edge3::SetIndex,
            " " , py::arg("index") )
        .def(
            "GetIndex",
            (unsigned int(Edge3::*)() const ) &Edge3::GetIndex,
            " "  )
        .def(
            "GetMapIndex",
            (::std::pair<unsigned int, unsigned int>(Edge3::*)()) &Edge3::GetMapIndex,
            " "  )
        .def(
            "RemoveNodes",
            (void(Edge3::*)()) &Edge3::RemoveNodes,
            " "  )
        .def(
            "SetNodes",
            (void(Edge3::*)(::Node<3> *, ::Node<3> *)) &Edge3::SetNodes,
            " " , py::arg("pNodeA"), py::arg("pNodeB") )
        .def(
            "ReplaceNode",
            (void(Edge3::*)(::Node<3> *, ::Node<3> *)) &Edge3::ReplaceNode,
            " " , py::arg("pOldNode"), py::arg("pNewNode") )
        .def(
            "GetNode",
            (::Node<3> *(Edge3::*)(unsigned int) const ) &Edge3::GetNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNumNodes",
            (unsigned int(Edge3::*)()) &Edge3::GetNumNodes,
            " "  )
        .def(
            "ContainsNode",
            (bool(Edge3::*)(::Node<3> *) const ) &Edge3::ContainsNode,
            " " , py::arg("pNode") )
        .def(
            "rGetCentreLocation",
            (::boost::numeric::ublas::c_vector<double, 3>(Edge3::*)()) &Edge3::rGetCentreLocation,
            " "  )
        .def(
            "rGetLength",
            (double(Edge3::*)()) &Edge3::rGetLength,
            " "  )
        .def(
            "GetOtherElements",
            (::std::set<unsigned int>(Edge3::*)(unsigned int)) &Edge3::GetOtherElements,
            " " , py::arg("elementIndex") )
        .def(
            "AddElement",
            (void(Edge3::*)(unsigned int)) &Edge3::AddElement,
            " " , py::arg("elementIndex") )
        .def(
            "RemoveElement",
            (void(Edge3::*)(unsigned int)) &Edge3::RemoveElement,
            " " , py::arg("elementIndex") )
        .def(
            "GetNeighbouringElementIndices",
            (::std::set<unsigned int>(Edge3::*)()) &Edge3::GetNeighbouringElementIndices,
            " "  )
        .def(
            "GetNumElements",
            (unsigned int(Edge3::*)()) &Edge3::GetNumElements,
            " "  )
        .def(
            "IsBoundaryEdge",
            (bool(Edge3::*)() const ) &Edge3::IsBoundaryEdge,
            " "  )
    ;
}
