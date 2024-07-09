#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "MutableElement.hpp"

#include "MutableElement1_2.cppwg.hpp"

namespace py = pybind11;
typedef MutableElement<1,2 > MutableElement1_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class MutableElement1_2_Overrides : public MutableElement1_2{
    public:
    using MutableElement1_2::MutableElement;
    void UpdateNode(unsigned int const & rIndex, ::Node<2> * pNode) override {
        PYBIND11_OVERRIDE(
            void,
            MutableElement1_2,
            UpdateNode,
                    rIndex,
        pNode);
    }
    void RegisterWithNodes() override {
        PYBIND11_OVERRIDE(
            void,
            MutableElement1_2,
            RegisterWithNodes,
            );
    }
    void MarkAsDeleted() override {
        PYBIND11_OVERRIDE(
            void,
            MutableElement1_2,
            MarkAsDeleted,
            );
    }
    bool IsElementOnBoundary() const  override {
        PYBIND11_OVERRIDE(
            bool,
            MutableElement1_2,
            IsElementOnBoundary,
            );
    }

};
void register_MutableElement1_2_class(py::module &m){
py::class_<MutableElement1_2 , MutableElement1_2_Overrides , boost::shared_ptr<MutableElement1_2 >  , AbstractElement<1, 2>  >(m, "MutableElement1_2")
        .def(py::init<unsigned int, ::std::vector<Node<2> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(
            "UpdateNode",
            (void(MutableElement1_2::*)(unsigned int const &, ::Node<2> *)) &MutableElement1_2::UpdateNode,
            " " , py::arg("rIndex"), py::arg("pNode") )
        .def(
            "RegisterWithNodes",
            (void(MutableElement1_2::*)()) &MutableElement1_2::RegisterWithNodes,
            " "  )
        .def(
            "MarkAsDeleted",
            (void(MutableElement1_2::*)()) &MutableElement1_2::MarkAsDeleted,
            " "  )
        .def(
            "ResetIndex",
            (void(MutableElement1_2::*)(unsigned int)) &MutableElement1_2::ResetIndex,
            " " , py::arg("index") )
        .def(
            "DeleteNode",
            (void(MutableElement1_2::*)(unsigned int const &)) &MutableElement1_2::DeleteNode,
            " " , py::arg("rIndex") )
        .def(
            "AddNode",
            (void(MutableElement1_2::*)(::Node<2> *, unsigned int const &)) &MutableElement1_2::AddNode,
            " " , py::arg("pNode"), py::arg("rIndex") )
        .def(
            "GetEdge",
            (::Edge<2> *(MutableElement1_2::*)(unsigned int) const ) &MutableElement1_2::GetEdge,
            " " , py::arg("localIndex") , py::return_value_policy::reference)
        .def(
            "ContainsEdge",
            (bool(MutableElement1_2::*)(::Edge<2> const *) const ) &MutableElement1_2::ContainsEdge,
            " " , py::arg("pEdge") )
        .def(
            "GetNumEdges",
            (unsigned int(MutableElement1_2::*)() const ) &MutableElement1_2::GetNumEdges,
            " "  )
        .def(
            "SetEdgeHelper",
            (void(MutableElement1_2::*)(::EdgeHelper<2> *)) &MutableElement1_2::SetEdgeHelper,
            " " , py::arg("pEdgeHelper") )
        .def(
            "BuildEdges",
            (void(MutableElement1_2::*)()) &MutableElement1_2::BuildEdges,
            " "  )
        .def(
            "ClearEdges",
            (void(MutableElement1_2::*)()) &MutableElement1_2::ClearEdges,
            " "  )
        .def(
            "GetEdgeGlobalIndex",
            (unsigned int(MutableElement1_2::*)(unsigned int) const ) &MutableElement1_2::GetEdgeGlobalIndex,
            " " , py::arg("localIndex") )
        .def(
            "GetNeighbouringElementAtEdgeIndex",
            (::std::set<unsigned int>(MutableElement1_2::*)(unsigned int)) &MutableElement1_2::GetNeighbouringElementAtEdgeIndex,
            " " , py::arg("localIndex") )
        .def(
            "GetNodeLocalIndex",
            (unsigned int(MutableElement1_2::*)(unsigned int) const ) &MutableElement1_2::GetNodeLocalIndex,
            " " , py::arg("globalIndex") )
        .def(
            "RegisterWithEdges",
            (void(MutableElement1_2::*)()) &MutableElement1_2::RegisterWithEdges,
            " "  )
        .def(
            "RebuildEdges",
            (void(MutableElement1_2::*)()) &MutableElement1_2::RebuildEdges,
            " "  )
        .def(
            "IsElementOnBoundary",
            (bool(MutableElement1_2::*)() const ) &MutableElement1_2::IsElementOnBoundary,
            " "  )
        .def(
            "GetLocalEdgeIndex",
            (long int(MutableElement1_2::*)(::Edge<2> const *) const ) &MutableElement1_2::GetLocalEdgeIndex,
            " " , py::arg("pEdge") )
    ;
}
