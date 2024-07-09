#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "MutableElement.hpp"

#include "MutableElement3_3.cppwg.hpp"

namespace py = pybind11;
typedef MutableElement<3,3 > MutableElement3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class MutableElement3_3_Overrides : public MutableElement3_3{
    public:
    using MutableElement3_3::MutableElement;
    void RegisterWithNodes() override {
        PYBIND11_OVERRIDE(
            void,
            MutableElement3_3,
            RegisterWithNodes,
            );
    }
    void MarkAsDeleted() override {
        PYBIND11_OVERRIDE(
            void,
            MutableElement3_3,
            MarkAsDeleted,
            );
    }
    void UpdateNode(unsigned int const & rIndex, ::Node<3> * pNode) override {
        PYBIND11_OVERRIDE(
            void,
            MutableElement3_3,
            UpdateNode,
                    rIndex,
        pNode);
    }
    bool IsElementOnBoundary() const  override {
        PYBIND11_OVERRIDE(
            bool,
            MutableElement3_3,
            IsElementOnBoundary,
            );
    }

};
void register_MutableElement3_3_class(py::module &m){
py::class_<MutableElement3_3 , MutableElement3_3_Overrides , boost::shared_ptr<MutableElement3_3 >  , AbstractElement<3, 3>  >(m, "MutableElement3_3")
        .def(py::init<unsigned int >(), py::arg("index"))
        .def(py::init<unsigned int, ::std::vector<Node<3> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(
            "RegisterWithNodes",
            (void(MutableElement3_3::*)()) &MutableElement3_3::RegisterWithNodes,
            " "  )
        .def(
            "MarkAsDeleted",
            (void(MutableElement3_3::*)()) &MutableElement3_3::MarkAsDeleted,
            " "  )
        .def(
            "ResetIndex",
            (void(MutableElement3_3::*)(unsigned int)) &MutableElement3_3::ResetIndex,
            " " , py::arg("index") )
        .def(
            "UpdateNode",
            (void(MutableElement3_3::*)(unsigned int const &, ::Node<3> *)) &MutableElement3_3::UpdateNode,
            " " , py::arg("rIndex"), py::arg("pNode") )
        .def(
            "DeleteNode",
            (void(MutableElement3_3::*)(unsigned int const &)) &MutableElement3_3::DeleteNode,
            " " , py::arg("rIndex") )
        .def(
            "AddNode",
            (void(MutableElement3_3::*)(::Node<3> *, unsigned int const &)) &MutableElement3_3::AddNode,
            " " , py::arg("pNode"), py::arg("rIndex") )
        .def(
            "GetNodeLocalIndex",
            (unsigned int(MutableElement3_3::*)(unsigned int) const ) &MutableElement3_3::GetNodeLocalIndex,
            " " , py::arg("globalIndex") )
        .def(
            "RegisterWithEdges",
            (void(MutableElement3_3::*)()) &MutableElement3_3::RegisterWithEdges,
            " "  )
        .def(
            "RebuildEdges",
            (void(MutableElement3_3::*)()) &MutableElement3_3::RebuildEdges,
            " "  )
        .def(
            "IsElementOnBoundary",
            (bool(MutableElement3_3::*)() const ) &MutableElement3_3::IsElementOnBoundary,
            " "  )
        .def(
            "SetEdgeHelper",
            (void(MutableElement3_3::*)(::EdgeHelper<3> *)) &MutableElement3_3::SetEdgeHelper,
            " " , py::arg("pEdgeHelper") )
        .def(
            "ClearEdges",
            (void(MutableElement3_3::*)()) &MutableElement3_3::ClearEdges,
            " "  )
        .def(
            "BuildEdges",
            (void(MutableElement3_3::*)()) &MutableElement3_3::BuildEdges,
            " "  )
        .def(
            "GetEdgeGlobalIndex",
            (unsigned int(MutableElement3_3::*)(unsigned int) const ) &MutableElement3_3::GetEdgeGlobalIndex,
            " " , py::arg("localIndex") )
        .def(
            "GetEdge",
            (::Edge<3> *(MutableElement3_3::*)(unsigned int) const ) &MutableElement3_3::GetEdge,
            " " , py::arg("localIndex") , py::return_value_policy::reference)
        .def(
            "GetNumEdges",
            (unsigned int(MutableElement3_3::*)() const ) &MutableElement3_3::GetNumEdges,
            " "  )
        .def(
            "GetNeighbouringElementAtEdgeIndex",
            (::std::set<unsigned int>(MutableElement3_3::*)(unsigned int)) &MutableElement3_3::GetNeighbouringElementAtEdgeIndex,
            " " , py::arg("localIndex") )
        .def(
            "ContainsEdge",
            (bool(MutableElement3_3::*)(::Edge<3> const *) const ) &MutableElement3_3::ContainsEdge,
            " " , py::arg("pEdge") )
        .def(
            "GetLocalEdgeIndex",
            (long int(MutableElement3_3::*)(::Edge<3> const *) const ) &MutableElement3_3::GetLocalEdgeIndex,
            " " , py::arg("pEdge") )
    ;
}
