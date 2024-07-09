#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "MutableElement.hpp"

#include "MutableElement2_2.cppwg.hpp"

namespace py = pybind11;
typedef MutableElement<2,2 > MutableElement2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class MutableElement2_2_Overrides : public MutableElement2_2{
    public:
    using MutableElement2_2::MutableElement;
    void RegisterWithNodes() override {
        PYBIND11_OVERRIDE(
            void,
            MutableElement2_2,
            RegisterWithNodes,
            );
    }
    void MarkAsDeleted() override {
        PYBIND11_OVERRIDE(
            void,
            MutableElement2_2,
            MarkAsDeleted,
            );
    }
    void UpdateNode(unsigned int const & rIndex, ::Node<2> * pNode) override {
        PYBIND11_OVERRIDE(
            void,
            MutableElement2_2,
            UpdateNode,
                    rIndex,
        pNode);
    }
    bool IsElementOnBoundary() const  override {
        PYBIND11_OVERRIDE(
            bool,
            MutableElement2_2,
            IsElementOnBoundary,
            );
    }

};
void register_MutableElement2_2_class(py::module &m){
py::class_<MutableElement2_2 , MutableElement2_2_Overrides , boost::shared_ptr<MutableElement2_2 >  , AbstractElement<2, 2>  >(m, "MutableElement2_2")
        .def(py::init<unsigned int >(), py::arg("index"))
        .def(py::init<unsigned int, ::std::vector<Node<2> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(
            "RegisterWithNodes",
            (void(MutableElement2_2::*)()) &MutableElement2_2::RegisterWithNodes,
            " "  )
        .def(
            "MarkAsDeleted",
            (void(MutableElement2_2::*)()) &MutableElement2_2::MarkAsDeleted,
            " "  )
        .def(
            "ResetIndex",
            (void(MutableElement2_2::*)(unsigned int)) &MutableElement2_2::ResetIndex,
            " " , py::arg("index") )
        .def(
            "UpdateNode",
            (void(MutableElement2_2::*)(unsigned int const &, ::Node<2> *)) &MutableElement2_2::UpdateNode,
            " " , py::arg("rIndex"), py::arg("pNode") )
        .def(
            "DeleteNode",
            (void(MutableElement2_2::*)(unsigned int const &)) &MutableElement2_2::DeleteNode,
            " " , py::arg("rIndex") )
        .def(
            "AddNode",
            (void(MutableElement2_2::*)(::Node<2> *, unsigned int const &)) &MutableElement2_2::AddNode,
            " " , py::arg("pNode"), py::arg("rIndex") )
        .def(
            "GetNodeLocalIndex",
            (unsigned int(MutableElement2_2::*)(unsigned int) const ) &MutableElement2_2::GetNodeLocalIndex,
            " " , py::arg("globalIndex") )
        .def(
            "RegisterWithEdges",
            (void(MutableElement2_2::*)()) &MutableElement2_2::RegisterWithEdges,
            " "  )
        .def(
            "RebuildEdges",
            (void(MutableElement2_2::*)()) &MutableElement2_2::RebuildEdges,
            " "  )
        .def(
            "IsElementOnBoundary",
            (bool(MutableElement2_2::*)() const ) &MutableElement2_2::IsElementOnBoundary,
            " "  )
        .def(
            "SetEdgeHelper",
            (void(MutableElement2_2::*)(::EdgeHelper<2> *)) &MutableElement2_2::SetEdgeHelper,
            " " , py::arg("pEdgeHelper") )
        .def(
            "ClearEdges",
            (void(MutableElement2_2::*)()) &MutableElement2_2::ClearEdges,
            " "  )
        .def(
            "BuildEdges",
            (void(MutableElement2_2::*)()) &MutableElement2_2::BuildEdges,
            " "  )
        .def(
            "GetEdgeGlobalIndex",
            (unsigned int(MutableElement2_2::*)(unsigned int) const ) &MutableElement2_2::GetEdgeGlobalIndex,
            " " , py::arg("localIndex") )
        .def(
            "GetEdge",
            (::Edge<2> *(MutableElement2_2::*)(unsigned int) const ) &MutableElement2_2::GetEdge,
            " " , py::arg("localIndex") , py::return_value_policy::reference)
        .def(
            "GetNumEdges",
            (unsigned int(MutableElement2_2::*)() const ) &MutableElement2_2::GetNumEdges,
            " "  )
        .def(
            "GetNeighbouringElementAtEdgeIndex",
            (::std::set<unsigned int>(MutableElement2_2::*)(unsigned int)) &MutableElement2_2::GetNeighbouringElementAtEdgeIndex,
            " " , py::arg("localIndex") )
        .def(
            "ContainsEdge",
            (bool(MutableElement2_2::*)(::Edge<2> const *) const ) &MutableElement2_2::ContainsEdge,
            " " , py::arg("pEdge") )
        .def(
            "GetLocalEdgeIndex",
            (long int(MutableElement2_2::*)(::Edge<2> const *) const ) &MutableElement2_2::GetLocalEdgeIndex,
            " " , py::arg("pEdge") )
    ;
}
