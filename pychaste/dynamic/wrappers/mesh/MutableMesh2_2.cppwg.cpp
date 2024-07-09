#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "MutableMesh.hpp"

#include "MutableMesh2_2.cppwg.hpp"

namespace py = pybind11;
typedef MutableMesh<2,2 > MutableMesh2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class MutableMesh2_2_Overrides : public MutableMesh2_2{
    public:
    using MutableMesh2_2::MutableMesh;
    void Clear() override {
        PYBIND11_OVERRIDE(
            void,
            MutableMesh2_2,
            Clear,
            );
    }
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableMesh2_2,
            GetNumNodes,
            );
    }
    unsigned int GetNumElements() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableMesh2_2,
            GetNumElements,
            );
    }
    unsigned int GetNumBoundaryElements() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableMesh2_2,
            GetNumBoundaryElements,
            );
    }
    unsigned int AddNode(::Node<2> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableMesh2_2,
            AddNode,
                    pNewNode);
    }
    void SetNode(unsigned int index, ::ChastePoint<2> point, bool concreteMove) override {
        PYBIND11_OVERRIDE(
            void,
            MutableMesh2_2,
            SetNode,
                    index,
        point,
        concreteMove);
    }
    void DeleteNode(unsigned int index) override {
        PYBIND11_OVERRIDE(
            void,
            MutableMesh2_2,
            DeleteNode,
                    index);
    }
    void DeleteElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            void,
            MutableMesh2_2,
            DeleteElement,
                    index);
    }
    void ReMesh(::NodeMap & map) override {
        PYBIND11_OVERRIDE(
            void,
            MutableMesh2_2,
            ReMesh,
                    map);
    }

};
void register_MutableMesh2_2_class(py::module &m){
py::class_<MutableMesh2_2 , MutableMesh2_2_Overrides , boost::shared_ptr<MutableMesh2_2 >  , TetrahedralMesh<2, 2>  >(m, "MutableMesh2_2")
        .def(py::init< >())
        .def(py::init<::std::vector<Node<2> *> >(), py::arg("nodes"))
        .def(
            "Clear",
            (void(MutableMesh2_2::*)()) &MutableMesh2_2::Clear,
            " "  )
        .def(
            "GetNumNodes",
            (unsigned int(MutableMesh2_2::*)() const ) &MutableMesh2_2::GetNumNodes,
            " "  )
        .def(
            "GetNumElements",
            (unsigned int(MutableMesh2_2::*)() const ) &MutableMesh2_2::GetNumElements,
            " "  )
        .def(
            "GetNumBoundaryElements",
            (unsigned int(MutableMesh2_2::*)() const ) &MutableMesh2_2::GetNumBoundaryElements,
            " "  )
        .def(
            "AddNode",
            (unsigned int(MutableMesh2_2::*)(::Node<2> *)) &MutableMesh2_2::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "AddElement",
            (unsigned int(MutableMesh2_2::*)(::Element<2, 2> *)) &MutableMesh2_2::AddElement,
            " " , py::arg("pNewElement") )
        .def(
            "SetNode",
            (void(MutableMesh2_2::*)(unsigned int, ::ChastePoint<2>, bool)) &MutableMesh2_2::SetNode,
            " " , py::arg("index"), py::arg("point"), py::arg("concreteMove") = true )
        .def(
            "MoveMergeNode",
            (void(MutableMesh2_2::*)(unsigned int, unsigned int, bool)) &MutableMesh2_2::MoveMergeNode,
            " " , py::arg("index"), py::arg("targetIndex"), py::arg("concreteMove") = true )
        .def(
            "DeleteNode",
            (void(MutableMesh2_2::*)(unsigned int)) &MutableMesh2_2::DeleteNode,
            " " , py::arg("index") )
        .def(
            "DeleteElement",
            (void(MutableMesh2_2::*)(unsigned int)) &MutableMesh2_2::DeleteElement,
            " " , py::arg("index") )
        .def(
            "DeleteNodePriorToReMesh",
            (void(MutableMesh2_2::*)(unsigned int)) &MutableMesh2_2::DeleteNodePriorToReMesh,
            " " , py::arg("index") )
        .def(
            "RefineElement",
            (unsigned int(MutableMesh2_2::*)(::Element<2, 2> *, ::ChastePoint<2>)) &MutableMesh2_2::RefineElement,
            " " , py::arg("pElement"), py::arg("point") )
        .def(
            "DeleteBoundaryNodeAt",
            (void(MutableMesh2_2::*)(unsigned int)) &MutableMesh2_2::DeleteBoundaryNodeAt,
            " " , py::arg("index") )
        .def(
            "ReIndex",
            (void(MutableMesh2_2::*)(::NodeMap &)) &MutableMesh2_2::ReIndex,
            " " , py::arg("map") )
        .def(
            "ReMesh",
            (void(MutableMesh2_2::*)(::NodeMap &)) &MutableMesh2_2::ReMesh,
            " " , py::arg("map") )
        .def(
            "ReMesh",
            (void(MutableMesh2_2::*)()) &MutableMesh2_2::ReMesh,
            " "  )
        .def(
            "SplitEdge",
            (::boost::numeric::ublas::c_vector<unsigned int, 3>(MutableMesh2_2::*)(::Node<2> *, ::Node<2> *)) &MutableMesh2_2::SplitEdge,
            " " , py::arg("pNodeA"), py::arg("pNodeB") )
        .def(
            "CheckIsVoronoi",
            (bool(MutableMesh2_2::*)(double)) &MutableMesh2_2::CheckIsVoronoi,
            " " , py::arg("maxPenetration") = 0. )
    ;
}
