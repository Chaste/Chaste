#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "MutableMesh.hpp"

#include "MutableMesh3_3.cppwg.hpp"

namespace py = pybind11;
typedef MutableMesh<3,3 > MutableMesh3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class MutableMesh3_3_Overrides : public MutableMesh3_3{
    public:
    using MutableMesh3_3::MutableMesh;
    void Clear() override {
        PYBIND11_OVERRIDE(
            void,
            MutableMesh3_3,
            Clear,
            );
    }
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableMesh3_3,
            GetNumNodes,
            );
    }
    unsigned int GetNumElements() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableMesh3_3,
            GetNumElements,
            );
    }
    unsigned int GetNumBoundaryElements() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableMesh3_3,
            GetNumBoundaryElements,
            );
    }
    unsigned int AddNode(::Node<3> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableMesh3_3,
            AddNode,
                    pNewNode);
    }
    void SetNode(unsigned int index, ::ChastePoint<3> point, bool concreteMove) override {
        PYBIND11_OVERRIDE(
            void,
            MutableMesh3_3,
            SetNode,
                    index,
        point,
        concreteMove);
    }
    void DeleteNode(unsigned int index) override {
        PYBIND11_OVERRIDE(
            void,
            MutableMesh3_3,
            DeleteNode,
                    index);
    }
    void DeleteElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            void,
            MutableMesh3_3,
            DeleteElement,
                    index);
    }
    void ReMesh(::NodeMap & map) override {
        PYBIND11_OVERRIDE(
            void,
            MutableMesh3_3,
            ReMesh,
                    map);
    }

};
void register_MutableMesh3_3_class(py::module &m){
py::class_<MutableMesh3_3 , MutableMesh3_3_Overrides , boost::shared_ptr<MutableMesh3_3 >  , TetrahedralMesh<3, 3>  >(m, "MutableMesh3_3")
        .def(py::init< >())
        .def(py::init<::std::vector<Node<3> *> >(), py::arg("nodes"))
        .def(
            "Clear",
            (void(MutableMesh3_3::*)()) &MutableMesh3_3::Clear,
            " "  )
        .def(
            "GetNumNodes",
            (unsigned int(MutableMesh3_3::*)() const ) &MutableMesh3_3::GetNumNodes,
            " "  )
        .def(
            "GetNumElements",
            (unsigned int(MutableMesh3_3::*)() const ) &MutableMesh3_3::GetNumElements,
            " "  )
        .def(
            "GetNumBoundaryElements",
            (unsigned int(MutableMesh3_3::*)() const ) &MutableMesh3_3::GetNumBoundaryElements,
            " "  )
        .def(
            "AddNode",
            (unsigned int(MutableMesh3_3::*)(::Node<3> *)) &MutableMesh3_3::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "AddElement",
            (unsigned int(MutableMesh3_3::*)(::Element<3, 3> *)) &MutableMesh3_3::AddElement,
            " " , py::arg("pNewElement") )
        .def(
            "SetNode",
            (void(MutableMesh3_3::*)(unsigned int, ::ChastePoint<3>, bool)) &MutableMesh3_3::SetNode,
            " " , py::arg("index"), py::arg("point"), py::arg("concreteMove") = true )
        .def(
            "MoveMergeNode",
            (void(MutableMesh3_3::*)(unsigned int, unsigned int, bool)) &MutableMesh3_3::MoveMergeNode,
            " " , py::arg("index"), py::arg("targetIndex"), py::arg("concreteMove") = true )
        .def(
            "DeleteNode",
            (void(MutableMesh3_3::*)(unsigned int)) &MutableMesh3_3::DeleteNode,
            " " , py::arg("index") )
        .def(
            "DeleteElement",
            (void(MutableMesh3_3::*)(unsigned int)) &MutableMesh3_3::DeleteElement,
            " " , py::arg("index") )
        .def(
            "DeleteNodePriorToReMesh",
            (void(MutableMesh3_3::*)(unsigned int)) &MutableMesh3_3::DeleteNodePriorToReMesh,
            " " , py::arg("index") )
        .def(
            "RefineElement",
            (unsigned int(MutableMesh3_3::*)(::Element<3, 3> *, ::ChastePoint<3>)) &MutableMesh3_3::RefineElement,
            " " , py::arg("pElement"), py::arg("point") )
        .def(
            "DeleteBoundaryNodeAt",
            (void(MutableMesh3_3::*)(unsigned int)) &MutableMesh3_3::DeleteBoundaryNodeAt,
            " " , py::arg("index") )
        .def(
            "ReIndex",
            (void(MutableMesh3_3::*)(::NodeMap &)) &MutableMesh3_3::ReIndex,
            " " , py::arg("map") )
        .def(
            "ReMesh",
            (void(MutableMesh3_3::*)(::NodeMap &)) &MutableMesh3_3::ReMesh,
            " " , py::arg("map") )
        .def(
            "ReMesh",
            (void(MutableMesh3_3::*)()) &MutableMesh3_3::ReMesh,
            " "  )
        .def(
            "SplitEdge",
            (::boost::numeric::ublas::c_vector<unsigned int, 3>(MutableMesh3_3::*)(::Node<3> *, ::Node<3> *)) &MutableMesh3_3::SplitEdge,
            " " , py::arg("pNodeA"), py::arg("pNodeB") )
        .def(
            "CheckIsVoronoi",
            (bool(MutableMesh3_3::*)(double)) &MutableMesh3_3::CheckIsVoronoi,
            " " , py::arg("maxPenetration") = 0. )
    ;
}
