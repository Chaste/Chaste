#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"

#include "Cylindrical2dNodesOnlyMesh.cppwg.hpp"

namespace py = pybind11;
typedef Cylindrical2dNodesOnlyMesh Cylindrical2dNodesOnlyMesh;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef unsigned int unsignedint;

class Cylindrical2dNodesOnlyMesh_Overrides : public Cylindrical2dNodesOnlyMesh{
    public:
    using Cylindrical2dNodesOnlyMesh::Cylindrical2dNodesOnlyMesh;
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 2> const & rLocation2) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            Cylindrical2dNodesOnlyMesh,
            GetVectorFromAtoB,
                    rLocation1,
        rLocation2);
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            Cylindrical2dNodesOnlyMesh,
            GetWidth,
                    rDimension);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> point, bool concreteMove) override {
        PYBIND11_OVERRIDE(
            void,
            Cylindrical2dNodesOnlyMesh,
            SetNode,
                    nodeIndex,
        point,
        concreteMove);
    }
    unsigned int AddNode(::Node<2> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            Cylindrical2dNodesOnlyMesh,
            AddNode,
                    pNewNode);
    }
    void RefreshMesh() override {
        PYBIND11_OVERRIDE(
            void,
            Cylindrical2dNodesOnlyMesh,
            RefreshMesh,
            );
    }

};
void register_Cylindrical2dNodesOnlyMesh_class(py::module &m){
py::class_<Cylindrical2dNodesOnlyMesh , Cylindrical2dNodesOnlyMesh_Overrides , boost::shared_ptr<Cylindrical2dNodesOnlyMesh >  , NodesOnlyMesh<2>  >(m, "Cylindrical2dNodesOnlyMesh")
        .def(py::init<double >(), py::arg("width"))
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(Cylindrical2dNodesOnlyMesh::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &Cylindrical2dNodesOnlyMesh::GetVectorFromAtoB,
            " " , py::arg("rLocation1"), py::arg("rLocation2") )
        .def(
            "GetWidth",
            (double(Cylindrical2dNodesOnlyMesh::*)(unsigned int const &) const ) &Cylindrical2dNodesOnlyMesh::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "SetNode",
            (void(Cylindrical2dNodesOnlyMesh::*)(unsigned int, ::ChastePoint<2>, bool)) &Cylindrical2dNodesOnlyMesh::SetNode,
            " " , py::arg("nodeIndex"), py::arg("point"), py::arg("concreteMove") = false )
        .def(
            "AddNode",
            (unsigned int(Cylindrical2dNodesOnlyMesh::*)(::Node<2> *)) &Cylindrical2dNodesOnlyMesh::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "RefreshMesh",
            (void(Cylindrical2dNodesOnlyMesh::*)()) &Cylindrical2dNodesOnlyMesh::RefreshMesh,
            " "  )
    ;
}
