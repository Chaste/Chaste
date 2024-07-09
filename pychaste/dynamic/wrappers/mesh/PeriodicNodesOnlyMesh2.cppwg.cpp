#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PeriodicNodesOnlyMesh.hpp"

#include "PeriodicNodesOnlyMesh2.cppwg.hpp"

namespace py = pybind11;
typedef PeriodicNodesOnlyMesh<2 > PeriodicNodesOnlyMesh2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef unsigned int unsignedint;

class PeriodicNodesOnlyMesh2_Overrides : public PeriodicNodesOnlyMesh2{
    public:
    using PeriodicNodesOnlyMesh2::PeriodicNodesOnlyMesh;
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 2> const & rLocation2) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            PeriodicNodesOnlyMesh2,
            GetVectorFromAtoB,
                    rLocation1,
        rLocation2);
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            PeriodicNodesOnlyMesh2,
            GetWidth,
                    rDimension);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> point, bool concreteMove) override {
        PYBIND11_OVERRIDE(
            void,
            PeriodicNodesOnlyMesh2,
            SetNode,
                    nodeIndex,
        point,
        concreteMove);
    }
    unsigned int AddNode(::Node<2> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            PeriodicNodesOnlyMesh2,
            AddNode,
                    pNewNode);
    }
    void RefreshMesh() override {
        PYBIND11_OVERRIDE(
            void,
            PeriodicNodesOnlyMesh2,
            RefreshMesh,
            );
    }

};
void register_PeriodicNodesOnlyMesh2_class(py::module &m){
py::class_<PeriodicNodesOnlyMesh2 , PeriodicNodesOnlyMesh2_Overrides , boost::shared_ptr<PeriodicNodesOnlyMesh2 >  , NodesOnlyMesh<2>  >(m, "PeriodicNodesOnlyMesh2")
        .def(py::init<::boost::numeric::ublas::c_vector<double, 2> >(), py::arg("width"))
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(PeriodicNodesOnlyMesh2::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &PeriodicNodesOnlyMesh2::GetVectorFromAtoB,
            " " , py::arg("rLocation1"), py::arg("rLocation2") )
        .def(
            "GetWidth",
            (double(PeriodicNodesOnlyMesh2::*)(unsigned int const &) const ) &PeriodicNodesOnlyMesh2::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "GetPeriodicWidths",
            (::boost::numeric::ublas::c_vector<double, 2>(PeriodicNodesOnlyMesh2::*)() const ) &PeriodicNodesOnlyMesh2::GetPeriodicWidths,
            " "  )
        .def(
            "SetNode",
            (void(PeriodicNodesOnlyMesh2::*)(unsigned int, ::ChastePoint<2>, bool)) &PeriodicNodesOnlyMesh2::SetNode,
            " " , py::arg("nodeIndex"), py::arg("point"), py::arg("concreteMove") = false )
        .def(
            "AddNode",
            (unsigned int(PeriodicNodesOnlyMesh2::*)(::Node<2> *)) &PeriodicNodesOnlyMesh2::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "RefreshMesh",
            (void(PeriodicNodesOnlyMesh2::*)()) &PeriodicNodesOnlyMesh2::RefreshMesh,
            " "  )
    ;
}
