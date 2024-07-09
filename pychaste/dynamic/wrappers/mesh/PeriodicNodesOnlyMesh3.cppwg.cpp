#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PeriodicNodesOnlyMesh.hpp"

#include "PeriodicNodesOnlyMesh3.cppwg.hpp"

namespace py = pybind11;
typedef PeriodicNodesOnlyMesh<3 > PeriodicNodesOnlyMesh3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef unsigned int unsignedint;

class PeriodicNodesOnlyMesh3_Overrides : public PeriodicNodesOnlyMesh3{
    public:
    using PeriodicNodesOnlyMesh3::PeriodicNodesOnlyMesh;
    ::boost::numeric::ublas::c_vector<double, 3> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 3> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 3> const & rLocation2) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            PeriodicNodesOnlyMesh3,
            GetVectorFromAtoB,
                    rLocation1,
        rLocation2);
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            PeriodicNodesOnlyMesh3,
            GetWidth,
                    rDimension);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<3> point, bool concreteMove) override {
        PYBIND11_OVERRIDE(
            void,
            PeriodicNodesOnlyMesh3,
            SetNode,
                    nodeIndex,
        point,
        concreteMove);
    }
    unsigned int AddNode(::Node<3> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            PeriodicNodesOnlyMesh3,
            AddNode,
                    pNewNode);
    }
    void RefreshMesh() override {
        PYBIND11_OVERRIDE(
            void,
            PeriodicNodesOnlyMesh3,
            RefreshMesh,
            );
    }

};
void register_PeriodicNodesOnlyMesh3_class(py::module &m){
py::class_<PeriodicNodesOnlyMesh3 , PeriodicNodesOnlyMesh3_Overrides , boost::shared_ptr<PeriodicNodesOnlyMesh3 >  , NodesOnlyMesh<3>  >(m, "PeriodicNodesOnlyMesh3")
        .def(py::init<::boost::numeric::ublas::c_vector<double, 3> >(), py::arg("width"))
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 3>(PeriodicNodesOnlyMesh3::*)(::boost::numeric::ublas::c_vector<double, 3> const &, ::boost::numeric::ublas::c_vector<double, 3> const &)) &PeriodicNodesOnlyMesh3::GetVectorFromAtoB,
            " " , py::arg("rLocation1"), py::arg("rLocation2") )
        .def(
            "GetWidth",
            (double(PeriodicNodesOnlyMesh3::*)(unsigned int const &) const ) &PeriodicNodesOnlyMesh3::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "GetPeriodicWidths",
            (::boost::numeric::ublas::c_vector<double, 3>(PeriodicNodesOnlyMesh3::*)() const ) &PeriodicNodesOnlyMesh3::GetPeriodicWidths,
            " "  )
        .def(
            "SetNode",
            (void(PeriodicNodesOnlyMesh3::*)(unsigned int, ::ChastePoint<3>, bool)) &PeriodicNodesOnlyMesh3::SetNode,
            " " , py::arg("nodeIndex"), py::arg("point"), py::arg("concreteMove") = false )
        .def(
            "AddNode",
            (unsigned int(PeriodicNodesOnlyMesh3::*)(::Node<3> *)) &PeriodicNodesOnlyMesh3::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "RefreshMesh",
            (void(PeriodicNodesOnlyMesh3::*)()) &PeriodicNodesOnlyMesh3::RefreshMesh,
            " "  )
    ;
}
