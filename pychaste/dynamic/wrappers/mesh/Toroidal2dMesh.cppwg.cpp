#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Toroidal2dMesh.hpp"

#include "Toroidal2dMesh.cppwg.hpp"

namespace py = pybind11;
typedef Toroidal2dMesh Toroidal2dMesh;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef unsigned int unsignedint;

class Toroidal2dMesh_Overrides : public Toroidal2dMesh{
    public:
    using Toroidal2dMesh::Toroidal2dMesh;
    void ReMesh(::NodeMap & rMap) override {
        PYBIND11_OVERRIDE(
            void,
            Toroidal2dMesh,
            ReMesh,
                    rMap);
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 2> const & rLocation2) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            Toroidal2dMesh,
            GetVectorFromAtoB,
                    rLocation1,
        rLocation2);
    }
    void SetNode(unsigned int index, ::ChastePoint<2> point, bool concreteMove) override {
        PYBIND11_OVERRIDE(
            void,
            Toroidal2dMesh,
            SetNode,
                    index,
        point,
        concreteMove);
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            Toroidal2dMesh,
            GetWidth,
                    rDimension);
    }
    unsigned int AddNode(::Node<2> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            Toroidal2dMesh,
            AddNode,
                    pNewNode);
    }
    void RefreshMesh() override {
        PYBIND11_OVERRIDE(
            void,
            Toroidal2dMesh,
            RefreshMesh,
            );
    }

};
void register_Toroidal2dMesh_class(py::module &m){
py::class_<Toroidal2dMesh , Toroidal2dMesh_Overrides , boost::shared_ptr<Toroidal2dMesh >  , MutableMesh<2, 2>  >(m, "Toroidal2dMesh")
        .def(py::init<double, double >(), py::arg("width"), py::arg("depth"))
        .def(py::init<double, double, ::std::vector<Node<2> *> >(), py::arg("width"), py::arg("depth"), py::arg("nodes"))
        .def(
            "ReMesh",
            (void(Toroidal2dMesh::*)(::NodeMap &)) &Toroidal2dMesh::ReMesh,
            " " , py::arg("rMap") )
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(Toroidal2dMesh::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &Toroidal2dMesh::GetVectorFromAtoB,
            " " , py::arg("rLocation1"), py::arg("rLocation2") )
        .def(
            "SetNode",
            (void(Toroidal2dMesh::*)(unsigned int, ::ChastePoint<2>, bool)) &Toroidal2dMesh::SetNode,
            " " , py::arg("index"), py::arg("point"), py::arg("concreteMove") )
        .def(
            "GetWidth",
            (double(Toroidal2dMesh::*)(unsigned int const &) const ) &Toroidal2dMesh::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "AddNode",
            (unsigned int(Toroidal2dMesh::*)(::Node<2> *)) &Toroidal2dMesh::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "RefreshMesh",
            (void(Toroidal2dMesh::*)()) &Toroidal2dMesh::RefreshMesh,
            " "  )
    ;
}
