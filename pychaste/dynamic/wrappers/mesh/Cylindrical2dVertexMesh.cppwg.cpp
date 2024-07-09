#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Cylindrical2dVertexMesh.hpp"

#include "Cylindrical2dVertexMesh.cppwg.hpp"

namespace py = pybind11;
typedef Cylindrical2dVertexMesh Cylindrical2dVertexMesh;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::VertexMesh<2, 2> * _VertexMesh_lt_2_2_gt_Ptr;

class Cylindrical2dVertexMesh_Overrides : public Cylindrical2dVertexMesh{
    public:
    using Cylindrical2dVertexMesh::Cylindrical2dVertexMesh;
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 2> const & rLocation2) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            Cylindrical2dVertexMesh,
            GetVectorFromAtoB,
                    rLocation1,
        rLocation2);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> point) override {
        PYBIND11_OVERRIDE(
            void,
            Cylindrical2dVertexMesh,
            SetNode,
                    nodeIndex,
        point);
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            Cylindrical2dVertexMesh,
            GetWidth,
                    rDimension);
    }
    void Scale(double const xScale, double const yScale, double const zScale) override {
        PYBIND11_OVERRIDE(
            void,
            Cylindrical2dVertexMesh,
            Scale,
                    xScale,
        yScale,
        zScale);
    }
    ::VertexMesh<2, 2> * GetMeshForVtk() override {
        PYBIND11_OVERRIDE(
            _VertexMesh_lt_2_2_gt_Ptr,
            Cylindrical2dVertexMesh,
            GetMeshForVtk,
            );
    }

};
void register_Cylindrical2dVertexMesh_class(py::module &m){
py::class_<Cylindrical2dVertexMesh , Cylindrical2dVertexMesh_Overrides , boost::shared_ptr<Cylindrical2dVertexMesh >  , MutableVertexMesh<2, 2>  >(m, "Cylindrical2dVertexMesh")
        .def(py::init<double, ::std::vector<Node<2> *>, ::std::vector<VertexElement<2, 2> *>, double, double >(), py::arg("width"), py::arg("nodes"), py::arg("vertexElements"), py::arg("cellRearrangementThreshold") = 0.01, py::arg("t2Threshold") = 0.001)
        .def(py::init<::Cylindrical2dMesh &, bool >(), py::arg("rMesh"), py::arg("isBounded") = false)
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(Cylindrical2dVertexMesh::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &Cylindrical2dVertexMesh::GetVectorFromAtoB,
            " " , py::arg("rLocation1"), py::arg("rLocation2") )
        .def(
            "SetNode",
            (void(Cylindrical2dVertexMesh::*)(unsigned int, ::ChastePoint<2>)) &Cylindrical2dVertexMesh::SetNode,
            " " , py::arg("nodeIndex"), py::arg("point") )
        .def(
            "GetWidth",
            (double(Cylindrical2dVertexMesh::*)(unsigned int const &) const ) &Cylindrical2dVertexMesh::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "AddNode",
            (unsigned int(Cylindrical2dVertexMesh::*)(::Node<2> *)) &Cylindrical2dVertexMesh::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "CheckNodeLocation",
            (void(Cylindrical2dVertexMesh::*)(::Node<2> *)) &Cylindrical2dVertexMesh::CheckNodeLocation,
            " " , py::arg("pNode") )
        .def(
            "Scale",
            (void(Cylindrical2dVertexMesh::*)(double const, double const, double const)) &Cylindrical2dVertexMesh::Scale,
            " " , py::arg("xScale") = 1., py::arg("yScale") = 1., py::arg("zScale") = 1. )
        .def(
            "GetMeshForVtk",
            (::VertexMesh<2, 2> *(Cylindrical2dVertexMesh::*)()) &Cylindrical2dVertexMesh::GetMeshForVtk,
            " "  , py::return_value_policy::reference)
    ;
}
