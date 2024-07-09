#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Toroidal2dVertexMesh.hpp"

#include "Toroidal2dVertexMesh.cppwg.hpp"

namespace py = pybind11;
typedef Toroidal2dVertexMesh Toroidal2dVertexMesh;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::VertexMesh<2, 2> * _VertexMesh_lt_2_2_gt_Ptr;

class Toroidal2dVertexMesh_Overrides : public Toroidal2dVertexMesh{
    public:
    using Toroidal2dVertexMesh::Toroidal2dVertexMesh;
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 2> const & rLocation2) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            Toroidal2dVertexMesh,
            GetVectorFromAtoB,
                    rLocation1,
        rLocation2);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> point) override {
        PYBIND11_OVERRIDE(
            void,
            Toroidal2dVertexMesh,
            SetNode,
                    nodeIndex,
        point);
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            Toroidal2dVertexMesh,
            GetWidth,
                    rDimension);
    }
    ::VertexMesh<2, 2> * GetMeshForVtk() override {
        PYBIND11_OVERRIDE(
            _VertexMesh_lt_2_2_gt_Ptr,
            Toroidal2dVertexMesh,
            GetMeshForVtk,
            );
    }

};
void register_Toroidal2dVertexMesh_class(py::module &m){
py::class_<Toroidal2dVertexMesh , Toroidal2dVertexMesh_Overrides , boost::shared_ptr<Toroidal2dVertexMesh >  , MutableVertexMesh<2, 2>  >(m, "Toroidal2dVertexMesh")
        .def(py::init<double, double, ::std::vector<Node<2> *>, ::std::vector<VertexElement<2, 2> *>, double, double >(), py::arg("width"), py::arg("height"), py::arg("nodes"), py::arg("vertexElements"), py::arg("cellRearrangementThreshold") = 0.01, py::arg("t2Threshold") = 0.001)
        .def(py::init< >())
        .def(py::init<::Toroidal2dMesh &, bool >(), py::arg("rMesh"), py::arg("isBounded") = false)
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(Toroidal2dVertexMesh::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &Toroidal2dVertexMesh::GetVectorFromAtoB,
            " " , py::arg("rLocation1"), py::arg("rLocation2") )
        .def(
            "SetNode",
            (void(Toroidal2dVertexMesh::*)(unsigned int, ::ChastePoint<2>)) &Toroidal2dVertexMesh::SetNode,
            " " , py::arg("nodeIndex"), py::arg("point") )
        .def(
            "GetWidth",
            (double(Toroidal2dVertexMesh::*)(unsigned int const &) const ) &Toroidal2dVertexMesh::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "SetHeight",
            (void(Toroidal2dVertexMesh::*)(double)) &Toroidal2dVertexMesh::SetHeight,
            " " , py::arg("height") )
        .def(
            "SetWidth",
            (void(Toroidal2dVertexMesh::*)(double)) &Toroidal2dVertexMesh::SetWidth,
            " " , py::arg("width") )
        .def(
            "AddNode",
            (unsigned int(Toroidal2dVertexMesh::*)(::Node<2> *)) &Toroidal2dVertexMesh::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "CheckNodeLocation",
            (void(Toroidal2dVertexMesh::*)(::Node<2> *)) &Toroidal2dVertexMesh::CheckNodeLocation,
            " " , py::arg("pNode") )
        .def(
            "GetMeshForVtk",
            (::VertexMesh<2, 2> *(Toroidal2dVertexMesh::*)()) &Toroidal2dVertexMesh::GetMeshForVtk,
            " "  , py::return_value_policy::reference)
        .def(
            "ConstructFromMeshReader",
            (void(Toroidal2dVertexMesh::*)(::AbstractMeshReader<2, 2> &, double, double)) &Toroidal2dVertexMesh::ConstructFromMeshReader,
            " " , py::arg("rMeshReader"), py::arg("width"), py::arg("height") )
    ;
}
