#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VertexMesh.hpp"

#include "VertexMesh3_3.cppwg.hpp"

namespace py = pybind11;
typedef VertexMesh<3,3 > VertexMesh3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::VertexMesh<3, 3> * _VertexMesh_lt_3_3_gt_Ptr;
typedef unsigned int unsignedint;

class VertexMesh3_3_Overrides : public VertexMesh3_3{
    public:
    using VertexMesh3_3::VertexMesh;
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexMesh3_3,
            GetNumNodes,
            );
    }
    unsigned int GetNumElements() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexMesh3_3,
            GetNumElements,
            );
    }
    unsigned int GetNumFaces() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexMesh3_3,
            GetNumFaces,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 3> GetCentroidOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            VertexMesh3_3,
            GetCentroidOfElement,
                    index);
    }
    void Clear() override {
        PYBIND11_OVERRIDE(
            void,
            VertexMesh3_3,
            Clear,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 3> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 3> const & rLocationA, ::boost::numeric::ublas::c_vector<double, 3> const & rLocationB) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            VertexMesh3_3,
            GetVectorFromAtoB,
                    rLocationA,
        rLocationB);
    }
    double GetVolumeOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            double,
            VertexMesh3_3,
            GetVolumeOfElement,
                    index);
    }
    double GetSurfaceAreaOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            double,
            VertexMesh3_3,
            GetSurfaceAreaOfElement,
                    index);
    }
    ::boost::numeric::ublas::c_vector<double, 3> CalculateMomentsOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            VertexMesh3_3,
            CalculateMomentsOfElement,
                    index);
    }
    double CalculateAreaOfFace(::VertexElement<2, 3> * pFace) override {
        PYBIND11_OVERRIDE(
            double,
            VertexMesh3_3,
            CalculateAreaOfFace,
                    pFace);
    }
    ::VertexMesh<3, 3> * GetMeshForVtk() override {
        PYBIND11_OVERRIDE(
            _VertexMesh_lt_3_3_gt_Ptr,
            VertexMesh3_3,
            GetMeshForVtk,
            );
    }
    unsigned int SolveNodeMapping(unsigned int index) const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexMesh3_3,
            SolveNodeMapping,
                    index);
    }

};
void register_VertexMesh3_3_class(py::module &m){
py::class_<VertexMesh3_3 , VertexMesh3_3_Overrides , boost::shared_ptr<VertexMesh3_3 >  , AbstractMesh<3, 3>  >(m, "VertexMesh3_3")
        .def(py::init<::std::vector<Node<3> *>, ::std::vector<VertexElement<3, 3> *> >(), py::arg("nodes"), py::arg("vertexElements"))
        .def(py::init<::std::vector<Node<3> *>, ::std::vector<VertexElement<2, 3> *>, ::std::vector<VertexElement<3, 3> *> >(), py::arg("nodes"), py::arg("faces"), py::arg("vertexElements"))
        .def(py::init< >())
        .def(
            "GetElementIteratorBegin",
            (::VertexMesh<3, 3>::VertexElementIterator(VertexMesh3_3::*)(bool)) &VertexMesh3_3::GetElementIteratorBegin,
            " " , py::arg("skipDeletedElements") = true )
        .def(
            "GetElementIteratorEnd",
            (::VertexMesh<3, 3>::VertexElementIterator(VertexMesh3_3::*)()) &VertexMesh3_3::GetElementIteratorEnd,
            " "  )
        .def(
            "GetNumEdges",
            (unsigned int(VertexMesh3_3::*)() const ) &VertexMesh3_3::GetNumEdges,
            " "  )
        .def(
            "GetEdge",
            (::Edge<3> *(VertexMesh3_3::*)(unsigned int) const ) &VertexMesh3_3::GetEdge,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNumNodes",
            (unsigned int(VertexMesh3_3::*)() const ) &VertexMesh3_3::GetNumNodes,
            " "  )
        .def(
            "GetNumElements",
            (unsigned int(VertexMesh3_3::*)() const ) &VertexMesh3_3::GetNumElements,
            " "  )
        .def(
            "GetNumAllElements",
            (unsigned int(VertexMesh3_3::*)() const ) &VertexMesh3_3::GetNumAllElements,
            " "  )
        .def(
            "GetNumFaces",
            (unsigned int(VertexMesh3_3::*)() const ) &VertexMesh3_3::GetNumFaces,
            " "  )
        .def(
            "GetElement",
            (::VertexElement<3, 3> *(VertexMesh3_3::*)(unsigned int) const ) &VertexMesh3_3::GetElement,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetCentroidOfElement",
            (::boost::numeric::ublas::c_vector<double, 3>(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::GetCentroidOfElement,
            " " , py::arg("index") )
        .def(
            "ConstructFromMeshReader",
            (void(VertexMesh3_3::*)(::AbstractMeshReader<3, 3> &)) &VertexMesh3_3::ConstructFromMeshReader,
            " " , py::arg("rMeshReader") )
        .def(
            "Clear",
            (void(VertexMesh3_3::*)()) &VertexMesh3_3::Clear,
            " "  )
        .def(
            "GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex",
            (unsigned int(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex,
            " " , py::arg("elementIndex") )
        .def(
            "GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex",
            (unsigned int(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex,
            " " , py::arg("nodeIndex") )
        .def(
            "GetRosetteRankOfElement",
            (unsigned int(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::GetRosetteRankOfElement,
            " " , py::arg("index") )
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 3>(VertexMesh3_3::*)(::boost::numeric::ublas::c_vector<double, 3> const &, ::boost::numeric::ublas::c_vector<double, 3> const &)) &VertexMesh3_3::GetVectorFromAtoB,
            " " , py::arg("rLocationA"), py::arg("rLocationB") )
        .def(
            "GetVolumeOfElement",
            (double(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::GetVolumeOfElement,
            " " , py::arg("index") )
        .def(
            "GetSurfaceAreaOfElement",
            (double(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::GetSurfaceAreaOfElement,
            " " , py::arg("index") )
        .def(
            "GetAreaGradientOfElementAtNode",
            (::boost::numeric::ublas::c_vector<double, 3>(VertexMesh3_3::*)(::VertexElement<3, 3> *, unsigned int)) &VertexMesh3_3::GetAreaGradientOfElementAtNode,
            " " , py::arg("pElement"), py::arg("localIndex") )
        .def(
            "GetPreviousEdgeGradientOfElementAtNode",
            (::boost::numeric::ublas::c_vector<double, 3>(VertexMesh3_3::*)(::VertexElement<3, 3> *, unsigned int)) &VertexMesh3_3::GetPreviousEdgeGradientOfElementAtNode,
            " " , py::arg("pElement"), py::arg("localIndex") )
        .def(
            "GetNextEdgeGradientOfElementAtNode",
            (::boost::numeric::ublas::c_vector<double, 3>(VertexMesh3_3::*)(::VertexElement<3, 3> *, unsigned int)) &VertexMesh3_3::GetNextEdgeGradientOfElementAtNode,
            " " , py::arg("pElement"), py::arg("localIndex") )
        .def(
            "GetPerimeterGradientOfElementAtNode",
            (::boost::numeric::ublas::c_vector<double, 3>(VertexMesh3_3::*)(::VertexElement<3, 3> *, unsigned int)) &VertexMesh3_3::GetPerimeterGradientOfElementAtNode,
            " " , py::arg("pElement"), py::arg("localIndex") )
        .def(
            "CalculateMomentsOfElement",
            (::boost::numeric::ublas::c_vector<double, 3>(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::CalculateMomentsOfElement,
            " " , py::arg("index") )
        .def(
            "GetEdgeLength",
            (double(VertexMesh3_3::*)(unsigned int, unsigned int)) &VertexMesh3_3::GetEdgeLength,
            " " , py::arg("elementIndex1"), py::arg("elementIndex2") )
        .def(
            "GetElongationShapeFactorOfElement",
            (double(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::GetElongationShapeFactorOfElement,
            " " , py::arg("elementIndex") )
        .def(
            "CalculateUnitNormalToFaceWithArea",
            (double(VertexMesh3_3::*)(::VertexElement<2, 3> *, ::boost::numeric::ublas::c_vector<double, 3> &)) &VertexMesh3_3::CalculateUnitNormalToFaceWithArea,
            " " , py::arg("pFace"), py::arg("rNormal") )
        .def(
            "CalculateAreaOfFace",
            (double(VertexMesh3_3::*)(::VertexElement<2, 3> *)) &VertexMesh3_3::CalculateAreaOfFace,
            " " , py::arg("pFace") )
        .def(
            "GetShortAxisOfElement",
            (::boost::numeric::ublas::c_vector<double, 3>(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::GetShortAxisOfElement,
            " " , py::arg("index") )
        .def(
            "GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::GetNeighbouringNodeIndices,
            " " , py::arg("nodeIndex") )
        .def(
            "GetNeighbouringNodeNotAlsoInElement",
            (::std::set<unsigned int>(VertexMesh3_3::*)(unsigned int, unsigned int)) &VertexMesh3_3::GetNeighbouringNodeNotAlsoInElement,
            " " , py::arg("nodeIndex"), py::arg("elemIndex") )
        .def(
            "GetNeighbouringElementIndices",
            (::std::set<unsigned int>(VertexMesh3_3::*)(unsigned int)) &VertexMesh3_3::GetNeighbouringElementIndices,
            " " , py::arg("elementIndex") )
        .def(
            "GetMeshForVtk",
            (::VertexMesh<3, 3> *(VertexMesh3_3::*)()) &VertexMesh3_3::GetMeshForVtk,
            " "  , py::return_value_policy::reference)
    ;
}
