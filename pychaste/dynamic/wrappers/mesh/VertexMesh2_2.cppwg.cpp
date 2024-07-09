#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VertexMesh.hpp"

#include "VertexMesh2_2.cppwg.hpp"

namespace py = pybind11;
typedef VertexMesh<2,2 > VertexMesh2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::VertexMesh<2, 2> * _VertexMesh_lt_2_2_gt_Ptr;
typedef unsigned int unsignedint;

class VertexMesh2_2_Overrides : public VertexMesh2_2{
    public:
    using VertexMesh2_2::VertexMesh;
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexMesh2_2,
            GetNumNodes,
            );
    }
    unsigned int GetNumElements() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexMesh2_2,
            GetNumElements,
            );
    }
    unsigned int GetNumFaces() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexMesh2_2,
            GetNumFaces,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetCentroidOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            VertexMesh2_2,
            GetCentroidOfElement,
                    index);
    }
    void Clear() override {
        PYBIND11_OVERRIDE(
            void,
            VertexMesh2_2,
            Clear,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocationA, ::boost::numeric::ublas::c_vector<double, 2> const & rLocationB) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            VertexMesh2_2,
            GetVectorFromAtoB,
                    rLocationA,
        rLocationB);
    }
    double GetVolumeOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            double,
            VertexMesh2_2,
            GetVolumeOfElement,
                    index);
    }
    double GetSurfaceAreaOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            double,
            VertexMesh2_2,
            GetSurfaceAreaOfElement,
                    index);
    }
    ::boost::numeric::ublas::c_vector<double, 3> CalculateMomentsOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            VertexMesh2_2,
            CalculateMomentsOfElement,
                    index);
    }
    double CalculateAreaOfFace(::VertexElement<1, 2> * pFace) override {
        PYBIND11_OVERRIDE(
            double,
            VertexMesh2_2,
            CalculateAreaOfFace,
                    pFace);
    }
    ::VertexMesh<2, 2> * GetMeshForVtk() override {
        PYBIND11_OVERRIDE(
            _VertexMesh_lt_2_2_gt_Ptr,
            VertexMesh2_2,
            GetMeshForVtk,
            );
    }
    unsigned int SolveNodeMapping(unsigned int index) const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexMesh2_2,
            SolveNodeMapping,
                    index);
    }

};
void register_VertexMesh2_2_class(py::module &m){
py::class_<VertexMesh2_2 , VertexMesh2_2_Overrides , boost::shared_ptr<VertexMesh2_2 >  , AbstractMesh<2, 2>  >(m, "VertexMesh2_2")
        .def(py::init<::std::vector<Node<2> *>, ::std::vector<VertexElement<2, 2> *> >(), py::arg("nodes"), py::arg("vertexElements"))
        .def(py::init<::std::vector<Node<2> *>, ::std::vector<VertexElement<1, 2> *>, ::std::vector<VertexElement<2, 2> *> >(), py::arg("nodes"), py::arg("faces"), py::arg("vertexElements"))
        .def(py::init< >())
        .def(
            "GetElementIteratorBegin",
            (::VertexMesh<2, 2>::VertexElementIterator(VertexMesh2_2::*)(bool)) &VertexMesh2_2::GetElementIteratorBegin,
            " " , py::arg("skipDeletedElements") = true )
        .def(
            "GetElementIteratorEnd",
            (::VertexMesh<2, 2>::VertexElementIterator(VertexMesh2_2::*)()) &VertexMesh2_2::GetElementIteratorEnd,
            " "  )
        .def(
            "GetNumEdges",
            (unsigned int(VertexMesh2_2::*)() const ) &VertexMesh2_2::GetNumEdges,
            " "  )
        .def(
            "GetEdge",
            (::Edge<2> *(VertexMesh2_2::*)(unsigned int) const ) &VertexMesh2_2::GetEdge,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNumNodes",
            (unsigned int(VertexMesh2_2::*)() const ) &VertexMesh2_2::GetNumNodes,
            " "  )
        .def(
            "GetNumElements",
            (unsigned int(VertexMesh2_2::*)() const ) &VertexMesh2_2::GetNumElements,
            " "  )
        .def(
            "GetNumAllElements",
            (unsigned int(VertexMesh2_2::*)() const ) &VertexMesh2_2::GetNumAllElements,
            " "  )
        .def(
            "GetNumFaces",
            (unsigned int(VertexMesh2_2::*)() const ) &VertexMesh2_2::GetNumFaces,
            " "  )
        .def(
            "GetElement",
            (::VertexElement<2, 2> *(VertexMesh2_2::*)(unsigned int) const ) &VertexMesh2_2::GetElement,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetCentroidOfElement",
            (::boost::numeric::ublas::c_vector<double, 2>(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::GetCentroidOfElement,
            " " , py::arg("index") )
        .def(
            "ConstructFromMeshReader",
            (void(VertexMesh2_2::*)(::AbstractMeshReader<2, 2> &)) &VertexMesh2_2::ConstructFromMeshReader,
            " " , py::arg("rMeshReader") )
        .def(
            "Clear",
            (void(VertexMesh2_2::*)()) &VertexMesh2_2::Clear,
            " "  )
        .def(
            "GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex",
            (unsigned int(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex,
            " " , py::arg("elementIndex") )
        .def(
            "GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex",
            (unsigned int(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex,
            " " , py::arg("nodeIndex") )
        .def(
            "GetRosetteRankOfElement",
            (unsigned int(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::GetRosetteRankOfElement,
            " " , py::arg("index") )
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(VertexMesh2_2::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &VertexMesh2_2::GetVectorFromAtoB,
            " " , py::arg("rLocationA"), py::arg("rLocationB") )
        .def(
            "GetVolumeOfElement",
            (double(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::GetVolumeOfElement,
            " " , py::arg("index") )
        .def(
            "GetSurfaceAreaOfElement",
            (double(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::GetSurfaceAreaOfElement,
            " " , py::arg("index") )
        .def(
            "GetAreaGradientOfElementAtNode",
            (::boost::numeric::ublas::c_vector<double, 2>(VertexMesh2_2::*)(::VertexElement<2, 2> *, unsigned int)) &VertexMesh2_2::GetAreaGradientOfElementAtNode,
            " " , py::arg("pElement"), py::arg("localIndex") )
        .def(
            "GetPreviousEdgeGradientOfElementAtNode",
            (::boost::numeric::ublas::c_vector<double, 2>(VertexMesh2_2::*)(::VertexElement<2, 2> *, unsigned int)) &VertexMesh2_2::GetPreviousEdgeGradientOfElementAtNode,
            " " , py::arg("pElement"), py::arg("localIndex") )
        .def(
            "GetNextEdgeGradientOfElementAtNode",
            (::boost::numeric::ublas::c_vector<double, 2>(VertexMesh2_2::*)(::VertexElement<2, 2> *, unsigned int)) &VertexMesh2_2::GetNextEdgeGradientOfElementAtNode,
            " " , py::arg("pElement"), py::arg("localIndex") )
        .def(
            "GetPerimeterGradientOfElementAtNode",
            (::boost::numeric::ublas::c_vector<double, 2>(VertexMesh2_2::*)(::VertexElement<2, 2> *, unsigned int)) &VertexMesh2_2::GetPerimeterGradientOfElementAtNode,
            " " , py::arg("pElement"), py::arg("localIndex") )
        .def(
            "CalculateMomentsOfElement",
            (::boost::numeric::ublas::c_vector<double, 3>(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::CalculateMomentsOfElement,
            " " , py::arg("index") )
        .def(
            "GetEdgeLength",
            (double(VertexMesh2_2::*)(unsigned int, unsigned int)) &VertexMesh2_2::GetEdgeLength,
            " " , py::arg("elementIndex1"), py::arg("elementIndex2") )
        .def(
            "GetElongationShapeFactorOfElement",
            (double(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::GetElongationShapeFactorOfElement,
            " " , py::arg("elementIndex") )
        .def(
            "CalculateUnitNormalToFaceWithArea",
            (double(VertexMesh2_2::*)(::VertexElement<1, 2> *, ::boost::numeric::ublas::c_vector<double, 2> &)) &VertexMesh2_2::CalculateUnitNormalToFaceWithArea,
            " " , py::arg("pFace"), py::arg("rNormal") )
        .def(
            "CalculateAreaOfFace",
            (double(VertexMesh2_2::*)(::VertexElement<1, 2> *)) &VertexMesh2_2::CalculateAreaOfFace,
            " " , py::arg("pFace") )
        .def(
            "GetShortAxisOfElement",
            (::boost::numeric::ublas::c_vector<double, 2>(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::GetShortAxisOfElement,
            " " , py::arg("index") )
        .def(
            "GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::GetNeighbouringNodeIndices,
            " " , py::arg("nodeIndex") )
        .def(
            "GetNeighbouringNodeNotAlsoInElement",
            (::std::set<unsigned int>(VertexMesh2_2::*)(unsigned int, unsigned int)) &VertexMesh2_2::GetNeighbouringNodeNotAlsoInElement,
            " " , py::arg("nodeIndex"), py::arg("elemIndex") )
        .def(
            "GetNeighbouringElementIndices",
            (::std::set<unsigned int>(VertexMesh2_2::*)(unsigned int)) &VertexMesh2_2::GetNeighbouringElementIndices,
            " " , py::arg("elementIndex") )
        .def(
            "GetMeshForVtk",
            (::VertexMesh<2, 2> *(VertexMesh2_2::*)()) &VertexMesh2_2::GetMeshForVtk,
            " "  , py::return_value_policy::reference)
    ;
}
