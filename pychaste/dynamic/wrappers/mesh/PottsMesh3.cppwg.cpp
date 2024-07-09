#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PottsMesh.hpp"

#include "PottsMesh3.cppwg.hpp"

namespace py = pybind11;
typedef PottsMesh<3 > PottsMesh3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef unsigned int unsignedint;

class PottsMesh3_Overrides : public PottsMesh3{
    public:
    using PottsMesh3::PottsMesh;
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            PottsMesh3,
            GetNumNodes,
            );
    }
    unsigned int GetNumElements() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            PottsMesh3,
            GetNumElements,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 3> GetCentroidOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            PottsMesh3,
            GetCentroidOfElement,
                    index);
    }
    void Clear() override {
        PYBIND11_OVERRIDE(
            void,
            PottsMesh3,
            Clear,
            );
    }
    double GetVolumeOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            double,
            PottsMesh3,
            GetVolumeOfElement,
                    index);
    }
    double GetSurfaceAreaOfElement(unsigned int index) override {
        PYBIND11_OVERRIDE(
            double,
            PottsMesh3,
            GetSurfaceAreaOfElement,
                    index);
    }
    unsigned int SolveNodeMapping(unsigned int index) const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            PottsMesh3,
            SolveNodeMapping,
                    index);
    }

};
void register_PottsMesh3_class(py::module &m){
py::class_<PottsMesh3 , PottsMesh3_Overrides , boost::shared_ptr<PottsMesh3 >  , AbstractMesh<3, 3>  >(m, "PottsMesh3")
        .def(py::init<::std::vector<Node<3> *>, ::std::vector<PottsElement<3> *>, ::std::vector<std::set<unsigned int>>, ::std::vector<std::set<unsigned int>> >(), py::arg("nodes"), py::arg("pottsElements"), py::arg("vonNeumannNeighbouringNodeIndices"), py::arg("mooreNeighbouringNodeIndices"))
        .def(py::init< >())
        .def(
            "GetElementIteratorBegin",
            (::PottsMesh<3>::PottsElementIterator(PottsMesh3::*)(bool)) &PottsMesh3::GetElementIteratorBegin,
            " " , py::arg("skipDeletedElements") = true )
        .def(
            "GetElementIteratorEnd",
            (::PottsMesh<3>::PottsElementIterator(PottsMesh3::*)()) &PottsMesh3::GetElementIteratorEnd,
            " "  )
        .def(
            "GetNumNodes",
            (unsigned int(PottsMesh3::*)() const ) &PottsMesh3::GetNumNodes,
            " "  )
        .def(
            "GetNumElements",
            (unsigned int(PottsMesh3::*)() const ) &PottsMesh3::GetNumElements,
            " "  )
        .def(
            "GetNumAllElements",
            (unsigned int(PottsMesh3::*)() const ) &PottsMesh3::GetNumAllElements,
            " "  )
        .def(
            "GetElement",
            (::PottsElement<3> *(PottsMesh3::*)(unsigned int) const ) &PottsMesh3::GetElement,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetCentroidOfElement",
            (::boost::numeric::ublas::c_vector<double, 3>(PottsMesh3::*)(unsigned int)) &PottsMesh3::GetCentroidOfElement,
            " " , py::arg("index") )
        .def(
            "ConstructFromMeshReader",
            (void(PottsMesh3::*)(::AbstractMeshReader<3, 3> &)) &PottsMesh3::ConstructFromMeshReader,
            " " , py::arg("rMeshReader") )
        .def(
            "Clear",
            (void(PottsMesh3::*)()) &PottsMesh3::Clear,
            " "  )
        .def(
            "GetVolumeOfElement",
            (double(PottsMesh3::*)(unsigned int)) &PottsMesh3::GetVolumeOfElement,
            " " , py::arg("index") )
        .def(
            "GetSurfaceAreaOfElement",
            (double(PottsMesh3::*)(unsigned int)) &PottsMesh3::GetSurfaceAreaOfElement,
            " " , py::arg("index") )
        .def(
            "GetMooreNeighbouringNodeIndices",
            (::std::set<unsigned int>(PottsMesh3::*)(unsigned int)) &PottsMesh3::GetMooreNeighbouringNodeIndices,
            " " , py::arg("nodeIndex") )
        .def(
            "GetVonNeumannNeighbouringNodeIndices",
            (::std::set<unsigned int>(PottsMesh3::*)(unsigned int)) &PottsMesh3::GetVonNeumannNeighbouringNodeIndices,
            " " , py::arg("nodeIndex") )
        .def(
            "DeleteNode",
            (void(PottsMesh3::*)(unsigned int)) &PottsMesh3::DeleteNode,
            " " , py::arg("index") )
        .def(
            "DeleteElement",
            (void(PottsMesh3::*)(unsigned int)) &PottsMesh3::DeleteElement,
            " " , py::arg("index") )
        .def(
            "RemoveDeletedElements",
            (void(PottsMesh3::*)()) &PottsMesh3::RemoveDeletedElements,
            " "  )
        .def(
            "DivideElement",
            (unsigned int(PottsMesh3::*)(::PottsElement<3> *, bool)) &PottsMesh3::DivideElement,
            " " , py::arg("pElement"), py::arg("placeOriginalElementBelow") = false )
        .def(
            "AddElement",
            (unsigned int(PottsMesh3::*)(::PottsElement<3> *)) &PottsMesh3::AddElement,
            " " , py::arg("pNewElement") )
        .def(
            "GetNeighbouringElementIndices",
            (::std::set<unsigned int>(PottsMesh3::*)(unsigned int)) &PottsMesh3::GetNeighbouringElementIndices,
            " " , py::arg("elementIndex") )
    ;
}
