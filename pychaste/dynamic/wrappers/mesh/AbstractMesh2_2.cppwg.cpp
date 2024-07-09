#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractMesh.hpp"

#include "AbstractMesh2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractMesh<2,2 > AbstractMesh2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::Node<2> * _Node_lt_2_gt_Ptr;
typedef ::DistributedVectorFactory * _DistributedVectorFactoryPtr;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::ChasteCuboid<2> _ChasteCuboid_lt_2_gt_;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class AbstractMesh2_2_Overrides : public AbstractMesh2_2{
    public:
    using AbstractMesh2_2::AbstractMesh;
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractMesh2_2,
            GetNumNodes,
            );
    }
    unsigned int GetNumAllNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractMesh2_2,
            GetNumAllNodes,
            );
    }
    ::Node<2> * GetNodeOrHaloNode(unsigned int index) const  override {
        PYBIND11_OVERRIDE(
            _Node_lt_2_gt_Ptr,
            AbstractMesh2_2,
            GetNodeOrHaloNode,
                    index);
    }
    void ReadNodesPerProcessorFile(::std::string const & rNodesPerProcessorFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh2_2,
            ReadNodesPerProcessorFile,
                    rNodesPerProcessorFile);
    }
    ::DistributedVectorFactory * GetDistributedVectorFactory() override {
        PYBIND11_OVERRIDE(
            _DistributedVectorFactoryPtr,
            AbstractMesh2_2,
            GetDistributedVectorFactory,
            );
    }
    void SetDistributedVectorFactory(::DistributedVectorFactory * pFactory) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh2_2,
            SetDistributedVectorFactory,
                    pFactory);
    }
    void PermuteNodes() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh2_2,
            PermuteNodes,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocationA, ::boost::numeric::ublas::c_vector<double, 2> const & rLocationB) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            AbstractMesh2_2,
            GetVectorFromAtoB,
                    rLocationA,
        rLocationB);
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            AbstractMesh2_2,
            GetWidth,
                    rDimension);
    }
    ::ChasteCuboid<2> CalculateBoundingBox() const  override {
        PYBIND11_OVERRIDE(
            _ChasteCuboid_lt_2_gt_,
            AbstractMesh2_2,
            CalculateBoundingBox,
            );
    }
    unsigned int GetNearestNodeIndex(::ChastePoint<2> const & rTestPoint) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractMesh2_2,
            GetNearestNodeIndex,
                    rTestPoint);
    }
    void Scale(double const xFactor, double const yFactor, double const zFactor) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh2_2,
            Scale,
                    xFactor,
        yFactor,
        zFactor);
    }
    void Translate(::boost::numeric::ublas::c_vector<double, 2> const & rDisplacement) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh2_2,
            Translate,
                    rDisplacement);
    }
    void Rotate(::boost::numeric::ublas::c_matrix<double, 2, 2> rotationMatrix) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh2_2,
            Rotate,
                    rotationMatrix);
    }
    void RefreshMesh() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh2_2,
            RefreshMesh,
            );
    }
    void SetElementOwnerships() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh2_2,
            SetElementOwnerships,
            );
    }

};
void register_AbstractMesh2_2_class(py::module &m){
py::class_<AbstractMesh2_2 , AbstractMesh2_2_Overrides , boost::shared_ptr<AbstractMesh2_2 >   >(m, "AbstractMesh2_2")
        .def(
            "GetNodeIteratorBegin",
            (::AbstractMesh<2, 2>::NodeIterator(AbstractMesh2_2::*)(bool)) &AbstractMesh2_2::GetNodeIteratorBegin,
            " " , py::arg("skipDeletedNodes") = true )
        .def(
            "GetNodeIteratorEnd",
            (::AbstractMesh<2, 2>::NodeIterator(AbstractMesh2_2::*)()) &AbstractMesh2_2::GetNodeIteratorEnd,
            " "  )
        .def(
            "GetNumNodes",
            (unsigned int(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::GetNumNodes,
            " "  )
        .def(
            "GetNumBoundaryNodes",
            (unsigned int(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::GetNumBoundaryNodes,
            " "  )
        .def(
            "GetNumAllNodes",
            (unsigned int(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::GetNumAllNodes,
            " "  )
        .def(
            "GetNumNodeAttributes",
            (unsigned int(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::GetNumNodeAttributes,
            " "  )
        .def(
            "GetNode",
            (::Node<2> *(AbstractMesh2_2::*)(unsigned int) const ) &AbstractMesh2_2::GetNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNodeOrHaloNode",
            (::Node<2> *(AbstractMesh2_2::*)(unsigned int) const ) &AbstractMesh2_2::GetNodeOrHaloNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNodeFromPrePermutationIndex",
            (::Node<2> *(AbstractMesh2_2::*)(unsigned int) const ) &AbstractMesh2_2::GetNodeFromPrePermutationIndex,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "ReadNodesPerProcessorFile",
            (void(AbstractMesh2_2::*)(::std::string const &)) &AbstractMesh2_2::ReadNodesPerProcessorFile,
            " " , py::arg("rNodesPerProcessorFile") )
        .def(
            "GetDistributedVectorFactory",
            (::DistributedVectorFactory *(AbstractMesh2_2::*)()) &AbstractMesh2_2::GetDistributedVectorFactory,
            " "  , py::return_value_policy::reference)
        .def(
            "SetDistributedVectorFactory",
            (void(AbstractMesh2_2::*)(::DistributedVectorFactory *)) &AbstractMesh2_2::SetDistributedVectorFactory,
            " " , py::arg("pFactory") )
        .def(
            "PermuteNodes",
            (void(AbstractMesh2_2::*)()) &AbstractMesh2_2::PermuteNodes,
            " "  )
        .def(
            "GetBoundaryNodeIteratorBegin",
            (::AbstractMesh<2, 2>::BoundaryNodeIterator(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::GetBoundaryNodeIteratorBegin,
            " "  )
        .def(
            "GetBoundaryNodeIteratorEnd",
            (::AbstractMesh<2, 2>::BoundaryNodeIterator(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::GetBoundaryNodeIteratorEnd,
            " "  )
        .def(
            "GetMeshFileBaseName",
            (::std::string(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::GetMeshFileBaseName,
            " "  )
        .def(
            "IsMeshOnDisk",
            (bool(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::IsMeshOnDisk,
            " "  )
        .def(
            "rGetNodePermutation",
            (::std::vector<unsigned int> const &(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::rGetNodePermutation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(AbstractMesh2_2::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &AbstractMesh2_2::GetVectorFromAtoB,
            " " , py::arg("rLocationA"), py::arg("rLocationB") )
        .def(
            "GetDistanceBetweenNodes",
            (double(AbstractMesh2_2::*)(unsigned int, unsigned int)) &AbstractMesh2_2::GetDistanceBetweenNodes,
            " " , py::arg("indexA"), py::arg("indexB") )
        .def(
            "GetWidth",
            (double(AbstractMesh2_2::*)(unsigned int const &) const ) &AbstractMesh2_2::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "CalculateBoundingBox",
            (::ChasteCuboid<2>(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::CalculateBoundingBox,
            " "  )
        .def(
            "GetNearestNodeIndex",
            (unsigned int(AbstractMesh2_2::*)(::ChastePoint<2> const &)) &AbstractMesh2_2::GetNearestNodeIndex,
            " " , py::arg("rTestPoint") )
        .def(
            "Scale",
            (void(AbstractMesh2_2::*)(double const, double const, double const)) &AbstractMesh2_2::Scale,
            " " , py::arg("xFactor") = 1., py::arg("yFactor") = 1., py::arg("zFactor") = 1. )
        .def(
            "Translate",
            (void(AbstractMesh2_2::*)(::boost::numeric::ublas::c_vector<double, 2> const &)) &AbstractMesh2_2::Translate,
            " " , py::arg("rDisplacement") )
        .def(
            "Translate",
            (void(AbstractMesh2_2::*)(double const, double const, double const)) &AbstractMesh2_2::Translate,
            " " , py::arg("xMovement") = 0., py::arg("yMovement") = 0., py::arg("zMovement") = 0. )
        .def(
            "Rotate",
            (void(AbstractMesh2_2::*)(::boost::numeric::ublas::c_matrix<double, 2, 2>)) &AbstractMesh2_2::Rotate,
            " " , py::arg("rotationMatrix") )
        .def(
            "Rotate",
            (void(AbstractMesh2_2::*)(::boost::numeric::ublas::c_vector<double, 3>, double)) &AbstractMesh2_2::Rotate,
            " " , py::arg("axis"), py::arg("angle") )
        .def(
            "RotateX",
            (void(AbstractMesh2_2::*)(double const)) &AbstractMesh2_2::RotateX,
            " " , py::arg("theta") )
        .def(
            "RotateY",
            (void(AbstractMesh2_2::*)(double const)) &AbstractMesh2_2::RotateY,
            " " , py::arg("theta") )
        .def(
            "RotateZ",
            (void(AbstractMesh2_2::*)(double const)) &AbstractMesh2_2::RotateZ,
            " " , py::arg("theta") )
        .def(
            "Rotate",
            (void(AbstractMesh2_2::*)(double)) &AbstractMesh2_2::Rotate,
            " " , py::arg("theta") )
        .def(
            "RefreshMesh",
            (void(AbstractMesh2_2::*)()) &AbstractMesh2_2::RefreshMesh,
            " "  )
        .def(
            "IsMeshChanging",
            (bool(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::IsMeshChanging,
            " "  )
        .def(
            "CalculateMaximumContainingElementsPerProcess",
            (unsigned int(AbstractMesh2_2::*)() const ) &AbstractMesh2_2::CalculateMaximumContainingElementsPerProcess,
            " "  )
        .def(
            "SetMeshHasChangedSinceLoading",
            (void(AbstractMesh2_2::*)()) &AbstractMesh2_2::SetMeshHasChangedSinceLoading,
            " "  )
    ;
}
