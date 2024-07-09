#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractMesh.hpp"

#include "AbstractMesh3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractMesh<3,3 > AbstractMesh3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::Node<3> * _Node_lt_3_gt_Ptr;
typedef ::DistributedVectorFactory * _DistributedVectorFactoryPtr;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::ChasteCuboid<3> _ChasteCuboid_lt_3_gt_;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class AbstractMesh3_3_Overrides : public AbstractMesh3_3{
    public:
    using AbstractMesh3_3::AbstractMesh;
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractMesh3_3,
            GetNumNodes,
            );
    }
    unsigned int GetNumAllNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractMesh3_3,
            GetNumAllNodes,
            );
    }
    ::Node<3> * GetNodeOrHaloNode(unsigned int index) const  override {
        PYBIND11_OVERRIDE(
            _Node_lt_3_gt_Ptr,
            AbstractMesh3_3,
            GetNodeOrHaloNode,
                    index);
    }
    void ReadNodesPerProcessorFile(::std::string const & rNodesPerProcessorFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh3_3,
            ReadNodesPerProcessorFile,
                    rNodesPerProcessorFile);
    }
    ::DistributedVectorFactory * GetDistributedVectorFactory() override {
        PYBIND11_OVERRIDE(
            _DistributedVectorFactoryPtr,
            AbstractMesh3_3,
            GetDistributedVectorFactory,
            );
    }
    void SetDistributedVectorFactory(::DistributedVectorFactory * pFactory) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh3_3,
            SetDistributedVectorFactory,
                    pFactory);
    }
    void PermuteNodes() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh3_3,
            PermuteNodes,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 3> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 3> const & rLocationA, ::boost::numeric::ublas::c_vector<double, 3> const & rLocationB) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            AbstractMesh3_3,
            GetVectorFromAtoB,
                    rLocationA,
        rLocationB);
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            AbstractMesh3_3,
            GetWidth,
                    rDimension);
    }
    ::ChasteCuboid<3> CalculateBoundingBox() const  override {
        PYBIND11_OVERRIDE(
            _ChasteCuboid_lt_3_gt_,
            AbstractMesh3_3,
            CalculateBoundingBox,
            );
    }
    unsigned int GetNearestNodeIndex(::ChastePoint<3> const & rTestPoint) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractMesh3_3,
            GetNearestNodeIndex,
                    rTestPoint);
    }
    void Scale(double const xFactor, double const yFactor, double const zFactor) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh3_3,
            Scale,
                    xFactor,
        yFactor,
        zFactor);
    }
    void Translate(::boost::numeric::ublas::c_vector<double, 3> const & rDisplacement) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh3_3,
            Translate,
                    rDisplacement);
    }
    void Rotate(::boost::numeric::ublas::c_matrix<double, 3, 3> rotationMatrix) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh3_3,
            Rotate,
                    rotationMatrix);
    }
    void RefreshMesh() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh3_3,
            RefreshMesh,
            );
    }
    void SetElementOwnerships() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh3_3,
            SetElementOwnerships,
            );
    }

};
void register_AbstractMesh3_3_class(py::module &m){
py::class_<AbstractMesh3_3 , AbstractMesh3_3_Overrides , boost::shared_ptr<AbstractMesh3_3 >   >(m, "AbstractMesh3_3")
        .def(
            "GetNodeIteratorBegin",
            (::AbstractMesh<3, 3>::NodeIterator(AbstractMesh3_3::*)(bool)) &AbstractMesh3_3::GetNodeIteratorBegin,
            " " , py::arg("skipDeletedNodes") = true )
        .def(
            "GetNodeIteratorEnd",
            (::AbstractMesh<3, 3>::NodeIterator(AbstractMesh3_3::*)()) &AbstractMesh3_3::GetNodeIteratorEnd,
            " "  )
        .def(
            "GetNumNodes",
            (unsigned int(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::GetNumNodes,
            " "  )
        .def(
            "GetNumBoundaryNodes",
            (unsigned int(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::GetNumBoundaryNodes,
            " "  )
        .def(
            "GetNumAllNodes",
            (unsigned int(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::GetNumAllNodes,
            " "  )
        .def(
            "GetNumNodeAttributes",
            (unsigned int(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::GetNumNodeAttributes,
            " "  )
        .def(
            "GetNode",
            (::Node<3> *(AbstractMesh3_3::*)(unsigned int) const ) &AbstractMesh3_3::GetNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNodeOrHaloNode",
            (::Node<3> *(AbstractMesh3_3::*)(unsigned int) const ) &AbstractMesh3_3::GetNodeOrHaloNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNodeFromPrePermutationIndex",
            (::Node<3> *(AbstractMesh3_3::*)(unsigned int) const ) &AbstractMesh3_3::GetNodeFromPrePermutationIndex,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "ReadNodesPerProcessorFile",
            (void(AbstractMesh3_3::*)(::std::string const &)) &AbstractMesh3_3::ReadNodesPerProcessorFile,
            " " , py::arg("rNodesPerProcessorFile") )
        .def(
            "GetDistributedVectorFactory",
            (::DistributedVectorFactory *(AbstractMesh3_3::*)()) &AbstractMesh3_3::GetDistributedVectorFactory,
            " "  , py::return_value_policy::reference)
        .def(
            "SetDistributedVectorFactory",
            (void(AbstractMesh3_3::*)(::DistributedVectorFactory *)) &AbstractMesh3_3::SetDistributedVectorFactory,
            " " , py::arg("pFactory") )
        .def(
            "PermuteNodes",
            (void(AbstractMesh3_3::*)()) &AbstractMesh3_3::PermuteNodes,
            " "  )
        .def(
            "GetBoundaryNodeIteratorBegin",
            (::AbstractMesh<3, 3>::BoundaryNodeIterator(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::GetBoundaryNodeIteratorBegin,
            " "  )
        .def(
            "GetBoundaryNodeIteratorEnd",
            (::AbstractMesh<3, 3>::BoundaryNodeIterator(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::GetBoundaryNodeIteratorEnd,
            " "  )
        .def(
            "GetMeshFileBaseName",
            (::std::string(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::GetMeshFileBaseName,
            " "  )
        .def(
            "IsMeshOnDisk",
            (bool(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::IsMeshOnDisk,
            " "  )
        .def(
            "rGetNodePermutation",
            (::std::vector<unsigned int> const &(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::rGetNodePermutation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractMesh3_3::*)(::boost::numeric::ublas::c_vector<double, 3> const &, ::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractMesh3_3::GetVectorFromAtoB,
            " " , py::arg("rLocationA"), py::arg("rLocationB") )
        .def(
            "GetDistanceBetweenNodes",
            (double(AbstractMesh3_3::*)(unsigned int, unsigned int)) &AbstractMesh3_3::GetDistanceBetweenNodes,
            " " , py::arg("indexA"), py::arg("indexB") )
        .def(
            "GetWidth",
            (double(AbstractMesh3_3::*)(unsigned int const &) const ) &AbstractMesh3_3::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "CalculateBoundingBox",
            (::ChasteCuboid<3>(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::CalculateBoundingBox,
            " "  )
        .def(
            "GetNearestNodeIndex",
            (unsigned int(AbstractMesh3_3::*)(::ChastePoint<3> const &)) &AbstractMesh3_3::GetNearestNodeIndex,
            " " , py::arg("rTestPoint") )
        .def(
            "Scale",
            (void(AbstractMesh3_3::*)(double const, double const, double const)) &AbstractMesh3_3::Scale,
            " " , py::arg("xFactor") = 1., py::arg("yFactor") = 1., py::arg("zFactor") = 1. )
        .def(
            "Translate",
            (void(AbstractMesh3_3::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractMesh3_3::Translate,
            " " , py::arg("rDisplacement") )
        .def(
            "Translate",
            (void(AbstractMesh3_3::*)(double const, double const, double const)) &AbstractMesh3_3::Translate,
            " " , py::arg("xMovement") = 0., py::arg("yMovement") = 0., py::arg("zMovement") = 0. )
        .def(
            "Rotate",
            (void(AbstractMesh3_3::*)(::boost::numeric::ublas::c_matrix<double, 3, 3>)) &AbstractMesh3_3::Rotate,
            " " , py::arg("rotationMatrix") )
        .def(
            "Rotate",
            (void(AbstractMesh3_3::*)(::boost::numeric::ublas::c_vector<double, 3>, double)) &AbstractMesh3_3::Rotate,
            " " , py::arg("axis"), py::arg("angle") )
        .def(
            "RotateX",
            (void(AbstractMesh3_3::*)(double const)) &AbstractMesh3_3::RotateX,
            " " , py::arg("theta") )
        .def(
            "RotateY",
            (void(AbstractMesh3_3::*)(double const)) &AbstractMesh3_3::RotateY,
            " " , py::arg("theta") )
        .def(
            "RotateZ",
            (void(AbstractMesh3_3::*)(double const)) &AbstractMesh3_3::RotateZ,
            " " , py::arg("theta") )
        .def(
            "Rotate",
            (void(AbstractMesh3_3::*)(double)) &AbstractMesh3_3::Rotate,
            " " , py::arg("theta") )
        .def(
            "RefreshMesh",
            (void(AbstractMesh3_3::*)()) &AbstractMesh3_3::RefreshMesh,
            " "  )
        .def(
            "IsMeshChanging",
            (bool(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::IsMeshChanging,
            " "  )
        .def(
            "CalculateMaximumContainingElementsPerProcess",
            (unsigned int(AbstractMesh3_3::*)() const ) &AbstractMesh3_3::CalculateMaximumContainingElementsPerProcess,
            " "  )
        .def(
            "SetMeshHasChangedSinceLoading",
            (void(AbstractMesh3_3::*)()) &AbstractMesh3_3::SetMeshHasChangedSinceLoading,
            " "  )
    ;
}
