/*

Copyright (c) 2005-2024, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

// This file is auto-generated; manual changes will be overwritten.
// To make enduring changes, see pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractMesh.hpp"

#include "AbstractMesh_2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractMesh<2, 2> AbstractMesh_2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::Node<2> * _Node_lt_2_gt_Ptr;
typedef ::DistributedVectorFactory * _DistributedVectorFactoryPtr;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::ChasteCuboid<2> _ChasteCuboid_lt_2_gt_;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class AbstractMesh_2_2_Overrides : public AbstractMesh_2_2
{
public:
    using AbstractMesh_2_2::AbstractMesh;
    unsigned int GetNumNodes() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractMesh_2_2,
            GetNumNodes,
            );
    }
    unsigned int GetNumAllNodes() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractMesh_2_2,
            GetNumAllNodes,
            );
    }
    ::Node<2> * GetNodeOrHaloNode(unsigned int index) const override
    {
        PYBIND11_OVERRIDE(
            _Node_lt_2_gt_Ptr,
            AbstractMesh_2_2,
            GetNodeOrHaloNode,
            index);
    }
    void ReadNodesPerProcessorFile(::std::string const & rNodesPerProcessorFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh_2_2,
            ReadNodesPerProcessorFile,
            rNodesPerProcessorFile);
    }
    ::DistributedVectorFactory * GetDistributedVectorFactory() override
    {
        PYBIND11_OVERRIDE(
            _DistributedVectorFactoryPtr,
            AbstractMesh_2_2,
            GetDistributedVectorFactory,
            );
    }
    void SetDistributedVectorFactory(::DistributedVectorFactory * pFactory) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh_2_2,
            SetDistributedVectorFactory,
            pFactory);
    }
    void PermuteNodes() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh_2_2,
            PermuteNodes,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocationA, ::boost::numeric::ublas::c_vector<double, 2> const & rLocationB) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            AbstractMesh_2_2,
            GetVectorFromAtoB,
            rLocationA,
            rLocationB);
    }
    double GetWidth(unsigned int const & rDimension) const override
    {
        PYBIND11_OVERRIDE(
            double,
            AbstractMesh_2_2,
            GetWidth,
            rDimension);
    }
    ::ChasteCuboid<2> CalculateBoundingBox() const override
    {
        PYBIND11_OVERRIDE(
            _ChasteCuboid_lt_2_gt_,
            AbstractMesh_2_2,
            CalculateBoundingBox,
            );
    }
    unsigned int GetNearestNodeIndex(::ChastePoint<2> const & rTestPoint) override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractMesh_2_2,
            GetNearestNodeIndex,
            rTestPoint);
    }
    void Scale(double const xFactor, double const yFactor, double const zFactor) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh_2_2,
            Scale,
            xFactor,
            yFactor,
            zFactor);
    }
    void Translate(::boost::numeric::ublas::c_vector<double, 2> const & rDisplacement) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh_2_2,
            Translate,
            rDisplacement);
    }
    void Rotate(::boost::numeric::ublas::c_matrix<double, 2, 2> rotationMatrix) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh_2_2,
            Rotate,
            rotationMatrix);
    }
    void RefreshMesh() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh_2_2,
            RefreshMesh,
            );
    }
    void SetElementOwnerships() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractMesh_2_2,
            SetElementOwnerships,
            );
    }
};

void register_AbstractMesh_2_2_class(py::module &m)
{
    py::class_<AbstractMesh_2_2, AbstractMesh_2_2_Overrides, boost::shared_ptr<AbstractMesh_2_2>>(m, "AbstractMesh_2_2")
        .def("GetNodeIteratorBegin",
            (::AbstractMesh<2, 2>::NodeIterator(AbstractMesh_2_2::*)(bool)) &AbstractMesh_2_2::GetNodeIteratorBegin,
            " ", py::arg("skipDeletedNodes") = true)
        .def("GetNodeIteratorEnd",
            (::AbstractMesh<2, 2>::NodeIterator(AbstractMesh_2_2::*)()) &AbstractMesh_2_2::GetNodeIteratorEnd,
            " ")
        .def("GetNumNodes",
            (unsigned int(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::GetNumNodes,
            " ")
        .def("GetNumBoundaryNodes",
            (unsigned int(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::GetNumBoundaryNodes,
            " ")
        .def("GetNumAllNodes",
            (unsigned int(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::GetNumAllNodes,
            " ")
        .def("GetNumNodeAttributes",
            (unsigned int(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::GetNumNodeAttributes,
            " ")
        .def("GetNode",
            (::Node<2> *(AbstractMesh_2_2::*)(unsigned int) const) &AbstractMesh_2_2::GetNode,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("GetNodeOrHaloNode",
            (::Node<2> *(AbstractMesh_2_2::*)(unsigned int) const) &AbstractMesh_2_2::GetNodeOrHaloNode,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("GetNodeFromPrePermutationIndex",
            (::Node<2> *(AbstractMesh_2_2::*)(unsigned int) const) &AbstractMesh_2_2::GetNodeFromPrePermutationIndex,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("ReadNodesPerProcessorFile",
            (void(AbstractMesh_2_2::*)(::std::string const &)) &AbstractMesh_2_2::ReadNodesPerProcessorFile,
            " ", py::arg("rNodesPerProcessorFile"))
        .def("GetDistributedVectorFactory",
            (::DistributedVectorFactory *(AbstractMesh_2_2::*)()) &AbstractMesh_2_2::GetDistributedVectorFactory,
            " ", py::return_value_policy::reference)
        .def("SetDistributedVectorFactory",
            (void(AbstractMesh_2_2::*)(::DistributedVectorFactory *)) &AbstractMesh_2_2::SetDistributedVectorFactory,
            " ", py::arg("pFactory"))
        .def("PermuteNodes",
            (void(AbstractMesh_2_2::*)()) &AbstractMesh_2_2::PermuteNodes,
            " ")
        .def("GetBoundaryNodeIteratorBegin",
            (::AbstractMesh<2, 2>::BoundaryNodeIterator(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::GetBoundaryNodeIteratorBegin,
            " ")
        .def("GetBoundaryNodeIteratorEnd",
            (::AbstractMesh<2, 2>::BoundaryNodeIterator(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::GetBoundaryNodeIteratorEnd,
            " ")
        .def("GetMeshFileBaseName",
            (::std::string(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::GetMeshFileBaseName,
            " ")
        .def("IsMeshOnDisk",
            (bool(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::IsMeshOnDisk,
            " ")
        .def("rGetNodePermutation",
            (::std::vector<unsigned int> const &(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::rGetNodePermutation,
            " ", py::return_value_policy::reference_internal)
        .def("GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(AbstractMesh_2_2::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &AbstractMesh_2_2::GetVectorFromAtoB,
            " ", py::arg("rLocationA"), py::arg("rLocationB"))
        .def("GetDistanceBetweenNodes",
            (double(AbstractMesh_2_2::*)(unsigned int, unsigned int)) &AbstractMesh_2_2::GetDistanceBetweenNodes,
            " ", py::arg("indexA"), py::arg("indexB"))
        .def("GetWidth",
            (double(AbstractMesh_2_2::*)(unsigned int const &) const) &AbstractMesh_2_2::GetWidth,
            " ", py::arg("rDimension"))
        .def("CalculateBoundingBox",
            (::ChasteCuboid<2>(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::CalculateBoundingBox,
            " ")
        .def("GetNearestNodeIndex",
            (unsigned int(AbstractMesh_2_2::*)(::ChastePoint<2> const &)) &AbstractMesh_2_2::GetNearestNodeIndex,
            " ", py::arg("rTestPoint"))
        .def("Scale",
            (void(AbstractMesh_2_2::*)(double const, double const, double const)) &AbstractMesh_2_2::Scale,
            " ", py::arg("xFactor") = 1., py::arg("yFactor") = 1., py::arg("zFactor") = 1.)
        .def("Translate",
            (void(AbstractMesh_2_2::*)(::boost::numeric::ublas::c_vector<double, 2> const &)) &AbstractMesh_2_2::Translate,
            " ", py::arg("rDisplacement"))
        .def("Translate",
            (void(AbstractMesh_2_2::*)(double const, double const, double const)) &AbstractMesh_2_2::Translate,
            " ", py::arg("xMovement") = 0., py::arg("yMovement") = 0., py::arg("zMovement") = 0.)
        .def("Rotate",
            (void(AbstractMesh_2_2::*)(::boost::numeric::ublas::c_matrix<double, 2, 2>)) &AbstractMesh_2_2::Rotate,
            " ", py::arg("rotationMatrix"))
        .def("Rotate",
            (void(AbstractMesh_2_2::*)(::boost::numeric::ublas::c_vector<double, 3>, double)) &AbstractMesh_2_2::Rotate,
            " ", py::arg("axis"), py::arg("angle"))
        .def("RotateX",
            (void(AbstractMesh_2_2::*)(double const)) &AbstractMesh_2_2::RotateX,
            " ", py::arg("theta"))
        .def("RotateY",
            (void(AbstractMesh_2_2::*)(double const)) &AbstractMesh_2_2::RotateY,
            " ", py::arg("theta"))
        .def("RotateZ",
            (void(AbstractMesh_2_2::*)(double const)) &AbstractMesh_2_2::RotateZ,
            " ", py::arg("theta"))
        .def("Rotate",
            (void(AbstractMesh_2_2::*)(double)) &AbstractMesh_2_2::Rotate,
            " ", py::arg("theta"))
        .def("RefreshMesh",
            (void(AbstractMesh_2_2::*)()) &AbstractMesh_2_2::RefreshMesh,
            " ")
        .def("IsMeshChanging",
            (bool(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::IsMeshChanging,
            " ")
        .def("CalculateMaximumContainingElementsPerProcess",
            (unsigned int(AbstractMesh_2_2::*)() const) &AbstractMesh_2_2::CalculateMaximumContainingElementsPerProcess,
            " ")
        .def("SetMeshHasChangedSinceLoading",
            (void(AbstractMesh_2_2::*)()) &AbstractMesh_2_2::SetMeshHasChangedSinceLoading,
            " ")
    ;
}
