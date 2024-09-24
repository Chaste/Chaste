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

// This file is auto-generated. Manual changes will be overwritten. For changes
// to persist, update the configuration in pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NodesOnlyMesh.hpp"

#include "NodesOnlyMesh_2.cppwg.hpp"

namespace py = pybind11;
typedef NodesOnlyMesh<2> NodesOnlyMesh_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef ::Node<2> * _Node_lt_2_gt_Ptr;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class NodesOnlyMesh_2_Overrides : public NodesOnlyMesh_2
{
public:
    using NodesOnlyMesh_2::NodesOnlyMesh;
    void Clear() override
    {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh_2,
            Clear,
            );
    }
    unsigned int SolveNodeMapping(unsigned int index) const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh_2,
            SolveNodeMapping,
            index);
    }
    ::Node<2> * GetNodeOrHaloNode(unsigned int index) const override
    {
        PYBIND11_OVERRIDE(
            _Node_lt_2_gt_Ptr,
            NodesOnlyMesh_2,
            GetNodeOrHaloNode,
            index);
    }
    unsigned int GetNumNodes() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh_2,
            GetNumNodes,
            );
    }
    unsigned int GetMaximumNodeIndex() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh_2,
            GetMaximumNodeIndex,
            );
    }
    double GetWidth(unsigned int const & rDimension) const override
    {
        PYBIND11_OVERRIDE(
            double,
            NodesOnlyMesh_2,
            GetWidth,
            rDimension);
    }
    void ReMesh(::NodeMap & rMap) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh_2,
            ReMesh,
            rMap);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> point, bool concreteMove) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh_2,
            SetNode,
            nodeIndex,
            point,
            concreteMove);
    }
    unsigned int AddNode(::Node<2> * pNewNode) override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh_2,
            AddNode,
            pNewNode);
    }
    void DeleteNode(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh_2,
            DeleteNode,
            index);
    }
    void ConstructFromMeshReader(::AbstractMeshReader<2, 2> & rMeshReader) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh_2,
            ConstructFromMeshReader,
            rMeshReader);
    }
    void SetUpBoxCollection(double cutOffLength, ::boost::numeric::ublas::c_vector<double, 4> domainSize, int numLocalRows, ::boost::numeric::ublas::c_vector<bool, 2> isDimPeriodic) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh_2,
            SetUpBoxCollection,
            cutOffLength,
            domainSize,
            numLocalRows,
            isDimPeriodic);
    }
};

void register_NodesOnlyMesh_2_class(py::module &m)
{
    py::class_<NodesOnlyMesh_2, NodesOnlyMesh_2_Overrides, boost::shared_ptr<NodesOnlyMesh_2>, MutableMesh<2, 2>>(m, "NodesOnlyMesh_2")
        .def(py::init<>())
        .def("ConstructNodesWithoutMesh",
            (void(NodesOnlyMesh_2::*)(::std::vector<Node<2> *> const &, double)) &NodesOnlyMesh_2::ConstructNodesWithoutMesh,
            " ", py::arg("rNodes"), py::arg("maxInteractionDistance"))
        .def("ConstructNodesWithoutMesh",
            (void(NodesOnlyMesh_2::*)(::std::vector<boost::shared_ptr<Node<2>>> const &, double)) &NodesOnlyMesh_2::ConstructNodesWithoutMesh,
            " ", py::arg("rNodes"), py::arg("maxInteractionDistance"))
        .def("ConstructNodesWithoutMesh",
            (void(NodesOnlyMesh_2::*)(::AbstractMesh<2, 2> const &, double)) &NodesOnlyMesh_2::ConstructNodesWithoutMesh,
            " ", py::arg("rGeneratingMesh"), py::arg("maxInteractionDistance"))
        .def("rGetInitiallyOwnedNodes",
            (::std::vector<bool> &(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::rGetInitiallyOwnedNodes,
            " ", py::return_value_policy::reference_internal)
        .def("Clear",
            (void(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::Clear,
            " ")
        .def("SolveNodeMapping",
            (unsigned int(NodesOnlyMesh_2::*)(unsigned int) const) &NodesOnlyMesh_2::SolveNodeMapping,
            " ", py::arg("index"))
        .def("GetNodeOrHaloNode",
            (::Node<2> *(NodesOnlyMesh_2::*)(unsigned int) const) &NodesOnlyMesh_2::GetNodeOrHaloNode,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("IsOwned",
            (bool(NodesOnlyMesh_2::*)(::boost::numeric::ublas::c_vector<double, 2> &)) &NodesOnlyMesh_2::IsOwned,
            " ", py::arg("location"))
        .def("GetNumNodes",
            (unsigned int(NodesOnlyMesh_2::*)() const) &NodesOnlyMesh_2::GetNumNodes,
            " ")
        .def("GetMaximumNodeIndex",
            (unsigned int(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::GetMaximumNodeIndex,
            " ")
        .def("SetMaximumInteractionDistance",
            (void(NodesOnlyMesh_2::*)(double)) &NodesOnlyMesh_2::SetMaximumInteractionDistance,
            " ", py::arg("maxDistance"))
        .def("GetMaximumInteractionDistance",
            (double(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::GetMaximumInteractionDistance,
            " ")
        .def("GetWidth",
            (double(NodesOnlyMesh_2::*)(unsigned int const &) const) &NodesOnlyMesh_2::GetWidth,
            " ", py::arg("rDimension"))
        .def("SetCalculateNodeNeighbours",
            (void(NodesOnlyMesh_2::*)(bool)) &NodesOnlyMesh_2::SetCalculateNodeNeighbours,
            " ", py::arg("calculateNodeNeighbours"))
        .def("CalculateInteriorNodePairs",
            (void(NodesOnlyMesh_2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &)) &NodesOnlyMesh_2::CalculateInteriorNodePairs,
            " ", py::arg("rNodePairs"))
        .def("CalculateBoundaryNodePairs",
            (void(NodesOnlyMesh_2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &)) &NodesOnlyMesh_2::CalculateBoundaryNodePairs,
            " ", py::arg("rNodePairs"))
        .def("ReMesh",
            (void(NodesOnlyMesh_2::*)(::NodeMap &)) &NodesOnlyMesh_2::ReMesh,
            " ", py::arg("rMap"))
        .def("SetInitialBoxCollection",
            (void(NodesOnlyMesh_2::*)(::boost::numeric::ublas::c_vector<double, 4> const, double)) &NodesOnlyMesh_2::SetInitialBoxCollection,
            " ", py::arg("domainSize"), py::arg("maxInteractionDistance"))
        .def("UpdateBoxCollection",
            (void(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::UpdateBoxCollection,
            " ")
        .def("ResizeBoxCollection",
            (void(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::ResizeBoxCollection,
            " ")
        .def("GetIsPeriodicAcrossProcsFromBoxCollection",
            (bool(NodesOnlyMesh_2::*)() const) &NodesOnlyMesh_2::GetIsPeriodicAcrossProcsFromBoxCollection,
            " ")
        .def("AddNodesToBoxes",
            (void(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::AddNodesToBoxes,
            " ")
        .def("AddHaloNodesToBoxes",
            (void(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::AddHaloNodesToBoxes,
            " ")
        .def("CalculateNodesOutsideLocalDomain",
            (void(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::CalculateNodesOutsideLocalDomain,
            " ")
        .def("rGetNodesToSendLeft",
            (::std::vector<unsigned int> &(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::rGetNodesToSendLeft,
            " ", py::return_value_policy::reference_internal)
        .def("rGetNodesToSendRight",
            (::std::vector<unsigned int> &(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::rGetNodesToSendRight,
            " ", py::return_value_policy::reference_internal)
        .def("rGetHaloNodesToSendRight",
            (::std::vector<unsigned int> &(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::rGetHaloNodesToSendRight,
            " ", py::return_value_policy::reference_internal)
        .def("rGetHaloNodesToSendLeft",
            (::std::vector<unsigned int> &(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::rGetHaloNodesToSendLeft,
            " ", py::return_value_policy::reference_internal)
        .def("AddHaloNode",
            (void(NodesOnlyMesh_2::*)(::boost::shared_ptr<Node<2>>)) &NodesOnlyMesh_2::AddHaloNode,
            " ", py::arg("pNewNode"))
        .def("ClearHaloNodes",
            (void(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::ClearHaloNodes,
            " ")
        .def("SetNode",
            (void(NodesOnlyMesh_2::*)(unsigned int, ::ChastePoint<2>, bool)) &NodesOnlyMesh_2::SetNode,
            " ", py::arg("nodeIndex"), py::arg("point"), py::arg("concreteMove") = false)
        .def("AddNode",
            (unsigned int(NodesOnlyMesh_2::*)(::Node<2> *)) &NodesOnlyMesh_2::AddNode,
            " ", py::arg("pNewNode"))
        .def("AddMovedNode",
            (void(NodesOnlyMesh_2::*)(::boost::shared_ptr<Node<2>>)) &NodesOnlyMesh_2::AddMovedNode,
            " ", py::arg("pMovedNode"))
        .def("DeleteNode",
            (void(NodesOnlyMesh_2::*)(unsigned int)) &NodesOnlyMesh_2::DeleteNode,
            " ", py::arg("index"))
        .def("DeleteMovedNode",
            (void(NodesOnlyMesh_2::*)(unsigned int)) &NodesOnlyMesh_2::DeleteMovedNode,
            " ", py::arg("index"))
        .def("SetMinimumNodeDomainBoundarySeparation",
            (void(NodesOnlyMesh_2::*)(double)) &NodesOnlyMesh_2::SetMinimumNodeDomainBoundarySeparation,
            " ", py::arg("separation"))
        .def("LoadBalanceMesh",
            (void(NodesOnlyMesh_2::*)()) &NodesOnlyMesh_2::LoadBalanceMesh,
            " ")
        .def("ConstructFromMeshReader",
            (void(NodesOnlyMesh_2::*)(::AbstractMeshReader<2, 2> &)) &NodesOnlyMesh_2::ConstructFromMeshReader,
            " ", py::arg("rMeshReader"))
        .def("GetAllNodeIndices",
            (::std::vector<unsigned int>(NodesOnlyMesh_2::*)() const) &NodesOnlyMesh_2::GetAllNodeIndices,
            " ")
    ;
}
