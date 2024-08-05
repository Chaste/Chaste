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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NodesOnlyMesh.hpp"

#include "NodesOnlyMesh2.cppwg.hpp"

namespace py = pybind11;
typedef NodesOnlyMesh<2 > NodesOnlyMesh2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef ::Node<2> * _Node_lt_2_gt_Ptr;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class NodesOnlyMesh2_Overrides : public NodesOnlyMesh2{
    public:
    using NodesOnlyMesh2::NodesOnlyMesh;
    void Clear() override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh2,
            Clear,
            );
    }
    unsigned int SolveNodeMapping(unsigned int index) const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh2,
            SolveNodeMapping,
                    index);
    }
    ::Node<2> * GetNodeOrHaloNode(unsigned int index) const  override {
        PYBIND11_OVERRIDE(
            _Node_lt_2_gt_Ptr,
            NodesOnlyMesh2,
            GetNodeOrHaloNode,
                    index);
    }
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh2,
            GetNumNodes,
            );
    }
    unsigned int GetMaximumNodeIndex() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh2,
            GetMaximumNodeIndex,
            );
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            NodesOnlyMesh2,
            GetWidth,
                    rDimension);
    }
    void ReMesh(::NodeMap & rMap) override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh2,
            ReMesh,
                    rMap);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> point, bool concreteMove) override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh2,
            SetNode,
                    nodeIndex,
        point,
        concreteMove);
    }
    unsigned int AddNode(::Node<2> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh2,
            AddNode,
                    pNewNode);
    }
    void DeleteNode(unsigned int index) override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh2,
            DeleteNode,
                    index);
    }
    void ConstructFromMeshReader(::AbstractMeshReader<2, 2> & rMeshReader) override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh2,
            ConstructFromMeshReader,
                    rMeshReader);
    }
    void SetUpBoxCollection(double cutOffLength, ::boost::numeric::ublas::c_vector<double, 4> domainSize, int numLocalRows, ::boost::numeric::ublas::c_vector<bool, 2> isDimPeriodic) override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh2,
            SetUpBoxCollection,
                    cutOffLength,
        domainSize,
        numLocalRows,
        isDimPeriodic);
    }

};
void register_NodesOnlyMesh2_class(py::module &m){
py::class_<NodesOnlyMesh2 , NodesOnlyMesh2_Overrides , boost::shared_ptr<NodesOnlyMesh2 >  , MutableMesh<2, 2>  >(m, "NodesOnlyMesh2")
        .def(py::init< >())
        .def(
            "ConstructNodesWithoutMesh",
            (void(NodesOnlyMesh2::*)(::std::vector<Node<2> *> const &, double)) &NodesOnlyMesh2::ConstructNodesWithoutMesh,
            " " , py::arg("rNodes"), py::arg("maxInteractionDistance") )
        .def(
            "ConstructNodesWithoutMesh",
            (void(NodesOnlyMesh2::*)(::std::vector<boost::shared_ptr<Node<2>>> const &, double)) &NodesOnlyMesh2::ConstructNodesWithoutMesh,
            " " , py::arg("rNodes"), py::arg("maxInteractionDistance") )
        .def(
            "ConstructNodesWithoutMesh",
            (void(NodesOnlyMesh2::*)(::AbstractMesh<2, 2> const &, double)) &NodesOnlyMesh2::ConstructNodesWithoutMesh,
            " " , py::arg("rGeneratingMesh"), py::arg("maxInteractionDistance") )
        .def(
            "rGetInitiallyOwnedNodes",
            (::std::vector<bool> &(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::rGetInitiallyOwnedNodes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "Clear",
            (void(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::Clear,
            " "  )
        .def(
            "SolveNodeMapping",
            (unsigned int(NodesOnlyMesh2::*)(unsigned int) const ) &NodesOnlyMesh2::SolveNodeMapping,
            " " , py::arg("index") )
        .def(
            "GetNodeOrHaloNode",
            (::Node<2> *(NodesOnlyMesh2::*)(unsigned int) const ) &NodesOnlyMesh2::GetNodeOrHaloNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "IsOwned",
            (bool(NodesOnlyMesh2::*)(::boost::numeric::ublas::c_vector<double, 2> &)) &NodesOnlyMesh2::IsOwned,
            " " , py::arg("location") )
        .def(
            "GetNumNodes",
            (unsigned int(NodesOnlyMesh2::*)() const ) &NodesOnlyMesh2::GetNumNodes,
            " "  )
        .def(
            "GetMaximumNodeIndex",
            (unsigned int(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::GetMaximumNodeIndex,
            " "  )
        .def(
            "SetMaximumInteractionDistance",
            (void(NodesOnlyMesh2::*)(double)) &NodesOnlyMesh2::SetMaximumInteractionDistance,
            " " , py::arg("maxDistance") )
        .def(
            "GetMaximumInteractionDistance",
            (double(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::GetMaximumInteractionDistance,
            " "  )
        .def(
            "GetWidth",
            (double(NodesOnlyMesh2::*)(unsigned int const &) const ) &NodesOnlyMesh2::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "SetCalculateNodeNeighbours",
            (void(NodesOnlyMesh2::*)(bool)) &NodesOnlyMesh2::SetCalculateNodeNeighbours,
            " " , py::arg("calculateNodeNeighbours") )
        .def(
            "CalculateInteriorNodePairs",
            (void(NodesOnlyMesh2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &)) &NodesOnlyMesh2::CalculateInteriorNodePairs,
            " " , py::arg("rNodePairs") )
        .def(
            "CalculateBoundaryNodePairs",
            (void(NodesOnlyMesh2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &)) &NodesOnlyMesh2::CalculateBoundaryNodePairs,
            " " , py::arg("rNodePairs") )
        .def(
            "ReMesh",
            (void(NodesOnlyMesh2::*)(::NodeMap &)) &NodesOnlyMesh2::ReMesh,
            " " , py::arg("rMap") )
        .def(
            "SetInitialBoxCollection",
            (void(NodesOnlyMesh2::*)(::boost::numeric::ublas::c_vector<double, 4> const, double)) &NodesOnlyMesh2::SetInitialBoxCollection,
            " " , py::arg("domainSize"), py::arg("maxInteractionDistance") )
        .def(
            "UpdateBoxCollection",
            (void(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::UpdateBoxCollection,
            " "  )
        .def(
            "ResizeBoxCollection",
            (void(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::ResizeBoxCollection,
            " "  )
        .def(
            "GetIsPeriodicAcrossProcsFromBoxCollection",
            (bool(NodesOnlyMesh2::*)() const ) &NodesOnlyMesh2::GetIsPeriodicAcrossProcsFromBoxCollection,
            " "  )
        .def(
            "AddNodesToBoxes",
            (void(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::AddNodesToBoxes,
            " "  )
        .def(
            "AddHaloNodesToBoxes",
            (void(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::AddHaloNodesToBoxes,
            " "  )
        .def(
            "CalculateNodesOutsideLocalDomain",
            (void(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::CalculateNodesOutsideLocalDomain,
            " "  )
        .def(
            "rGetNodesToSendLeft",
            (::std::vector<unsigned int> &(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::rGetNodesToSendLeft,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetNodesToSendRight",
            (::std::vector<unsigned int> &(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::rGetNodesToSendRight,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetHaloNodesToSendRight",
            (::std::vector<unsigned int> &(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::rGetHaloNodesToSendRight,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetHaloNodesToSendLeft",
            (::std::vector<unsigned int> &(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::rGetHaloNodesToSendLeft,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "AddHaloNode",
            (void(NodesOnlyMesh2::*)(::boost::shared_ptr<Node<2>>)) &NodesOnlyMesh2::AddHaloNode,
            " " , py::arg("pNewNode") )
        .def(
            "ClearHaloNodes",
            (void(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::ClearHaloNodes,
            " "  )
        .def(
            "SetNode",
            (void(NodesOnlyMesh2::*)(unsigned int, ::ChastePoint<2>, bool)) &NodesOnlyMesh2::SetNode,
            " " , py::arg("nodeIndex"), py::arg("point"), py::arg("concreteMove") = false )
        .def(
            "AddNode",
            (unsigned int(NodesOnlyMesh2::*)(::Node<2> *)) &NodesOnlyMesh2::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "AddMovedNode",
            (void(NodesOnlyMesh2::*)(::boost::shared_ptr<Node<2>>)) &NodesOnlyMesh2::AddMovedNode,
            " " , py::arg("pMovedNode") )
        .def(
            "DeleteNode",
            (void(NodesOnlyMesh2::*)(unsigned int)) &NodesOnlyMesh2::DeleteNode,
            " " , py::arg("index") )
        .def(
            "DeleteMovedNode",
            (void(NodesOnlyMesh2::*)(unsigned int)) &NodesOnlyMesh2::DeleteMovedNode,
            " " , py::arg("index") )
        .def(
            "SetMinimumNodeDomainBoundarySeparation",
            (void(NodesOnlyMesh2::*)(double)) &NodesOnlyMesh2::SetMinimumNodeDomainBoundarySeparation,
            " " , py::arg("separation") )
        .def(
            "LoadBalanceMesh",
            (void(NodesOnlyMesh2::*)()) &NodesOnlyMesh2::LoadBalanceMesh,
            " "  )
        .def(
            "ConstructFromMeshReader",
            (void(NodesOnlyMesh2::*)(::AbstractMeshReader<2, 2> &)) &NodesOnlyMesh2::ConstructFromMeshReader,
            " " , py::arg("rMeshReader") )
        .def(
            "GetAllNodeIndices",
            (::std::vector<unsigned int>(NodesOnlyMesh2::*)() const ) &NodesOnlyMesh2::GetAllNodeIndices,
            " "  )
    ;
}
