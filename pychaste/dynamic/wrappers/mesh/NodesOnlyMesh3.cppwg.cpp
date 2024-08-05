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

#include "NodesOnlyMesh3.cppwg.hpp"

namespace py = pybind11;
typedef NodesOnlyMesh<3 > NodesOnlyMesh3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef ::Node<3> * _Node_lt_3_gt_Ptr;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class NodesOnlyMesh3_Overrides : public NodesOnlyMesh3{
    public:
    using NodesOnlyMesh3::NodesOnlyMesh;
    void Clear() override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh3,
            Clear,
            );
    }
    unsigned int SolveNodeMapping(unsigned int index) const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh3,
            SolveNodeMapping,
                    index);
    }
    ::Node<3> * GetNodeOrHaloNode(unsigned int index) const  override {
        PYBIND11_OVERRIDE(
            _Node_lt_3_gt_Ptr,
            NodesOnlyMesh3,
            GetNodeOrHaloNode,
                    index);
    }
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh3,
            GetNumNodes,
            );
    }
    unsigned int GetMaximumNodeIndex() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh3,
            GetMaximumNodeIndex,
            );
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            NodesOnlyMesh3,
            GetWidth,
                    rDimension);
    }
    void ReMesh(::NodeMap & rMap) override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh3,
            ReMesh,
                    rMap);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<3> point, bool concreteMove) override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh3,
            SetNode,
                    nodeIndex,
        point,
        concreteMove);
    }
    unsigned int AddNode(::Node<3> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodesOnlyMesh3,
            AddNode,
                    pNewNode);
    }
    void DeleteNode(unsigned int index) override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh3,
            DeleteNode,
                    index);
    }
    void ConstructFromMeshReader(::AbstractMeshReader<3, 3> & rMeshReader) override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh3,
            ConstructFromMeshReader,
                    rMeshReader);
    }
    void SetUpBoxCollection(double cutOffLength, ::boost::numeric::ublas::c_vector<double, 6> domainSize, int numLocalRows, ::boost::numeric::ublas::c_vector<bool, 3> isDimPeriodic) override {
        PYBIND11_OVERRIDE(
            void,
            NodesOnlyMesh3,
            SetUpBoxCollection,
                    cutOffLength,
        domainSize,
        numLocalRows,
        isDimPeriodic);
    }

};
void register_NodesOnlyMesh3_class(py::module &m){
py::class_<NodesOnlyMesh3 , NodesOnlyMesh3_Overrides , boost::shared_ptr<NodesOnlyMesh3 >  , MutableMesh<3, 3>  >(m, "NodesOnlyMesh3")
        .def(py::init< >())
        .def(
            "ConstructNodesWithoutMesh",
            (void(NodesOnlyMesh3::*)(::std::vector<Node<3> *> const &, double)) &NodesOnlyMesh3::ConstructNodesWithoutMesh,
            " " , py::arg("rNodes"), py::arg("maxInteractionDistance") )
        .def(
            "ConstructNodesWithoutMesh",
            (void(NodesOnlyMesh3::*)(::std::vector<boost::shared_ptr<Node<3>>> const &, double)) &NodesOnlyMesh3::ConstructNodesWithoutMesh,
            " " , py::arg("rNodes"), py::arg("maxInteractionDistance") )
        .def(
            "ConstructNodesWithoutMesh",
            (void(NodesOnlyMesh3::*)(::AbstractMesh<3, 3> const &, double)) &NodesOnlyMesh3::ConstructNodesWithoutMesh,
            " " , py::arg("rGeneratingMesh"), py::arg("maxInteractionDistance") )
        .def(
            "rGetInitiallyOwnedNodes",
            (::std::vector<bool> &(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::rGetInitiallyOwnedNodes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "Clear",
            (void(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::Clear,
            " "  )
        .def(
            "SolveNodeMapping",
            (unsigned int(NodesOnlyMesh3::*)(unsigned int) const ) &NodesOnlyMesh3::SolveNodeMapping,
            " " , py::arg("index") )
        .def(
            "GetNodeOrHaloNode",
            (::Node<3> *(NodesOnlyMesh3::*)(unsigned int) const ) &NodesOnlyMesh3::GetNodeOrHaloNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "IsOwned",
            (bool(NodesOnlyMesh3::*)(::boost::numeric::ublas::c_vector<double, 3> &)) &NodesOnlyMesh3::IsOwned,
            " " , py::arg("location") )
        .def(
            "GetNumNodes",
            (unsigned int(NodesOnlyMesh3::*)() const ) &NodesOnlyMesh3::GetNumNodes,
            " "  )
        .def(
            "GetMaximumNodeIndex",
            (unsigned int(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::GetMaximumNodeIndex,
            " "  )
        .def(
            "SetMaximumInteractionDistance",
            (void(NodesOnlyMesh3::*)(double)) &NodesOnlyMesh3::SetMaximumInteractionDistance,
            " " , py::arg("maxDistance") )
        .def(
            "GetMaximumInteractionDistance",
            (double(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::GetMaximumInteractionDistance,
            " "  )
        .def(
            "GetWidth",
            (double(NodesOnlyMesh3::*)(unsigned int const &) const ) &NodesOnlyMesh3::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "SetCalculateNodeNeighbours",
            (void(NodesOnlyMesh3::*)(bool)) &NodesOnlyMesh3::SetCalculateNodeNeighbours,
            " " , py::arg("calculateNodeNeighbours") )
        .def(
            "CalculateInteriorNodePairs",
            (void(NodesOnlyMesh3::*)(::std::vector<std::pair<Node<3> *, Node<3> *>> &)) &NodesOnlyMesh3::CalculateInteriorNodePairs,
            " " , py::arg("rNodePairs") )
        .def(
            "CalculateBoundaryNodePairs",
            (void(NodesOnlyMesh3::*)(::std::vector<std::pair<Node<3> *, Node<3> *>> &)) &NodesOnlyMesh3::CalculateBoundaryNodePairs,
            " " , py::arg("rNodePairs") )
        .def(
            "ReMesh",
            (void(NodesOnlyMesh3::*)(::NodeMap &)) &NodesOnlyMesh3::ReMesh,
            " " , py::arg("rMap") )
        .def(
            "SetInitialBoxCollection",
            (void(NodesOnlyMesh3::*)(::boost::numeric::ublas::c_vector<double, 6> const, double)) &NodesOnlyMesh3::SetInitialBoxCollection,
            " " , py::arg("domainSize"), py::arg("maxInteractionDistance") )
        .def(
            "UpdateBoxCollection",
            (void(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::UpdateBoxCollection,
            " "  )
        .def(
            "ResizeBoxCollection",
            (void(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::ResizeBoxCollection,
            " "  )
        .def(
            "GetIsPeriodicAcrossProcsFromBoxCollection",
            (bool(NodesOnlyMesh3::*)() const ) &NodesOnlyMesh3::GetIsPeriodicAcrossProcsFromBoxCollection,
            " "  )
        .def(
            "AddNodesToBoxes",
            (void(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::AddNodesToBoxes,
            " "  )
        .def(
            "AddHaloNodesToBoxes",
            (void(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::AddHaloNodesToBoxes,
            " "  )
        .def(
            "CalculateNodesOutsideLocalDomain",
            (void(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::CalculateNodesOutsideLocalDomain,
            " "  )
        .def(
            "rGetNodesToSendLeft",
            (::std::vector<unsigned int> &(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::rGetNodesToSendLeft,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetNodesToSendRight",
            (::std::vector<unsigned int> &(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::rGetNodesToSendRight,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetHaloNodesToSendRight",
            (::std::vector<unsigned int> &(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::rGetHaloNodesToSendRight,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetHaloNodesToSendLeft",
            (::std::vector<unsigned int> &(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::rGetHaloNodesToSendLeft,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "AddHaloNode",
            (void(NodesOnlyMesh3::*)(::boost::shared_ptr<Node<3>>)) &NodesOnlyMesh3::AddHaloNode,
            " " , py::arg("pNewNode") )
        .def(
            "ClearHaloNodes",
            (void(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::ClearHaloNodes,
            " "  )
        .def(
            "SetNode",
            (void(NodesOnlyMesh3::*)(unsigned int, ::ChastePoint<3>, bool)) &NodesOnlyMesh3::SetNode,
            " " , py::arg("nodeIndex"), py::arg("point"), py::arg("concreteMove") = false )
        .def(
            "AddNode",
            (unsigned int(NodesOnlyMesh3::*)(::Node<3> *)) &NodesOnlyMesh3::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "AddMovedNode",
            (void(NodesOnlyMesh3::*)(::boost::shared_ptr<Node<3>>)) &NodesOnlyMesh3::AddMovedNode,
            " " , py::arg("pMovedNode") )
        .def(
            "DeleteNode",
            (void(NodesOnlyMesh3::*)(unsigned int)) &NodesOnlyMesh3::DeleteNode,
            " " , py::arg("index") )
        .def(
            "DeleteMovedNode",
            (void(NodesOnlyMesh3::*)(unsigned int)) &NodesOnlyMesh3::DeleteMovedNode,
            " " , py::arg("index") )
        .def(
            "SetMinimumNodeDomainBoundarySeparation",
            (void(NodesOnlyMesh3::*)(double)) &NodesOnlyMesh3::SetMinimumNodeDomainBoundarySeparation,
            " " , py::arg("separation") )
        .def(
            "LoadBalanceMesh",
            (void(NodesOnlyMesh3::*)()) &NodesOnlyMesh3::LoadBalanceMesh,
            " "  )
        .def(
            "ConstructFromMeshReader",
            (void(NodesOnlyMesh3::*)(::AbstractMeshReader<3, 3> &)) &NodesOnlyMesh3::ConstructFromMeshReader,
            " " , py::arg("rMeshReader") )
        .def(
            "GetAllNodeIndices",
            (::std::vector<unsigned int>(NodesOnlyMesh3::*)() const ) &NodesOnlyMesh3::GetAllNodeIndices,
            " "  )
    ;
}
