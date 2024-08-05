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
#include "MutableVertexMesh.hpp"

#include "MutableVertexMesh2_2.cppwg.hpp"

namespace py = pybind11;
typedef MutableVertexMesh<2,2 > MutableVertexMesh2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class MutableVertexMesh2_2_Overrides : public MutableVertexMesh2_2{
    public:
    using MutableVertexMesh2_2::MutableVertexMesh;
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> point) override {
        PYBIND11_OVERRIDE(
            void,
            MutableVertexMesh2_2,
            SetNode,
                    nodeIndex,
        point);
    }
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableVertexMesh2_2,
            GetNumNodes,
            );
    }
    unsigned int GetNumElements() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableVertexMesh2_2,
            GetNumElements,
            );
    }
    void Clear() override {
        PYBIND11_OVERRIDE(
            void,
            MutableVertexMesh2_2,
            Clear,
            );
    }
    void ReMesh(::VertexElementMap & rElementMap) override {
        PYBIND11_OVERRIDE(
            void,
            MutableVertexMesh2_2,
            ReMesh,
                    rElementMap);
    }
    bool CheckForSwapsFromShortEdges() override {
        PYBIND11_OVERRIDE(
            bool,
            MutableVertexMesh2_2,
            CheckForSwapsFromShortEdges,
            );
    }
    void IdentifySwapType(::Node<2> * pNodeA, ::Node<2> * pNodeB) override {
        PYBIND11_OVERRIDE(
            void,
            MutableVertexMesh2_2,
            IdentifySwapType,
                    pNodeA,
        pNodeB);
    }

};
void register_MutableVertexMesh2_2_class(py::module &m){
py::class_<MutableVertexMesh2_2 , MutableVertexMesh2_2_Overrides , boost::shared_ptr<MutableVertexMesh2_2 >  , VertexMesh<2, 2>  >(m, "MutableVertexMesh2_2")
        .def(py::init<::std::vector<Node<2> *>, ::std::vector<VertexElement<2, 2> *>, double, double, double, double, double, double >(), py::arg("nodes"), py::arg("vertexElements"), py::arg("cellRearrangementThreshold") = 0.01, py::arg("t2Threshold") = 0.001, py::arg("cellRearrangementRatio") = 1.5, py::arg("protorosetteFormationProbability") = 0., py::arg("protorosetteResolutionProbabilityPerTimestep") = 0., py::arg("rosetteResolutionProbabilityPerTimestep") = 0.)
        .def(py::init< >())
        .def(
            "PerformNodeMerge",
            (void(MutableVertexMesh2_2::*)(::Node<2> *, ::Node<2> *)) &MutableVertexMesh2_2::PerformNodeMerge,
            " " , py::arg("pNodeA"), py::arg("pNodeB") )
        .def(
            "SetCellRearrangementThreshold",
            (void(MutableVertexMesh2_2::*)(double)) &MutableVertexMesh2_2::SetCellRearrangementThreshold,
            " " , py::arg("cellRearrangementThreshold") )
        .def(
            "SetT2Threshold",
            (void(MutableVertexMesh2_2::*)(double)) &MutableVertexMesh2_2::SetT2Threshold,
            " " , py::arg("t2Threshold") )
        .def(
            "SetCellRearrangementRatio",
            (void(MutableVertexMesh2_2::*)(double)) &MutableVertexMesh2_2::SetCellRearrangementRatio,
            " " , py::arg("cellRearrangementRatio") )
        .def(
            "SetProtorosetteFormationProbability",
            (void(MutableVertexMesh2_2::*)(double)) &MutableVertexMesh2_2::SetProtorosetteFormationProbability,
            " " , py::arg("protorosetteFormationProbability") )
        .def(
            "SetProtorosetteResolutionProbabilityPerTimestep",
            (void(MutableVertexMesh2_2::*)(double)) &MutableVertexMesh2_2::SetProtorosetteResolutionProbabilityPerTimestep,
            " " , py::arg("protorosetteResolutionProbabilityPerTimestep") )
        .def(
            "SetRosetteResolutionProbabilityPerTimestep",
            (void(MutableVertexMesh2_2::*)(double)) &MutableVertexMesh2_2::SetRosetteResolutionProbabilityPerTimestep,
            " " , py::arg("rosetteResolutionProbabilityPerTimestep") )
        .def(
            "SetNode",
            (void(MutableVertexMesh2_2::*)(unsigned int, ::ChastePoint<2>)) &MutableVertexMesh2_2::SetNode,
            " " , py::arg("nodeIndex"), py::arg("point") )
        .def(
            "SetCheckForInternalIntersections",
            (void(MutableVertexMesh2_2::*)(bool)) &MutableVertexMesh2_2::SetCheckForInternalIntersections,
            " " , py::arg("checkForInternalIntersections") )
        .def(
            "SetCheckForT3Swaps",
            (void(MutableVertexMesh2_2::*)(bool)) &MutableVertexMesh2_2::SetCheckForT3Swaps,
            " " , py::arg("checkForT3Swaps") )
        .def(
            "GetCellRearrangementThreshold",
            (double(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetCellRearrangementThreshold,
            " "  )
        .def(
            "GetT2Threshold",
            (double(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetT2Threshold,
            " "  )
        .def(
            "GetCellRearrangementRatio",
            (double(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetCellRearrangementRatio,
            " "  )
        .def(
            "GetProtorosetteFormationProbability",
            (double(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetProtorosetteFormationProbability,
            " "  )
        .def(
            "GetProtorosetteResolutionProbabilityPerTimestep",
            (double(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetProtorosetteResolutionProbabilityPerTimestep,
            " "  )
        .def(
            "GetRosetteResolutionProbabilityPerTimestep",
            (double(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetRosetteResolutionProbabilityPerTimestep,
            " "  )
        .def(
            "SetDistanceForT3SwapChecking",
            (void(MutableVertexMesh2_2::*)(double)) &MutableVertexMesh2_2::SetDistanceForT3SwapChecking,
            " " , py::arg("distanceForT3SwapChecking") )
        .def(
            "GetDistanceForT3SwapChecking",
            (double(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetDistanceForT3SwapChecking,
            " "  )
        .def(
            "GetNumNodes",
            (unsigned int(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetNumNodes,
            " "  )
        .def(
            "GetNumElements",
            (unsigned int(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetNumElements,
            " "  )
        .def(
            "GetCheckForInternalIntersections",
            (bool(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetCheckForInternalIntersections,
            " "  )
        .def(
            "GetCheckForT3Swaps",
            (bool(MutableVertexMesh2_2::*)() const ) &MutableVertexMesh2_2::GetCheckForT3Swaps,
            " "  )
        .def(
            "GetLocationsOfT1Swaps",
            (::std::vector<boost::numeric::ublas::c_vector<double, 2>>(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::GetLocationsOfT1Swaps,
            " "  )
        .def(
            "GetLastT2SwapLocation",
            (::boost::numeric::ublas::c_vector<double, 2>(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::GetLastT2SwapLocation,
            " "  )
        .def(
            "GetLocationsOfT3Swaps",
            (::std::vector<boost::numeric::ublas::c_vector<double, 2>>(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::GetLocationsOfT3Swaps,
            " "  )
        .def(
            "GetLocationsOfIntersectionSwaps",
            (::std::vector<boost::numeric::ublas::c_vector<double, 2>>(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::GetLocationsOfIntersectionSwaps,
            " "  )
        .def(
            "ClearLocationsOfT1Swaps",
            (void(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::ClearLocationsOfT1Swaps,
            " "  )
        .def(
            "ClearLocationsOfT3Swaps",
            (void(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::ClearLocationsOfT3Swaps,
            " "  )
        .def(
            "ClearLocationsOfIntersectionSwaps",
            (void(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::ClearLocationsOfIntersectionSwaps,
            " "  )
        .def(
            "AddNode",
            (unsigned int(MutableVertexMesh2_2::*)(::Node<2> *)) &MutableVertexMesh2_2::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "DeleteElementPriorToReMesh",
            (void(MutableVertexMesh2_2::*)(unsigned int)) &MutableVertexMesh2_2::DeleteElementPriorToReMesh,
            " " , py::arg("index") )
        .def(
            "DeleteNodePriorToReMesh",
            (void(MutableVertexMesh2_2::*)(unsigned int)) &MutableVertexMesh2_2::DeleteNodePriorToReMesh,
            " " , py::arg("index") )
        .def(
            "DivideElementAlongShortAxis",
            (unsigned int(MutableVertexMesh2_2::*)(::VertexElement<2, 2> *, bool)) &MutableVertexMesh2_2::DivideElementAlongShortAxis,
            " " , py::arg("pElement"), py::arg("placeOriginalElementBelow") = false )
        .def(
            "DivideElementAlongGivenAxis",
            (unsigned int(MutableVertexMesh2_2::*)(::VertexElement<2, 2> *, ::boost::numeric::ublas::c_vector<double, 2>, bool)) &MutableVertexMesh2_2::DivideElementAlongGivenAxis,
            " " , py::arg("pElement"), py::arg("axisOfDivision"), py::arg("placeOriginalElementBelow") = false )
        .def(
            "AddElement",
            (unsigned int(MutableVertexMesh2_2::*)(::VertexElement<2, 2> *)) &MutableVertexMesh2_2::AddElement,
            " " , py::arg("pNewElement") )
        .def(
            "CheckForT2Swaps",
            (bool(MutableVertexMesh2_2::*)(::VertexElementMap &)) &MutableVertexMesh2_2::CheckForT2Swaps,
            " " , py::arg("rElementMap") )
        .def(
            "Clear",
            (void(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::Clear,
            " "  )
        .def(
            "DivideEdge",
            (void(MutableVertexMesh2_2::*)(::Node<2> *, ::Node<2> *)) &MutableVertexMesh2_2::DivideEdge,
            " " , py::arg("pNodeA"), py::arg("pNodeB") )
        .def(
            "RemoveDeletedNodesAndElements",
            (void(MutableVertexMesh2_2::*)(::VertexElementMap &)) &MutableVertexMesh2_2::RemoveDeletedNodesAndElements,
            " " , py::arg("rElementMap") )
        .def(
            "RemoveDeletedNodes",
            (void(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::RemoveDeletedNodes,
            " "  )
        .def(
            "ReMesh",
            (void(MutableVertexMesh2_2::*)(::VertexElementMap &)) &MutableVertexMesh2_2::ReMesh,
            " " , py::arg("rElementMap") )
        .def(
            "ReMesh",
            (void(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::ReMesh,
            " "  )
        .def(
            "SetMeshOperationTracking",
            (void(MutableVertexMesh2_2::*)(bool const)) &MutableVertexMesh2_2::SetMeshOperationTracking,
            " " , py::arg("track") )
        .def(
            "GetOperationRecorder",
            (::VertexMeshOperationRecorder<2, 2> *(MutableVertexMesh2_2::*)()) &MutableVertexMesh2_2::GetOperationRecorder,
            " "  , py::return_value_policy::reference)
    ;
}
