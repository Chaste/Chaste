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

#include "MutableVertexMesh3_3.cppwg.hpp"

namespace py = pybind11;
typedef MutableVertexMesh<3,3 > MutableVertexMesh3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class MutableVertexMesh3_3_Overrides : public MutableVertexMesh3_3{
    public:
    using MutableVertexMesh3_3::MutableVertexMesh;
    void SetNode(unsigned int nodeIndex, ::ChastePoint<3> point) override {
        PYBIND11_OVERRIDE(
            void,
            MutableVertexMesh3_3,
            SetNode,
                    nodeIndex,
        point);
    }
    unsigned int GetNumNodes() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableVertexMesh3_3,
            GetNumNodes,
            );
    }
    unsigned int GetNumElements() const  override {
        PYBIND11_OVERRIDE(
            unsignedint,
            MutableVertexMesh3_3,
            GetNumElements,
            );
    }
    void Clear() override {
        PYBIND11_OVERRIDE(
            void,
            MutableVertexMesh3_3,
            Clear,
            );
    }
    void ReMesh(::VertexElementMap & rElementMap) override {
        PYBIND11_OVERRIDE(
            void,
            MutableVertexMesh3_3,
            ReMesh,
                    rElementMap);
    }
    bool CheckForSwapsFromShortEdges() override {
        PYBIND11_OVERRIDE(
            bool,
            MutableVertexMesh3_3,
            CheckForSwapsFromShortEdges,
            );
    }
    void IdentifySwapType(::Node<3> * pNodeA, ::Node<3> * pNodeB) override {
        PYBIND11_OVERRIDE(
            void,
            MutableVertexMesh3_3,
            IdentifySwapType,
                    pNodeA,
        pNodeB);
    }

};
void register_MutableVertexMesh3_3_class(py::module &m){
py::class_<MutableVertexMesh3_3 , MutableVertexMesh3_3_Overrides , boost::shared_ptr<MutableVertexMesh3_3 >  , VertexMesh<3, 3>  >(m, "MutableVertexMesh3_3")
        .def(py::init<::std::vector<Node<3> *>, ::std::vector<VertexElement<3, 3> *>, double, double, double, double, double, double >(), py::arg("nodes"), py::arg("vertexElements"), py::arg("cellRearrangementThreshold") = 0.01, py::arg("t2Threshold") = 0.001, py::arg("cellRearrangementRatio") = 1.5, py::arg("protorosetteFormationProbability") = 0., py::arg("protorosetteResolutionProbabilityPerTimestep") = 0., py::arg("rosetteResolutionProbabilityPerTimestep") = 0.)
        .def(py::init< >())
        .def(
            "PerformNodeMerge",
            (void(MutableVertexMesh3_3::*)(::Node<3> *, ::Node<3> *)) &MutableVertexMesh3_3::PerformNodeMerge,
            " " , py::arg("pNodeA"), py::arg("pNodeB") )
        .def(
            "SetCellRearrangementThreshold",
            (void(MutableVertexMesh3_3::*)(double)) &MutableVertexMesh3_3::SetCellRearrangementThreshold,
            " " , py::arg("cellRearrangementThreshold") )
        .def(
            "SetT2Threshold",
            (void(MutableVertexMesh3_3::*)(double)) &MutableVertexMesh3_3::SetT2Threshold,
            " " , py::arg("t2Threshold") )
        .def(
            "SetCellRearrangementRatio",
            (void(MutableVertexMesh3_3::*)(double)) &MutableVertexMesh3_3::SetCellRearrangementRatio,
            " " , py::arg("cellRearrangementRatio") )
        .def(
            "SetProtorosetteFormationProbability",
            (void(MutableVertexMesh3_3::*)(double)) &MutableVertexMesh3_3::SetProtorosetteFormationProbability,
            " " , py::arg("protorosetteFormationProbability") )
        .def(
            "SetProtorosetteResolutionProbabilityPerTimestep",
            (void(MutableVertexMesh3_3::*)(double)) &MutableVertexMesh3_3::SetProtorosetteResolutionProbabilityPerTimestep,
            " " , py::arg("protorosetteResolutionProbabilityPerTimestep") )
        .def(
            "SetRosetteResolutionProbabilityPerTimestep",
            (void(MutableVertexMesh3_3::*)(double)) &MutableVertexMesh3_3::SetRosetteResolutionProbabilityPerTimestep,
            " " , py::arg("rosetteResolutionProbabilityPerTimestep") )
        .def(
            "SetNode",
            (void(MutableVertexMesh3_3::*)(unsigned int, ::ChastePoint<3>)) &MutableVertexMesh3_3::SetNode,
            " " , py::arg("nodeIndex"), py::arg("point") )
        .def(
            "SetCheckForInternalIntersections",
            (void(MutableVertexMesh3_3::*)(bool)) &MutableVertexMesh3_3::SetCheckForInternalIntersections,
            " " , py::arg("checkForInternalIntersections") )
        .def(
            "SetCheckForT3Swaps",
            (void(MutableVertexMesh3_3::*)(bool)) &MutableVertexMesh3_3::SetCheckForT3Swaps,
            " " , py::arg("checkForT3Swaps") )
        .def(
            "GetCellRearrangementThreshold",
            (double(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetCellRearrangementThreshold,
            " "  )
        .def(
            "GetT2Threshold",
            (double(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetT2Threshold,
            " "  )
        .def(
            "GetCellRearrangementRatio",
            (double(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetCellRearrangementRatio,
            " "  )
        .def(
            "GetProtorosetteFormationProbability",
            (double(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetProtorosetteFormationProbability,
            " "  )
        .def(
            "GetProtorosetteResolutionProbabilityPerTimestep",
            (double(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetProtorosetteResolutionProbabilityPerTimestep,
            " "  )
        .def(
            "GetRosetteResolutionProbabilityPerTimestep",
            (double(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetRosetteResolutionProbabilityPerTimestep,
            " "  )
        .def(
            "SetDistanceForT3SwapChecking",
            (void(MutableVertexMesh3_3::*)(double)) &MutableVertexMesh3_3::SetDistanceForT3SwapChecking,
            " " , py::arg("distanceForT3SwapChecking") )
        .def(
            "GetDistanceForT3SwapChecking",
            (double(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetDistanceForT3SwapChecking,
            " "  )
        .def(
            "GetNumNodes",
            (unsigned int(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetNumNodes,
            " "  )
        .def(
            "GetNumElements",
            (unsigned int(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetNumElements,
            " "  )
        .def(
            "GetCheckForInternalIntersections",
            (bool(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetCheckForInternalIntersections,
            " "  )
        .def(
            "GetCheckForT3Swaps",
            (bool(MutableVertexMesh3_3::*)() const ) &MutableVertexMesh3_3::GetCheckForT3Swaps,
            " "  )
        .def(
            "GetLocationsOfT1Swaps",
            (::std::vector<boost::numeric::ublas::c_vector<double, 3>>(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::GetLocationsOfT1Swaps,
            " "  )
        .def(
            "GetLastT2SwapLocation",
            (::boost::numeric::ublas::c_vector<double, 3>(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::GetLastT2SwapLocation,
            " "  )
        .def(
            "GetLocationsOfT3Swaps",
            (::std::vector<boost::numeric::ublas::c_vector<double, 3>>(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::GetLocationsOfT3Swaps,
            " "  )
        .def(
            "GetLocationsOfIntersectionSwaps",
            (::std::vector<boost::numeric::ublas::c_vector<double, 3>>(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::GetLocationsOfIntersectionSwaps,
            " "  )
        .def(
            "ClearLocationsOfT1Swaps",
            (void(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::ClearLocationsOfT1Swaps,
            " "  )
        .def(
            "ClearLocationsOfT3Swaps",
            (void(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::ClearLocationsOfT3Swaps,
            " "  )
        .def(
            "ClearLocationsOfIntersectionSwaps",
            (void(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::ClearLocationsOfIntersectionSwaps,
            " "  )
        .def(
            "AddNode",
            (unsigned int(MutableVertexMesh3_3::*)(::Node<3> *)) &MutableVertexMesh3_3::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "DeleteElementPriorToReMesh",
            (void(MutableVertexMesh3_3::*)(unsigned int)) &MutableVertexMesh3_3::DeleteElementPriorToReMesh,
            " " , py::arg("index") )
        .def(
            "DeleteNodePriorToReMesh",
            (void(MutableVertexMesh3_3::*)(unsigned int)) &MutableVertexMesh3_3::DeleteNodePriorToReMesh,
            " " , py::arg("index") )
        .def(
            "DivideElementAlongShortAxis",
            (unsigned int(MutableVertexMesh3_3::*)(::VertexElement<3, 3> *, bool)) &MutableVertexMesh3_3::DivideElementAlongShortAxis,
            " " , py::arg("pElement"), py::arg("placeOriginalElementBelow") = false )
        .def(
            "DivideElementAlongGivenAxis",
            (unsigned int(MutableVertexMesh3_3::*)(::VertexElement<3, 3> *, ::boost::numeric::ublas::c_vector<double, 3>, bool)) &MutableVertexMesh3_3::DivideElementAlongGivenAxis,
            " " , py::arg("pElement"), py::arg("axisOfDivision"), py::arg("placeOriginalElementBelow") = false )
        .def(
            "AddElement",
            (unsigned int(MutableVertexMesh3_3::*)(::VertexElement<3, 3> *)) &MutableVertexMesh3_3::AddElement,
            " " , py::arg("pNewElement") )
        .def(
            "CheckForT2Swaps",
            (bool(MutableVertexMesh3_3::*)(::VertexElementMap &)) &MutableVertexMesh3_3::CheckForT2Swaps,
            " " , py::arg("rElementMap") )
        .def(
            "Clear",
            (void(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::Clear,
            " "  )
        .def(
            "DivideEdge",
            (void(MutableVertexMesh3_3::*)(::Node<3> *, ::Node<3> *)) &MutableVertexMesh3_3::DivideEdge,
            " " , py::arg("pNodeA"), py::arg("pNodeB") )
        .def(
            "RemoveDeletedNodesAndElements",
            (void(MutableVertexMesh3_3::*)(::VertexElementMap &)) &MutableVertexMesh3_3::RemoveDeletedNodesAndElements,
            " " , py::arg("rElementMap") )
        .def(
            "RemoveDeletedNodes",
            (void(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::RemoveDeletedNodes,
            " "  )
        .def(
            "ReMesh",
            (void(MutableVertexMesh3_3::*)(::VertexElementMap &)) &MutableVertexMesh3_3::ReMesh,
            " " , py::arg("rElementMap") )
        .def(
            "ReMesh",
            (void(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::ReMesh,
            " "  )
        .def(
            "SetMeshOperationTracking",
            (void(MutableVertexMesh3_3::*)(bool const)) &MutableVertexMesh3_3::SetMeshOperationTracking,
            " " , py::arg("track") )
        .def(
            "GetOperationRecorder",
            (::VertexMeshOperationRecorder<3, 3> *(MutableVertexMesh3_3::*)()) &MutableVertexMesh3_3::GetOperationRecorder,
            " "  , py::return_value_policy::reference)
    ;
}
