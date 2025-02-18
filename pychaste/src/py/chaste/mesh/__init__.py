"""Mesh Module"""

__copyright__ = """Copyright (c) 2005-2025, University of Oxford.
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
"""

from chaste._pychaste_all import (
    ChasteCuboid_2,
    ChasteCuboid_3,
    ChasteEllipsoid_2,
    ChasteEllipsoid_3,
    ChastePoint_2,
    ChastePoint_3,
    Cylindrical2dMesh,
    Cylindrical2dNodesOnlyMesh,
    Cylindrical2dVertexMesh,
    CylindricalHoneycombMeshGenerator,
    CylindricalHoneycombVertexMeshGenerator,
    Edge_2,
    Edge_3,
    EdgeHelper_2,
    EdgeHelper_3,
    EdgeOperation,
    Element_2_2,
    Element_3_3,
    FluidSource_2,
    FluidSource_3,
    HoneycombMeshGenerator,
    HoneycombVertexMeshGenerator,
    ImmersedBoundaryElement_1_2,
    ImmersedBoundaryElement_2_2,
    ImmersedBoundaryElement_2_3,
    ImmersedBoundaryElement_3_3,
    ImmersedBoundaryHoneycombMeshGenerator,
    ImmersedBoundaryMesh_2_2,
    ImmersedBoundaryMesh_3_3,
    ImmersedBoundaryPalisadeMeshGenerator,
    MutableElement_1_2,
    MutableElement_2_2,
    MutableElement_2_3,
    MutableElement_3_3,
    MutableMesh_2_2,
    MutableMesh_3_3,
    MutableVertexMesh_2_2,
    MutableVertexMesh_3_3,
    Node_2,
    Node_3,
    NodeAttributes_2,
    NodeAttributes_3,
    NodesOnlyMesh_2,
    NodesOnlyMesh_3,
    PeriodicNodesOnlyMesh_2,
    PeriodicNodesOnlyMesh_3,
    PottsElement_2,
    PottsElement_3,
    PottsMesh_2,
    PottsMesh_3,
    PottsMeshGenerator_2,
    PottsMeshGenerator_3,
    PottsMeshWriter_2,
    PottsMeshWriter_3,
    TetrahedralMesh_2_2,
    TetrahedralMesh_3_3,
    Toroidal2dMesh,
    Toroidal2dVertexMesh,
    ToroidalHoneycombMeshGenerator,
    ToroidalHoneycombVertexMeshGenerator,
    VertexMesh_2_2,
    VertexMesh_3_3,
    VoronoiVertexMeshGenerator,
)
from chaste._syntax import DeprecatedClass, TemplateClassDict

# Template Class Syntax
ChasteCuboid = TemplateClassDict(
    {
        ("2",): ChasteCuboid_2,
        ("3",): ChasteCuboid_3,
    }
)

ChasteEllipsoid = TemplateClassDict(
    {
        ("2",): ChasteEllipsoid_2,
        ("3",): ChasteEllipsoid_3,
    }
)

ChastePoint = TemplateClassDict(
    {
        ("2",): ChastePoint_2,
        ("3",): ChastePoint_3,
    }
)

Edge = TemplateClassDict(
    {
        ("2",): Edge_2,
        ("3",): Edge_3,
    }
)

EdgeHelper = TemplateClassDict(
    {
        ("2",): EdgeHelper_2,
        ("3",): EdgeHelper_3,
    }
)

Element = TemplateClassDict(
    {
        ("2",): Element_2_2,
        ("2", "2"): Element_2_2,
        ("3",): Element_3_3,
        ("3", "3"): Element_3_3,
    }
)

FluidSource = TemplateClassDict(
    {
        ("2",): FluidSource_2,
        ("3",): FluidSource_3,
    }
)

ImmersedBoundaryElement = TemplateClassDict(
    {
        ("1", "2"): ImmersedBoundaryElement_1_2,
        ("2",): ImmersedBoundaryElement_2_2,
        ("2", "2"): ImmersedBoundaryElement_2_2,
        ("2", "3"): ImmersedBoundaryElement_2_3,
        ("3",): ImmersedBoundaryElement_3_3,
        ("3", "3"): ImmersedBoundaryElement_3_3,
    }
)

ImmersedBoundaryMesh = TemplateClassDict(
    {
        ("2",): ImmersedBoundaryMesh_2_2,
        ("2", "2"): ImmersedBoundaryMesh_2_2,
        ("3",): ImmersedBoundaryMesh_3_3,
        ("3", "3"): ImmersedBoundaryMesh_3_3,
    }
)

MutableElement = TemplateClassDict(
    {
        ("1", "2"): MutableElement_1_2,
        ("2",): MutableElement_2_2,
        ("2", "2"): MutableElement_2_2,
        ("2", "3"): MutableElement_2_3,
        ("3",): MutableElement_3_3,
        ("3", "3"): MutableElement_3_3,
    }
)

MutableMesh = TemplateClassDict(
    {
        ("2",): MutableMesh_2_2,
        ("2", "2"): MutableMesh_2_2,
        ("3",): MutableMesh_3_3,
        ("3", "3"): MutableMesh_3_3,
    }
)

MutableVertexMesh = TemplateClassDict(
    {
        ("2",): MutableVertexMesh_2_2,
        ("2", "2"): MutableVertexMesh_2_2,
        ("3",): MutableVertexMesh_3_3,
        ("3", "3"): MutableVertexMesh_3_3,
    }
)

Node = TemplateClassDict(
    {
        ("2",): Node_2,
        ("3",): Node_3,
    }
)

NodeAttributes = TemplateClassDict(
    {
        ("2",): NodeAttributes_2,
        ("3",): NodeAttributes_3,
    }
)

NodesOnlyMesh = TemplateClassDict(
    {
        ("2",): NodesOnlyMesh_2,
        ("3",): NodesOnlyMesh_3,
    }
)

PeriodicNodesOnlyMesh = TemplateClassDict(
    {
        ("2",): PeriodicNodesOnlyMesh_2,
        ("3",): PeriodicNodesOnlyMesh_3,
    }
)

PottsElement = TemplateClassDict(
    {
        ("2",): PottsElement_2,
        ("3",): PottsElement_3,
    }
)

PottsMesh = TemplateClassDict(
    {
        ("2",): PottsMesh_2,
        ("3",): PottsMesh_3,
    }
)

PottsMeshGenerator = TemplateClassDict(
    {
        ("2",): PottsMeshGenerator_2,
        ("3",): PottsMeshGenerator_3,
    }
)

PottsMeshWriter = TemplateClassDict(
    {
        ("2",): PottsMeshWriter_2,
        ("3",): PottsMeshWriter_3,
    }
)

TetrahedralMesh = TemplateClassDict(
    {
        ("2",): TetrahedralMesh_2_2,
        ("2", "2"): TetrahedralMesh_2_2,
        ("3",): TetrahedralMesh_3_3,
        ("3", "3"): TetrahedralMesh_3_3,
    }
)

VertexMesh = TemplateClassDict(
    {
        ("2",): VertexMesh_2_2,
        ("2", "2"): VertexMesh_2_2,
        ("3",): VertexMesh_3_3,
        ("3", "3"): VertexMesh_3_3,
    }
)

# Deprecated Class Syntax
ChasteCuboid2 = DeprecatedClass("ChasteCuboid2", ChasteCuboid_2)
ChasteCuboid3 = DeprecatedClass("ChasteCuboid3", ChasteCuboid_3)
ChasteEllipsoid2 = DeprecatedClass("ChasteEllipsoid2", ChasteEllipsoid_2)
ChasteEllipsoid3 = DeprecatedClass("ChasteEllipsoid3", ChasteEllipsoid_3)
ChastePoint2 = DeprecatedClass("ChastePoint2", ChastePoint_2)
ChastePoint3 = DeprecatedClass("ChastePoint3", ChastePoint_3)
Edge2 = DeprecatedClass("Edge2", Edge_2)
Edge3 = DeprecatedClass("Edge3", Edge_3)
EdgeHelper2 = DeprecatedClass("EdgeHelper2", EdgeHelper_2)
EdgeHelper3 = DeprecatedClass("EdgeHelper3", EdgeHelper_3)
Element2_2 = DeprecatedClass("Element2_2", Element_2_2)
Element3_3 = DeprecatedClass("Element3_3", Element_3_3)
FluidSource2 = DeprecatedClass("FluidSource2", FluidSource_2)
FluidSource3 = DeprecatedClass("FluidSource3", FluidSource_3)
ImmersedBoundaryElement1_2 = DeprecatedClass("ImmersedBoundaryElement1_2", ImmersedBoundaryElement_1_2)
ImmersedBoundaryElement2_2 = DeprecatedClass("ImmersedBoundaryElement2_2", ImmersedBoundaryElement_2_2)
ImmersedBoundaryElement2_3 = DeprecatedClass("ImmersedBoundaryElement2_3", ImmersedBoundaryElement_2_3)
ImmersedBoundaryElement3_3 = DeprecatedClass("ImmersedBoundaryElement3_3", ImmersedBoundaryElement_3_3)
ImmersedBoundaryMesh2_2 = DeprecatedClass("ImmersedBoundaryMesh2_2", ImmersedBoundaryMesh_2_2)
ImmersedBoundaryMesh3_3 = DeprecatedClass("ImmersedBoundaryMesh3_3", ImmersedBoundaryMesh_3_3)
MutableElement1_2 = DeprecatedClass("MutableElement1_2", MutableElement_1_2)
MutableElement2_2 = DeprecatedClass("MutableElement2_2", MutableElement_2_2)
MutableElement2_3 = DeprecatedClass("MutableElement2_3", MutableElement_2_3)
MutableElement3_3 = DeprecatedClass("MutableElement3_3", MutableElement_3_3)
MutableMesh2_2 = DeprecatedClass("MutableMesh2_2", MutableMesh_2_2)
MutableMesh3_3 = DeprecatedClass("MutableMesh3_3", MutableMesh_3_3)
MutableVertexMesh2_2 = DeprecatedClass("MutableVertexMesh2_2", MutableVertexMesh_2_2)
MutableVertexMesh3_3 = DeprecatedClass("MutableVertexMesh3_3", MutableVertexMesh_3_3)
Node2 = DeprecatedClass("Node2", Node_2)
Node3 = DeprecatedClass("Node3", Node_3)
NodeAttributes2 = DeprecatedClass("NodeAttributes2", NodeAttributes_2)
NodeAttributes3 = DeprecatedClass("NodeAttributes3", NodeAttributes_3)
NodesOnlyMesh2 = DeprecatedClass("NodesOnlyMesh2", NodesOnlyMesh_2)
NodesOnlyMesh3 = DeprecatedClass("NodesOnlyMesh3", NodesOnlyMesh_3)
PeriodicNodesOnlyMesh2 = DeprecatedClass("PeriodicNodesOnlyMesh2", PeriodicNodesOnlyMesh_2)
PeriodicNodesOnlyMesh3 = DeprecatedClass("PeriodicNodesOnlyMesh3", PeriodicNodesOnlyMesh_3)
PottsElement2 = DeprecatedClass("PottsElement2", PottsElement_2)
PottsElement3 = DeprecatedClass("PottsElement3", PottsElement_3)
PottsMesh2 = DeprecatedClass("PottsMesh2", PottsMesh_2)
PottsMesh3 = DeprecatedClass("PottsMesh3", PottsMesh_3)
PottsMeshGenerator2 = DeprecatedClass("PottsMeshGenerator2", PottsMeshGenerator_2)
PottsMeshGenerator3 = DeprecatedClass("PottsMeshGenerator3", PottsMeshGenerator_3)
PottsMeshWriter2 = DeprecatedClass("PottsMeshWriter2", PottsMeshWriter_2)
PottsMeshWriter3 = DeprecatedClass("PottsMeshWriter3", PottsMeshWriter_3)
TetrahedralMesh2_2 = DeprecatedClass("TetrahedralMesh2_2", TetrahedralMesh_2_2)
TetrahedralMesh3_3 = DeprecatedClass("TetrahedralMesh3_3", TetrahedralMesh_3_3)
VertexMesh2_2 = DeprecatedClass("VertexMesh2_2", VertexMesh_2_2)
VertexMesh3_3 = DeprecatedClass("VertexMesh3_3", VertexMesh_3_3)