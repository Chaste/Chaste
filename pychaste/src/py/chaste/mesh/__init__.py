"""Mesh Module"""

__copyright__ = """Copyright (c) 2005-2024, University of Oxford.
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

from chaste._pychaste_lib import (
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
from chaste._syntax import TemplatedClass

ChasteCuboid = TemplatedClass(
    {
        ("2",): ChasteCuboid_2,
        ("3",): ChasteCuboid_3,
    }
)

ChasteEllipsoid = TemplatedClass(
    {
        ("2",): ChasteEllipsoid_2,
        ("3",): ChasteEllipsoid_3,
    }
)

ChastePoint = TemplatedClass(
    {
        ("2",): ChastePoint_2,
        ("3",): ChastePoint_3,
    }
)

Edge = TemplatedClass(
    {
        ("2",): Edge_2,
        ("3",): Edge_3,
    }
)

EdgeHelper = TemplatedClass(
    {
        ("2",): EdgeHelper_2,
        ("3",): EdgeHelper_3,
    }
)

Element = TemplatedClass(
    {
        ("2",): Element_2_2,
        ("2", "2"): Element_2_2,
        ("3",): Element_3_3,
        ("3", "3"): Element_3_3,
    }
)

FluidSource = TemplatedClass(
    {
        ("2",): FluidSource_2,
        ("3",): FluidSource_3,
    }
)

ImmersedBoundaryElement = TemplatedClass(
    {
        ("1", "2"): ImmersedBoundaryElement_1_2,
        ("2",): ImmersedBoundaryElement_2_2,
        ("2", "2"): ImmersedBoundaryElement_2_2,
        ("2", "3"): ImmersedBoundaryElement_2_3,
        ("3",): ImmersedBoundaryElement_3_3,
        ("3", "3"): ImmersedBoundaryElement_3_3,
    }
)

ImmersedBoundaryMesh = TemplatedClass(
    {
        ("2",): ImmersedBoundaryMesh_2_2,
        ("2", "2"): ImmersedBoundaryMesh_2_2,
        ("3",): ImmersedBoundaryMesh_3_3,
        ("3", "3"): ImmersedBoundaryMesh_3_3,
    }
)

MutableElement = TemplatedClass(
    {
        ("1", "2"): MutableElement_1_2,
        ("2",): MutableElement_2_2,
        ("2", "2"): MutableElement_2_2,
        ("2", "3"): MutableElement_2_3,
        ("3",): MutableElement_3_3,
        ("3", "3"): MutableElement_3_3,
    }
)

MutableMesh = TemplatedClass(
    {
        ("2",): MutableMesh_2_2,
        ("2", "2"): MutableMesh_2_2,
        ("3",): MutableMesh_3_3,
        ("3", "3"): MutableMesh_3_3,
    }
)

MutableVertexMesh = TemplatedClass(
    {
        ("2",): MutableVertexMesh_2_2,
        ("2", "2"): MutableVertexMesh_2_2,
        ("3",): MutableVertexMesh_3_3,
        ("3", "3"): MutableVertexMesh_3_3,
    }
)

Node = TemplatedClass(
    {
        ("2",): Node_2,
        ("3",): Node_3,
    }
)

NodeAttributes = TemplatedClass(
    {
        ("2",): NodeAttributes_2,
        ("3",): NodeAttributes_3,
    }
)

NodesOnlyMesh = TemplatedClass(
    {
        ("2",): NodesOnlyMesh_2,
        ("3",): NodesOnlyMesh_3,
    }
)

PeriodicNodesOnlyMesh = TemplatedClass(
    {
        ("2",): PeriodicNodesOnlyMesh_2,
        ("3",): PeriodicNodesOnlyMesh_3,
    }
)

PottsElement = TemplatedClass(
    {
        ("2",): PottsElement_2,
        ("3",): PottsElement_3,
    }
)

PottsMesh = TemplatedClass(
    {
        ("2",): PottsMesh_2,
        ("3",): PottsMesh_3,
    }
)

PottsMeshGenerator = TemplatedClass(
    {
        ("2",): PottsMeshGenerator_2,
        ("3",): PottsMeshGenerator_3,
    }
)

PottsMeshWriter = TemplatedClass(
    {
        ("2",): PottsMeshWriter_2,
        ("3",): PottsMeshWriter_3,
    }
)

TetrahedralMesh = TemplatedClass(
    {
        ("2",): TetrahedralMesh_2_2,
        ("2", "2"): TetrahedralMesh_2_2,
        ("3",): TetrahedralMesh_3_3,
        ("3", "3"): TetrahedralMesh_3_3,
    }
)

VertexMesh = TemplatedClass(
    {
        ("2",): VertexMesh_2_2,
        ("2", "2"): VertexMesh_2_2,
        ("3",): VertexMesh_3_3,
        ("3", "3"): VertexMesh_3_3,
    }
)
