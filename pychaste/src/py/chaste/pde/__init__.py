"""PDE Module"""

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
    AveragedSourceEllipticPde_2,
    AveragedSourceEllipticPde_3,
    AveragedSourceParabolicPde_2,
    AveragedSourceParabolicPde_3,
    CellBasedEllipticPdeSolver_2,
    CellBasedEllipticPdeSolver_3,
    CellBasedParabolicPdeSolver_2,
    CellBasedParabolicPdeSolver_3,
    CellwiseSourceEllipticPde_2,
    CellwiseSourceEllipticPde_3,
    CellwiseSourceParabolicPde_2,
    CellwiseSourceParabolicPde_3,
    ConstBoundaryCondition_2,
    ConstBoundaryCondition_3,
    EllipticBoxDomainPdeModifier_2,
    EllipticBoxDomainPdeModifier_3,
    EllipticGrowingDomainPdeModifier_2,
    EllipticGrowingDomainPdeModifier_3,
    ParabolicBoxDomainPdeModifier_2,
    ParabolicBoxDomainPdeModifier_3,
    ParabolicGrowingDomainPdeModifier_2,
    ParabolicGrowingDomainPdeModifier_3,
    PdeSimulationTime,
    UniformSourceEllipticPde_2,
    UniformSourceEllipticPde_3,
    UniformSourceParabolicPde_2,
    UniformSourceParabolicPde_3,
    VolumeDependentAveragedSourceEllipticPde_2,
    VolumeDependentAveragedSourceEllipticPde_3,
)
from chaste._syntax import DeprecatedClass, TemplateClassDict

# Template Class Syntax
AveragedSourceEllipticPde = TemplateClassDict(
    {
        ("2",): AveragedSourceEllipticPde_2,
        ("3",): AveragedSourceEllipticPde_3,
    }
)

AveragedSourceParabolicPde = TemplateClassDict(
    {
        ("2",): AveragedSourceParabolicPde_2,
        ("3",): AveragedSourceParabolicPde_3,
    }
)

CellBasedEllipticPdeSolver = TemplateClassDict(
    {
        ("2",): CellBasedEllipticPdeSolver_2,
        ("3",): CellBasedEllipticPdeSolver_3,
    }
)

CellBasedParabolicPdeSolver = TemplateClassDict(
    {
        ("2",): CellBasedParabolicPdeSolver_2,
        ("3",): CellBasedParabolicPdeSolver_3,
    }
)

CellwiseSourceEllipticPde = TemplateClassDict(
    {
        ("2",): CellwiseSourceEllipticPde_2,
        ("3",): CellwiseSourceEllipticPde_3,
    }
)

CellwiseSourceParabolicPde = TemplateClassDict(
    {
        ("2",): CellwiseSourceParabolicPde_2,
        ("3",): CellwiseSourceParabolicPde_3,
    }
)

ConstBoundaryCondition = TemplateClassDict(
    {
        ("2",): ConstBoundaryCondition_2,
        ("3",): ConstBoundaryCondition_3,
    }
)

EllipticBoxDomainPdeModifier = TemplateClassDict(
    {
        ("2",): EllipticBoxDomainPdeModifier_2,
        ("3",): EllipticBoxDomainPdeModifier_3,
    }
)

EllipticGrowingDomainPdeModifier = TemplateClassDict(
    {
        ("2",): EllipticGrowingDomainPdeModifier_2,
        ("3",): EllipticGrowingDomainPdeModifier_3,
    }
)

ParabolicBoxDomainPdeModifier = TemplateClassDict(
    {
        ("2",): ParabolicBoxDomainPdeModifier_2,
        ("3",): ParabolicBoxDomainPdeModifier_3,
    }
)

ParabolicGrowingDomainPdeModifier = TemplateClassDict(
    {
        ("2",): ParabolicGrowingDomainPdeModifier_2,
        ("3",): ParabolicGrowingDomainPdeModifier_3,
    }
)


UniformSourceEllipticPde = TemplateClassDict(
    {
        ("2",): UniformSourceEllipticPde_2,
        ("3",): UniformSourceEllipticPde_3,
    }
)

UniformSourceParabolicPde = TemplateClassDict(
    {
        ("2",): UniformSourceParabolicPde_2,
        ("3",): UniformSourceParabolicPde_3,
    }
)

VolumeDependentAveragedSourceEllipticPde = TemplateClassDict(
    {
        ("2",): VolumeDependentAveragedSourceEllipticPde_2,
        ("3",): VolumeDependentAveragedSourceEllipticPde_3,
    }
)

# Deprecated Class Syntax
AveragedSourceEllipticPde2 = DeprecatedClass("AveragedSourceEllipticPde2", AveragedSourceEllipticPde_2)
AveragedSourceEllipticPde3 = DeprecatedClass("AveragedSourceEllipticPde3", AveragedSourceEllipticPde_3)
AveragedSourceParabolicPde2 = DeprecatedClass("AveragedSourceParabolicPde2", AveragedSourceParabolicPde_2)
AveragedSourceParabolicPde3 = DeprecatedClass("AveragedSourceParabolicPde3", AveragedSourceParabolicPde_3)
CellBasedEllipticPdeSolver2 = DeprecatedClass("CellBasedEllipticPdeSolver2", CellBasedEllipticPdeSolver_2)
CellBasedEllipticPdeSolver3 = DeprecatedClass("CellBasedEllipticPdeSolver3", CellBasedEllipticPdeSolver_3)
CellBasedParabolicPdeSolver2 = DeprecatedClass("CellBasedParabolicPdeSolver2", CellBasedParabolicPdeSolver_2)
CellBasedParabolicPdeSolver3 = DeprecatedClass("CellBasedParabolicPdeSolver3", CellBasedParabolicPdeSolver_3)
CellwiseSourceEllipticPde2 = DeprecatedClass("CellwiseSourceEllipticPde2", CellwiseSourceEllipticPde_2)
CellwiseSourceEllipticPde3 = DeprecatedClass("CellwiseSourceEllipticPde3", CellwiseSourceEllipticPde_3)
CellwiseSourceParabolicPde2 = DeprecatedClass("CellwiseSourceParabolicPde2", CellwiseSourceParabolicPde_2)
CellwiseSourceParabolicPde3 = DeprecatedClass("CellwiseSourceParabolicPde3", CellwiseSourceParabolicPde_3)
ConstBoundaryCondition2 = DeprecatedClass("ConstBoundaryCondition2", ConstBoundaryCondition_2)
ConstBoundaryCondition3 = DeprecatedClass("ConstBoundaryCondition3", ConstBoundaryCondition_3)
EllipticBoxDomainPdeModifier2 = DeprecatedClass("EllipticBoxDomainPdeModifier2", EllipticBoxDomainPdeModifier_2)
EllipticBoxDomainPdeModifier3 = DeprecatedClass("EllipticBoxDomainPdeModifier3", EllipticBoxDomainPdeModifier_3)
EllipticGrowingDomainPdeModifier2 = DeprecatedClass("EllipticGrowingDomainPdeModifier2", EllipticGrowingDomainPdeModifier_2)
EllipticGrowingDomainPdeModifier3 = DeprecatedClass("EllipticGrowingDomainPdeModifier3", EllipticGrowingDomainPdeModifier_3)
ParabolicBoxDomainPdeModifier2 = DeprecatedClass("ParabolicBoxDomainPdeModifier2", ParabolicBoxDomainPdeModifier_2)
ParabolicBoxDomainPdeModifier3 = DeprecatedClass("ParabolicBoxDomainPdeModifier3", ParabolicBoxDomainPdeModifier_3)
ParabolicGrowingDomainPdeModifier2 = DeprecatedClass("ParabolicGrowingDomainPdeModifier2", ParabolicGrowingDomainPdeModifier_2)
ParabolicGrowingDomainPdeModifier3 = DeprecatedClass("ParabolicGrowingDomainPdeModifier3", ParabolicGrowingDomainPdeModifier_3)
UniformSourceEllipticPde2 = DeprecatedClass("UniformSourceEllipticPde2", UniformSourceEllipticPde_2)
UniformSourceEllipticPde3 = DeprecatedClass("UniformSourceEllipticPde3", UniformSourceEllipticPde_3)
UniformSourceParabolicPde2 = DeprecatedClass("UniformSourceParabolicPde2", UniformSourceParabolicPde_2)
UniformSourceParabolicPde3 = DeprecatedClass("UniformSourceParabolicPde3", UniformSourceParabolicPde_3)
VolumeDependentAveragedSourceEllipticPde2 = DeprecatedClass("VolumeDependentAveragedSourceEllipticPde2", VolumeDependentAveragedSourceEllipticPde_2)
VolumeDependentAveragedSourceEllipticPde3 = DeprecatedClass("VolumeDependentAveragedSourceEllipticPde3", VolumeDependentAveragedSourceEllipticPde_3)
