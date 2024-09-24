"""PDE Module"""

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
from chaste._syntax import TemplatedClass

AveragedSourceEllipticPde = TemplatedClass(
    {
        ("2",): AveragedSourceEllipticPde_2,
        ("3",): AveragedSourceEllipticPde_3,
    }
)

AveragedSourceParabolicPde = TemplatedClass(
    {
        ("2",): AveragedSourceParabolicPde_2,
        ("3",): AveragedSourceParabolicPde_3,
    }
)

CellBasedEllipticPdeSolver = TemplatedClass(
    {
        ("2",): CellBasedEllipticPdeSolver_2,
        ("3",): CellBasedEllipticPdeSolver_3,
    }
)

CellBasedParabolicPdeSolver = TemplatedClass(
    {
        ("2",): CellBasedParabolicPdeSolver_2,
        ("3",): CellBasedParabolicPdeSolver_3,
    }
)

CellwiseSourceEllipticPde = TemplatedClass(
    {
        ("2",): CellwiseSourceEllipticPde_2,
        ("3",): CellwiseSourceEllipticPde_3,
    }
)

CellwiseSourceParabolicPde = TemplatedClass(
    {
        ("2",): CellwiseSourceParabolicPde_2,
        ("3",): CellwiseSourceParabolicPde_3,
    }
)

ConstBoundaryCondition = TemplatedClass(
    {
        ("2",): ConstBoundaryCondition_2,
        ("3",): ConstBoundaryCondition_3,
    }
)

EllipticBoxDomainPdeModifier = TemplatedClass(
    {
        ("2",): EllipticBoxDomainPdeModifier_2,
        ("3",): EllipticBoxDomainPdeModifier_3,
    }
)

EllipticGrowingDomainPdeModifier = TemplatedClass(
    {
        ("2",): EllipticGrowingDomainPdeModifier_2,
        ("3",): EllipticGrowingDomainPdeModifier_3,
    }
)

ParabolicBoxDomainPdeModifier = TemplatedClass(
    {
        ("2",): ParabolicBoxDomainPdeModifier_2,
        ("3",): ParabolicBoxDomainPdeModifier_3,
    }
)

ParabolicGrowingDomainPdeModifier = TemplatedClass(
    {
        ("2",): ParabolicGrowingDomainPdeModifier_2,
        ("3",): ParabolicGrowingDomainPdeModifier_3,
    }
)


UniformSourceEllipticPde = TemplatedClass(
    {
        ("2",): UniformSourceEllipticPde_2,
        ("3",): UniformSourceEllipticPde_3,
    }
)

UniformSourceParabolicPde = TemplatedClass(
    {
        ("2",): UniformSourceParabolicPde_2,
        ("3",): UniformSourceParabolicPde_3,
    }
)

VolumeDependentAveragedSourceEllipticPde = TemplatedClass(
    {
        ("2",): VolumeDependentAveragedSourceEllipticPde_2,
        ("3",): VolumeDependentAveragedSourceEllipticPde_3,
    }
)
