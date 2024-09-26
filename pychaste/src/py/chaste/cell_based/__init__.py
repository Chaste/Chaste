"""Cell Based Module"""

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
    AdhesionPottsUpdateRule_2,
    AdhesionPottsUpdateRule_3,
    Alarcon2004OxygenBasedCellCycleModel,
    AlwaysDivideCellCycleModel,
    ApcOneHitCellMutationState,
    ApcTwoHitCellMutationState,
    ApoptoticCellKiller_2,
    ApoptoticCellKiller_3,
    ApoptoticCellProperty,
    AttractingPlaneBoundaryCondition_2_2,
    AttractingPlaneBoundaryCondition_3_3,
    BernoulliTrialCellCycleModel,
    BetaCateninOneHitCellMutationState,
    BiasedBernoulliTrialCellCycleModel,
    BoundaryNodeWriter_2_2,
    BoundaryNodeWriter_3_3,
    BuskeAdhesiveForce_2,
    BuskeAdhesiveForce_3,
    BuskeCompressionForce_2,
    BuskeCompressionForce_3,
    BuskeElasticForce_2,
    BuskeElasticForce_3,
    CaBasedCellPopulation_2,
    CaBasedCellPopulation_3,
    Cell,
    CellAgesWriter_2_2,
    CellAgesWriter_3_3,
    CellAncestor,
    CellAncestorWriter_2_2,
    CellAncestorWriter_3_3,
    CellAppliedForceWriter_2_2,
    CellAppliedForceWriter_3_3,
    CellCycleModelProteinConcentrationsWriter_2_2,
    CellCycleModelProteinConcentrationsWriter_3_3,
    CellData,
    CellDataItemWriter_2_2,
    CellDataItemWriter_3_3,
    CellDeltaNotchWriter_2_2,
    CellDeltaNotchWriter_3_3,
    CellDivisionLocationsWriter_2_2,
    CellDivisionLocationsWriter_3_3,
    CellEdgeData,
    CellId,
    CellIdWriter_2_2,
    CellIdWriter_3_3,
    CellLabel,
    CellLabelWriter_2_2,
    CellLabelWriter_3_3,
    CellLocationIndexWriter_2_2,
    CellLocationIndexWriter_3_3,
    CellMutationStatesCountWriter_2_2,
    CellMutationStatesCountWriter_3_3,
    CellMutationStatesWriter_2_2,
    CellMutationStatesWriter_3_3,
    CellPopulationAdjacencyMatrixWriter_2_2,
    CellPopulationAdjacencyMatrixWriter_3_3,
    CellPopulationAreaWriter_2_2,
    CellPopulationAreaWriter_3_3,
    CellPopulationElementWriter_2_2,
    CellPopulationElementWriter_3_3,
    CellProliferativePhasesCountWriter_2_2,
    CellProliferativePhasesCountWriter_3_3,
    CellProliferativePhasesWriter_2_2,
    CellProliferativePhasesWriter_3_3,
    CellProliferativeTypesCountWriter_2_2,
    CellProliferativeTypesCountWriter_3_3,
    CellProliferativeTypesWriter_2_2,
    CellProliferativeTypesWriter_3_3,
    CellPropertyCollection,
    CellPropertyRegistry,
    CellRadiusWriter_2_2,
    CellRadiusWriter_3_3,
    CellRemovalLocationsWriter_2_2,
    CellRemovalLocationsWriter_3_3,
    CellRosetteRankWriter_2_2,
    CellRosetteRankWriter_3_3,
    CellsGenerator_Alarcon2004OxygenBasedCellCycleModel_2,
    CellsGenerator_Alarcon2004OxygenBasedCellCycleModel_3,
    CellsGenerator_AlwaysDivideCellCycleModel_2,
    CellsGenerator_AlwaysDivideCellCycleModel_3,
    CellsGenerator_BernoulliTrialCellCycleModel_2,
    CellsGenerator_BernoulliTrialCellCycleModel_3,
    CellsGenerator_BiasedBernoulliTrialCellCycleModel_2,
    CellsGenerator_BiasedBernoulliTrialCellCycleModel_3,
    CellsGenerator_ContactInhibitionCellCycleModel_2,
    CellsGenerator_ContactInhibitionCellCycleModel_3,
    CellsGenerator_ExponentialG1GenerationalCellCycleModel_2,
    CellsGenerator_ExponentialG1GenerationalCellCycleModel_3,
    CellsGenerator_FixedG1GenerationalCellCycleModel_2,
    CellsGenerator_FixedG1GenerationalCellCycleModel_3,
    CellsGenerator_FixedSequenceCellCycleModel_2,
    CellsGenerator_FixedSequenceCellCycleModel_3,
    CellsGenerator_GammaG1CellCycleModel_2,
    CellsGenerator_GammaG1CellCycleModel_3,
    CellsGenerator_LabelDependentBernoulliTrialCellCycleModel_2,
    CellsGenerator_LabelDependentBernoulliTrialCellCycleModel_3,
    CellsGenerator_NoCellCycleModel_2,
    CellsGenerator_NoCellCycleModel_3,
    CellsGenerator_SimpleOxygenBasedCellCycleModel_2,
    CellsGenerator_SimpleOxygenBasedCellCycleModel_3,
    CellsGenerator_StochasticOxygenBasedCellCycleModel_2,
    CellsGenerator_StochasticOxygenBasedCellCycleModel_3,
    CellsGenerator_TysonNovakCellCycleModel_2,
    CellsGenerator_TysonNovakCellCycleModel_3,
    CellsGenerator_UniformCellCycleModel_2,
    CellsGenerator_UniformCellCycleModel_3,
    CellsGenerator_UniformG1GenerationalCellCycleModel_2,
    CellsGenerator_UniformG1GenerationalCellCycleModel_3,
    CellSrnModel,
    CellVolumesWriter_2_2,
    CellVolumesWriter_3_3,
    ChemotacticForce_2,
    ChemotacticForce_3,
    ChemotaxisPottsUpdateRule_2,
    ChemotaxisPottsUpdateRule_3,
    ContactInhibitionCellCycleModel,
    DefaultCellProliferativeType,
    DeltaNotchEdgeInteriorTrackingModifier_2,
    DeltaNotchEdgeInteriorTrackingModifier_3,
    DeltaNotchEdgeSrnModel,
    DeltaNotchEdgeTrackingModifier_2,
    DeltaNotchEdgeTrackingModifier_3,
    DeltaNotchInteriorSrnModel,
    DeltaNotchSrnModel,
    DeltaNotchTrackingModifier_2,
    DeltaNotchTrackingModifier_3,
    DifferentialAdhesionGeneralisedLinearSpringForce_2_2,
    DifferentialAdhesionGeneralisedLinearSpringForce_3_3,
    DifferentialAdhesionPottsUpdateRule_2,
    DifferentialAdhesionPottsUpdateRule_3,
    DifferentiatedCellProliferativeType,
    DiffusionCaUpdateRule_2,
    DiffusionCaUpdateRule_3,
    DiffusionForce_2,
    DiffusionForce_3,
    DivisionBiasTrackingModifier_2,
    DivisionBiasTrackingModifier_3,
    ExclusionCaBasedDivisionRule_2,
    ExclusionCaBasedDivisionRule_3,
    ExponentialG1GenerationalCellCycleModel,
    ExtrinsicPullModifier_2,
    ExtrinsicPullModifier_3,
    FarhadifarForce_2,
    FarhadifarForce_3,
    FixedCentreBasedDivisionRule_2_2,
    FixedCentreBasedDivisionRule_3_3,
    FixedG1GenerationalCellCycleModel,
    FixedSequenceCellCycleModel,
    FixedVertexBasedDivisionRule_2,
    FixedVertexBasedDivisionRule_3,
    ForwardEulerNumericalMethod_2_2,
    ForwardEulerNumericalMethod_3_3,
    GammaG1CellCycleModel,
    GeneralisedLinearSpringForce_2_2,
    GeneralisedLinearSpringForce_3_3,
    Goldbeter1991SrnModel,
    HeterotypicBoundaryLengthWriter_2_2,
    HeterotypicBoundaryLengthWriter_3_3,
    ImmersedBoundaryBoundaryCellWriter_2_2,
    ImmersedBoundaryBoundaryCellWriter_3_3,
    ImmersedBoundaryCellPopulation_2,
    ImmersedBoundaryCellPopulation_3,
    ImmersedBoundaryKinematicFeedbackForce_2,
    ImmersedBoundaryKinematicFeedbackForce_3,
    ImmersedBoundaryLinearDifferentialAdhesionForce_2,
    ImmersedBoundaryLinearDifferentialAdhesionForce_3,
    ImmersedBoundaryLinearInteractionForce_2,
    ImmersedBoundaryLinearInteractionForce_3,
    ImmersedBoundaryLinearMembraneForce_2,
    ImmersedBoundaryLinearMembraneForce_3,
    ImmersedBoundaryMorseInteractionForce_2,
    ImmersedBoundaryMorseInteractionForce_3,
    ImmersedBoundaryMorseMembraneForce_2,
    ImmersedBoundaryMorseMembraneForce_3,
    ImmersedBoundaryNeighbourNumberWriter_2_2,
    ImmersedBoundaryNeighbourNumberWriter_3_3,
    ImmersedBoundarySimulationModifier_2,
    ImmersedBoundarySimulationModifier_3,
    ImmersedBoundarySvgWriter_2,
    ImmersedBoundarySvgWriter_3,
    IsolatedLabelledCellKiller_2,
    IsolatedLabelledCellKiller_3,
    LabelDependentBernoulliTrialCellCycleModel,
    LegacyCellProliferativeTypesWriter_2_2,
    LegacyCellProliferativeTypesWriter_3_3,
    MeshBasedCellPopulation_2_2,
    MeshBasedCellPopulation_3_3,
    MeshBasedCellPopulationWithGhostNodes_2,
    MeshBasedCellPopulationWithGhostNodes_3,
    NagaiHondaDifferentialAdhesionForce_2,
    NagaiHondaDifferentialAdhesionForce_3,
    NagaiHondaForce_2,
    NagaiHondaForce_3,
    NoCellCycleModel,
    NodeBasedCellPopulation_2,
    NodeBasedCellPopulation_3,
    NodeBasedCellPopulationWithBuskeUpdate_2,
    NodeBasedCellPopulationWithBuskeUpdate_3,
    NodeBasedCellPopulationWithParticles_2,
    NodeBasedCellPopulationWithParticles_3,
    NodeLocationWriter_2_2,
    NodeLocationWriter_3_3,
    NodeVelocityWriter_2_2,
    NodeVelocityWriter_3_3,
    NormallyDistributedTargetAreaModifier_2,
    NormallyDistributedTargetAreaModifier_3,
    NullSrnModel,
    OffLatticeSimulation_2_2,
    OffLatticeSimulation_3_3,
    OnLatticeSimulation_2,
    OnLatticeSimulation_3,
    PlanarPolarisedFarhadifarForce_2,
    PlanarPolarisedFarhadifarForce_3,
    PlaneBasedCellKiller_2,
    PlaneBasedCellKiller_3,
    PlaneBoundaryCondition_2_2,
    PlaneBoundaryCondition_3_3,
    PottsBasedCellPopulation_2,
    PottsBasedCellPopulation_3,
    PythonSimulationModifier_2,
    PythonSimulationModifier_3,
    RadialCellDataDistributionWriter_2_2,
    RadialCellDataDistributionWriter_3_3,
    RandomCaSwitchingUpdateRule_2,
    RandomCaSwitchingUpdateRule_3,
    RandomCellKiller_2,
    RandomCellKiller_3,
    RandomDirectionCentreBasedDivisionRule_2_2,
    RandomDirectionCentreBasedDivisionRule_3_3,
    RandomDirectionVertexBasedDivisionRule_2,
    RandomDirectionVertexBasedDivisionRule_3,
    RepulsionForce_2,
    RepulsionForce_3,
    ShortAxisImmersedBoundaryDivisionRule_2,
    ShortAxisImmersedBoundaryDivisionRule_3,
    ShortAxisVertexBasedDivisionRule_2,
    ShortAxisVertexBasedDivisionRule_3,
    ShovingCaBasedDivisionRule_2,
    ShovingCaBasedDivisionRule_3,
    SimpleOxygenBasedCellCycleModel,
    SimpleTargetAreaModifier_2,
    SimpleTargetAreaModifier_3,
    SimulationTime,
    SlidingBoundaryCondition_2,
    SlidingBoundaryCondition_3,
    SphereGeometryBoundaryCondition_2,
    SphereGeometryBoundaryCondition_3,
    StemCellProliferativeType,
    StochasticOxygenBasedCellCycleModel,
    SurfaceAreaConstraintPottsUpdateRule_2,
    SurfaceAreaConstraintPottsUpdateRule_3,
    T2SwapCellKiller_2,
    T2SwapCellKiller_3,
    TargetAreaLinearGrowthModifier_2,
    TargetAreaLinearGrowthModifier_3,
    TargetedCellKiller_2,
    TargetedCellKiller_3,
    TransitCellProliferativeType,
    TysonNovakCellCycleModel,
    UniformCellCycleModel,
    UniformG1GenerationalCellCycleModel,
    VertexBasedCellPopulation_2,
    VertexBasedCellPopulation_3,
    VertexBasedPopulationSrn_2,
    VertexBasedPopulationSrn_3,
    VertexIntersectionSwapLocationsWriter_2_2,
    VertexIntersectionSwapLocationsWriter_3_3,
    VertexT1SwapLocationsWriter_2_2,
    VertexT1SwapLocationsWriter_3_3,
    VertexT2SwapLocationsWriter_2_2,
    VertexT2SwapLocationsWriter_3_3,
    VertexT3SwapLocationsWriter_2_2,
    VertexT3SwapLocationsWriter_3_3,
    VolumeConstraintPottsUpdateRule_2,
    VolumeConstraintPottsUpdateRule_3,
    VolumeTrackingModifier_2,
    VolumeTrackingModifier_3,
    VonMisesVertexBasedDivisionRule_2,
    VonMisesVertexBasedDivisionRule_3,
    VoronoiDataWriter_2_2,
    VoronoiDataWriter_3_3,
    VtkSceneModifier_2,
    VtkSceneModifier_3,
    WelikyOsterForce_2,
    WelikyOsterForce_3,
    WildTypeCellMutationState,
)
from chaste._syntax import TemplatedClass
from chaste.cell_based._fortests import (
    AbstractCellBasedTestSuite,
    AbstractCellBasedWithTimingsTestSuite,
    SetupNotebookTest,
    TearDownNotebookTest,
)

AdhesionPottsUpdateRule = TemplatedClass(
    {
        ("2",): AdhesionPottsUpdateRule_2,
        ("3",): AdhesionPottsUpdateRule_3,
    }
)

ApoptoticCellKiller = TemplatedClass(
    {
        ("2",): ApoptoticCellKiller_2,
        ("3",): ApoptoticCellKiller_3,
    }
)

AttractingPlaneBoundaryCondition = TemplatedClass(
    {
        ("2",): AttractingPlaneBoundaryCondition_2_2,
        ("2", "2"): AttractingPlaneBoundaryCondition_2_2,
        ("3",): AttractingPlaneBoundaryCondition_3_3,
        ("3", "3"): AttractingPlaneBoundaryCondition_3_3,
    }
)

BoundaryNodeWriter = TemplatedClass(
    {
        ("2",): BoundaryNodeWriter_2_2,
        ("2", "2"): BoundaryNodeWriter_2_2,
        ("3",): BoundaryNodeWriter_3_3,
        ("3", "3"): BoundaryNodeWriter_3_3,
    }
)

BuskeAdhesiveForce = TemplatedClass(
    {
        ("2",): BuskeAdhesiveForce_2,
        ("3",): BuskeAdhesiveForce_3,
    }
)

BuskeCompressionForce = TemplatedClass(
    {
        ("2",): BuskeCompressionForce_2,
        ("3",): BuskeCompressionForce_3,
    }
)

BuskeElasticForce = TemplatedClass(
    {
        ("2",): BuskeElasticForce_2,
        ("3",): BuskeElasticForce_3,
    }
)

CaBasedCellPopulation = TemplatedClass(
    {
        ("2",): CaBasedCellPopulation_2,
        ("3",): CaBasedCellPopulation_3,
    }
)

CellAgesWriter = TemplatedClass(
    {
        ("2",): CellAgesWriter_2_2,
        ("2", "2"): CellAgesWriter_2_2,
        ("3",): CellAgesWriter_3_3,
        ("3", "3"): CellAgesWriter_3_3,
    }
)

CellAncestorWriter = TemplatedClass(
    {
        ("2",): CellAncestorWriter_2_2,
        ("2", "2"): CellAncestorWriter_2_2,
        ("3",): CellAncestorWriter_3_3,
        ("3", "3"): CellAncestorWriter_3_3,
    }
)

CellAppliedForceWriter = TemplatedClass(
    {
        ("2",): CellAppliedForceWriter_2_2,
        ("2", "2"): CellAppliedForceWriter_2_2,
        ("3",): CellAppliedForceWriter_3_3,
        ("3", "3"): CellAppliedForceWriter_3_3,
    }
)

CellCycleModelProteinConcentrationsWriter = TemplatedClass(
    {
        ("2",): CellCycleModelProteinConcentrationsWriter_2_2,
        ("2", "2"): CellCycleModelProteinConcentrationsWriter_2_2,
        ("3",): CellCycleModelProteinConcentrationsWriter_3_3,
        ("3", "3"): CellCycleModelProteinConcentrationsWriter_3_3,
    }
)

CellDataItemWriter = TemplatedClass(
    {
        ("2",): CellDataItemWriter_2_2,
        ("2", "2"): CellDataItemWriter_2_2,
        ("3",): CellDataItemWriter_3_3,
        ("3", "3"): CellDataItemWriter_3_3,
    }
)

CellDeltaNotchWriter = TemplatedClass(
    {
        ("2",): CellDeltaNotchWriter_2_2,
        ("2", "2"): CellDeltaNotchWriter_2_2,
        ("3",): CellDeltaNotchWriter_3_3,
        ("3", "3"): CellDeltaNotchWriter_3_3,
    }
)

CellDivisionLocationsWriter = TemplatedClass(
    {
        ("2",): CellDivisionLocationsWriter_2_2,
        ("2", "2"): CellDivisionLocationsWriter_2_2,
        ("3",): CellDivisionLocationsWriter_3_3,
        ("3", "3"): CellDivisionLocationsWriter_3_3,
    }
)

CellIdWriter = TemplatedClass(
    {
        ("2",): CellIdWriter_2_2,
        ("2", "2"): CellIdWriter_2_2,
        ("3",): CellIdWriter_3_3,
        ("3", "3"): CellIdWriter_3_3,
    }
)

CellLabelWriter = TemplatedClass(
    {
        ("2",): CellLabelWriter_2_2,
        ("2", "2"): CellLabelWriter_2_2,
        ("3",): CellLabelWriter_3_3,
        ("3", "3"): CellLabelWriter_3_3,
    }
)

CellLocationIndexWriter = TemplatedClass(
    {
        ("2",): CellLocationIndexWriter_2_2,
        ("2", "2"): CellLocationIndexWriter_2_2,
        ("3",): CellLocationIndexWriter_3_3,
        ("3", "3"): CellLocationIndexWriter_3_3,
    }
)

CellMutationStatesCountWriter = TemplatedClass(
    {
        ("2",): CellMutationStatesCountWriter_2_2,
        ("2", "2"): CellMutationStatesCountWriter_2_2,
        ("3",): CellMutationStatesCountWriter_3_3,
        ("3", "3"): CellMutationStatesCountWriter_3_3,
    }
)

CellMutationStatesWriter = TemplatedClass(
    {
        ("2",): CellMutationStatesWriter_2_2,
        ("2", "2"): CellMutationStatesWriter_2_2,
        ("3",): CellMutationStatesWriter_3_3,
        ("3", "3"): CellMutationStatesWriter_3_3,
    }
)

CellPopulationAdjacencyMatrixWriter = TemplatedClass(
    {
        ("2",): CellPopulationAdjacencyMatrixWriter_2_2,
        ("2", "2"): CellPopulationAdjacencyMatrixWriter_2_2,
        ("3",): CellPopulationAdjacencyMatrixWriter_3_3,
        ("3", "3"): CellPopulationAdjacencyMatrixWriter_3_3,
    }
)

CellPopulationAreaWriter = TemplatedClass(
    {
        ("2",): CellPopulationAreaWriter_2_2,
        ("2", "2"): CellPopulationAreaWriter_2_2,
        ("3",): CellPopulationAreaWriter_3_3,
        ("3", "3"): CellPopulationAreaWriter_3_3,
    }
)

CellPopulationElementWriter = TemplatedClass(
    {
        ("2",): CellPopulationElementWriter_2_2,
        ("2", "2"): CellPopulationElementWriter_2_2,
        ("3",): CellPopulationElementWriter_3_3,
        ("3", "3"): CellPopulationElementWriter_3_3,
    }
)

CellProliferativePhasesCountWriter = TemplatedClass(
    {
        ("2",): CellProliferativePhasesCountWriter_2_2,
        ("2", "2"): CellProliferativePhasesCountWriter_2_2,
        ("3",): CellProliferativePhasesCountWriter_3_3,
        ("3", "3"): CellProliferativePhasesCountWriter_3_3,
    }
)

CellProliferativePhasesWriter = TemplatedClass(
    {
        ("2",): CellProliferativePhasesWriter_2_2,
        ("2", "2"): CellProliferativePhasesWriter_2_2,
        ("3",): CellProliferativePhasesWriter_3_3,
        ("3", "3"): CellProliferativePhasesWriter_3_3,
    }
)

CellProliferativeTypesCountWriter = TemplatedClass(
    {
        ("2",): CellProliferativeTypesCountWriter_2_2,
        ("2", "2"): CellProliferativeTypesCountWriter_2_2,
        ("3",): CellProliferativeTypesCountWriter_3_3,
        ("3", "3"): CellProliferativeTypesCountWriter_3_3,
    }
)

CellProliferativeTypesWriter = TemplatedClass(
    {
        ("2",): CellProliferativeTypesWriter_2_2,
        ("2", "2"): CellProliferativeTypesWriter_2_2,
        ("3",): CellProliferativeTypesWriter_3_3,
        ("3", "3"): CellProliferativeTypesWriter_3_3,
    }
)

CellRadiusWriter = TemplatedClass(
    {
        ("2",): CellRadiusWriter_2_2,
        ("2", "2"): CellRadiusWriter_2_2,
        ("3",): CellRadiusWriter_3_3,
        ("3", "3"): CellRadiusWriter_3_3,
    }
)

CellRemovalLocationsWriter = TemplatedClass(
    {
        ("2",): CellRemovalLocationsWriter_2_2,
        ("2", "2"): CellRemovalLocationsWriter_2_2,
        ("3",): CellRemovalLocationsWriter_3_3,
        ("3", "3"): CellRemovalLocationsWriter_3_3,
    }
)

CellRosetteRankWriter = TemplatedClass(
    {
        ("2",): CellRosetteRankWriter_2_2,
        ("2", "2"): CellRosetteRankWriter_2_2,
        ("3",): CellRosetteRankWriter_3_3,
        ("3", "3"): CellRosetteRankWriter_3_3,
    }
)

CellsGenerator = TemplatedClass(
    {
        (
            "Alarcon2004OxygenBasedCellCycleModel",
            "2",
        ): CellsGenerator_Alarcon2004OxygenBasedCellCycleModel_2,
        (
            "Alarcon2004OxygenBasedCellCycleModel",
            "3",
        ): CellsGenerator_Alarcon2004OxygenBasedCellCycleModel_3,
        (
            "AlwaysDivideCellCycleModel",
            "2",
        ): CellsGenerator_AlwaysDivideCellCycleModel_2,
        (
            "AlwaysDivideCellCycleModel",
            "3",
        ): CellsGenerator_AlwaysDivideCellCycleModel_3,
        (
            "BernoulliTrialCellCycleModel",
            "2",
        ): CellsGenerator_BernoulliTrialCellCycleModel_2,
        (
            "BernoulliTrialCellCycleModel",
            "3",
        ): CellsGenerator_BernoulliTrialCellCycleModel_3,
        (
            "BiasedBernoulliTrialCellCycleModel",
            "2",
        ): CellsGenerator_BiasedBernoulliTrialCellCycleModel_2,
        (
            "BiasedBernoulliTrialCellCycleModel",
            "3",
        ): CellsGenerator_BiasedBernoulliTrialCellCycleModel_3,
        (
            "ContactInhibitionCellCycleModel",
            "2",
        ): CellsGenerator_ContactInhibitionCellCycleModel_2,
        (
            "ContactInhibitionCellCycleModel",
            "3",
        ): CellsGenerator_ContactInhibitionCellCycleModel_3,
        (
            "ExponentialG1GenerationalCellCycleModel",
            "2",
        ): CellsGenerator_ExponentialG1GenerationalCellCycleModel_2,
        (
            "ExponentialG1GenerationalCellCycleModel",
            "3",
        ): CellsGenerator_ExponentialG1GenerationalCellCycleModel_3,
        (
            "FixedG1GenerationalCellCycleModel",
            "2",
        ): CellsGenerator_FixedG1GenerationalCellCycleModel_2,
        (
            "FixedG1GenerationalCellCycleModel",
            "3",
        ): CellsGenerator_FixedG1GenerationalCellCycleModel_3,
        (
            "FixedSequenceCellCycleModel",
            "2",
        ): CellsGenerator_FixedSequenceCellCycleModel_2,
        (
            "FixedSequenceCellCycleModel",
            "3",
        ): CellsGenerator_FixedSequenceCellCycleModel_3,
        ("GammaG1CellCycleModel", "2"): CellsGenerator_GammaG1CellCycleModel_2,
        ("GammaG1CellCycleModel", "3"): CellsGenerator_GammaG1CellCycleModel_3,
        (
            "LabelDependentBernoulliTrialCellCycleModel",
            "2",
        ): CellsGenerator_LabelDependentBernoulliTrialCellCycleModel_2,
        (
            "LabelDependentBernoulliTrialCellCycleModel",
            "3",
        ): CellsGenerator_LabelDependentBernoulliTrialCellCycleModel_3,
        ("NoCellCycleModel", "2"): CellsGenerator_NoCellCycleModel_2,
        ("NoCellCycleModel", "3"): CellsGenerator_NoCellCycleModel_3,
        (
            "SimpleOxygenBasedCellCycleModel",
            "2",
        ): CellsGenerator_SimpleOxygenBasedCellCycleModel_2,
        (
            "SimpleOxygenBasedCellCycleModel",
            "3",
        ): CellsGenerator_SimpleOxygenBasedCellCycleModel_3,
        (
            "StochasticOxygenBasedCellCycleModel",
            "2",
        ): CellsGenerator_StochasticOxygenBasedCellCycleModel_2,
        (
            "StochasticOxygenBasedCellCycleModel",
            "3",
        ): CellsGenerator_StochasticOxygenBasedCellCycleModel_3,
        (
            "TysonNovakCellCycleModel",
            "2",
        ): CellsGenerator_TysonNovakCellCycleModel_2,
        (
            "TysonNovakCellCycleModel",
            "3",
        ): CellsGenerator_TysonNovakCellCycleModel_3,
        ("UniformCellCycleModel", "2"): CellsGenerator_UniformCellCycleModel_2,
        ("UniformCellCycleModel", "3"): CellsGenerator_UniformCellCycleModel_3,
        (
            "UniformG1GenerationalCellCycleModel",
            "2",
        ): CellsGenerator_UniformG1GenerationalCellCycleModel_2,
        (
            "UniformG1GenerationalCellCycleModel",
            "3",
        ): CellsGenerator_UniformG1GenerationalCellCycleModel_3,
    }
)

CellVolumesWriter = TemplatedClass(
    {
        ("2",): CellVolumesWriter_2_2,
        ("2", "2"): CellVolumesWriter_2_2,
        ("3",): CellVolumesWriter_3_3,
        ("3", "3"): CellVolumesWriter_3_3,
    }
)

ChemotacticForce = TemplatedClass(
    {
        ("2",): ChemotacticForce_2,
        ("3",): ChemotacticForce_3,
    }
)

ChemotaxisPottsUpdateRule = TemplatedClass(
    {
        ("2",): ChemotaxisPottsUpdateRule_2,
        ("3",): ChemotaxisPottsUpdateRule_3,
    }
)

DeltaNotchEdgeInteriorTrackingModifier = TemplatedClass(
    {
        ("2",): DeltaNotchEdgeInteriorTrackingModifier_2,
        ("3",): DeltaNotchEdgeInteriorTrackingModifier_3,
    }
)

DeltaNotchTrackingModifier = TemplatedClass(
    {
        ("2",): DeltaNotchTrackingModifier_2,
        ("3",): DeltaNotchTrackingModifier_3,
    }
)

DeltaNotchEdgeTrackingModifier = TemplatedClass(
    {
        ("2",): DeltaNotchEdgeTrackingModifier_2,
        ("3",): DeltaNotchEdgeTrackingModifier_3,
    }
)

DifferentialAdhesionGeneralisedLinearSpringForce = TemplatedClass(
    {
        ("2",): DifferentialAdhesionGeneralisedLinearSpringForce_2_2,
        ("2", "2"): DifferentialAdhesionGeneralisedLinearSpringForce_2_2,
        ("3",): DifferentialAdhesionGeneralisedLinearSpringForce_3_3,
        ("3", "3"): DifferentialAdhesionGeneralisedLinearSpringForce_3_3,
    }
)

DifferentialAdhesionPottsUpdateRule = TemplatedClass(
    {
        ("2",): DifferentialAdhesionPottsUpdateRule_2,
        ("3",): DifferentialAdhesionPottsUpdateRule_3,
    }
)

DiffusionCaUpdateRule = TemplatedClass(
    {
        ("2",): DiffusionCaUpdateRule_2,
        ("3",): DiffusionCaUpdateRule_3,
    }
)

DiffusionForce = TemplatedClass(
    {
        ("2",): DiffusionForce_2,
        ("3",): DiffusionForce_3,
    }
)

DivisionBiasTrackingModifier = TemplatedClass(
    {
        ("2",): DivisionBiasTrackingModifier_2,
        ("3",): DivisionBiasTrackingModifier_3,
    }
)

ExclusionCaBasedDivisionRule = TemplatedClass(
    {
        ("2",): ExclusionCaBasedDivisionRule_2,
        ("3",): ExclusionCaBasedDivisionRule_3,
    }
)

ExtrinsicPullModifier = TemplatedClass(
    {
        ("2",): ExtrinsicPullModifier_2,
        ("3",): ExtrinsicPullModifier_3,
    }
)

FarhadifarForce = TemplatedClass(
    {
        ("2",): FarhadifarForce_2,
        ("3",): FarhadifarForce_3,
    }
)

FixedCentreBasedDivisionRule = TemplatedClass(
    {
        ("2",): FixedCentreBasedDivisionRule_2_2,
        ("2", "2"): FixedCentreBasedDivisionRule_2_2,
        ("3",): FixedCentreBasedDivisionRule_3_3,
        ("3", "3"): FixedCentreBasedDivisionRule_3_3,
    }
)

FixedVertexBasedDivisionRule = TemplatedClass(
    {
        ("2",): FixedVertexBasedDivisionRule_2,
        ("3",): FixedVertexBasedDivisionRule_3,
    }
)

ForwardEulerNumericalMethod = TemplatedClass(
    {
        ("2",): ForwardEulerNumericalMethod_2_2,
        ("2", "2"): ForwardEulerNumericalMethod_2_2,
        ("3",): ForwardEulerNumericalMethod_3_3,
        ("3", "3"): ForwardEulerNumericalMethod_3_3,
    }
)

GeneralisedLinearSpringForce = TemplatedClass(
    {
        ("2",): GeneralisedLinearSpringForce_2_2,
        ("2", "2"): GeneralisedLinearSpringForce_2_2,
        ("3",): GeneralisedLinearSpringForce_3_3,
        ("3", "3"): GeneralisedLinearSpringForce_3_3,
    }
)

HeterotypicBoundaryLengthWriter = TemplatedClass(
    {
        ("2",): HeterotypicBoundaryLengthWriter_2_2,
        ("2", "2"): HeterotypicBoundaryLengthWriter_2_2,
        ("3",): HeterotypicBoundaryLengthWriter_3_3,
        ("3", "3"): HeterotypicBoundaryLengthWriter_3_3,
    }
)

ImmersedBoundaryBoundaryCellWriter = TemplatedClass(
    {
        ("2",): ImmersedBoundaryBoundaryCellWriter_2_2,
        ("2", "2"): ImmersedBoundaryBoundaryCellWriter_2_2,
        ("3",): ImmersedBoundaryBoundaryCellWriter_3_3,
        ("3", "3"): ImmersedBoundaryBoundaryCellWriter_3_3,
    }
)

ImmersedBoundaryCellPopulation = TemplatedClass(
    {
        ("2",): ImmersedBoundaryCellPopulation_2,
        ("3",): ImmersedBoundaryCellPopulation_3,
    }
)

ImmersedBoundaryKinematicFeedbackForce = TemplatedClass(
    {
        ("2",): ImmersedBoundaryKinematicFeedbackForce_2,
        ("3",): ImmersedBoundaryKinematicFeedbackForce_3,
    }
)

ImmersedBoundaryLinearDifferentialAdhesionForce = TemplatedClass(
    {
        ("2",): ImmersedBoundaryLinearDifferentialAdhesionForce_2,
        ("3",): ImmersedBoundaryLinearDifferentialAdhesionForce_3,
    }
)

ImmersedBoundaryLinearInteractionForce = TemplatedClass(
    {
        ("2",): ImmersedBoundaryLinearInteractionForce_2,
        ("3",): ImmersedBoundaryLinearInteractionForce_3,
    }
)

ImmersedBoundaryLinearMembraneForce = TemplatedClass(
    {
        ("2",): ImmersedBoundaryLinearMembraneForce_2,
        ("3",): ImmersedBoundaryLinearMembraneForce_3,
    }
)

ImmersedBoundaryMorseInteractionForce = TemplatedClass(
    {
        ("2",): ImmersedBoundaryMorseInteractionForce_2,
        ("3",): ImmersedBoundaryMorseInteractionForce_3,
    }
)

ImmersedBoundaryMorseMembraneForce = TemplatedClass(
    {
        ("2",): ImmersedBoundaryMorseMembraneForce_2,
        ("3",): ImmersedBoundaryMorseMembraneForce_3,
    }
)

ImmersedBoundaryNeighbourNumberWriter = TemplatedClass(
    {
        ("2",): ImmersedBoundaryNeighbourNumberWriter_2_2,
        ("2", "2"): ImmersedBoundaryNeighbourNumberWriter_2_2,
        ("3",): ImmersedBoundaryNeighbourNumberWriter_3_3,
        ("3", "3"): ImmersedBoundaryNeighbourNumberWriter_3_3,
    }
)

ImmersedBoundarySimulationModifier = TemplatedClass(
    {
        ("2",): ImmersedBoundarySimulationModifier_2,
        ("3",): ImmersedBoundarySimulationModifier_3,
    }
)

ImmersedBoundarySvgWriter = TemplatedClass(
    {
        ("2",): ImmersedBoundarySvgWriter_2,
        ("3",): ImmersedBoundarySvgWriter_3,
    }
)

IsolatedLabelledCellKiller = TemplatedClass(
    {
        ("2",): IsolatedLabelledCellKiller_2,
        ("3",): IsolatedLabelledCellKiller_3,
    }
)

LegacyCellProliferativeTypesWriter = TemplatedClass(
    {
        ("2",): LegacyCellProliferativeTypesWriter_2_2,
        ("2", "2"): LegacyCellProliferativeTypesWriter_2_2,
        ("3",): LegacyCellProliferativeTypesWriter_3_3,
        ("3", "3"): LegacyCellProliferativeTypesWriter_3_3,
    }
)

MeshBasedCellPopulation = TemplatedClass(
    {
        ("2",): MeshBasedCellPopulation_2_2,
        ("2", "2"): MeshBasedCellPopulation_2_2,
        ("3",): MeshBasedCellPopulation_3_3,
        ("3", "3"): MeshBasedCellPopulation_3_3,
    }
)

MeshBasedCellPopulationWithGhostNodes = TemplatedClass(
    {
        ("2",): MeshBasedCellPopulationWithGhostNodes_2,
        ("3",): MeshBasedCellPopulationWithGhostNodes_3,
    }
)

NagaiHondaDifferentialAdhesionForce = TemplatedClass(
    {
        ("2",): NagaiHondaDifferentialAdhesionForce_2,
        ("3",): NagaiHondaDifferentialAdhesionForce_3,
    }
)

NagaiHondaForce = TemplatedClass(
    {
        ("2",): NagaiHondaForce_2,
        ("3",): NagaiHondaForce_3,
    }
)

NodeBasedCellPopulation = TemplatedClass(
    {
        ("2",): NodeBasedCellPopulation_2,
        ("3",): NodeBasedCellPopulation_3,
    }
)

NodeBasedCellPopulationWithBuskeUpdate = TemplatedClass(
    {
        ("2",): NodeBasedCellPopulationWithBuskeUpdate_2,
        ("3",): NodeBasedCellPopulationWithBuskeUpdate_3,
    }
)

NodeBasedCellPopulationWithParticles = TemplatedClass(
    {
        ("2",): NodeBasedCellPopulationWithParticles_2,
        ("3",): NodeBasedCellPopulationWithParticles_3,
    }
)

NodeLocationWriter = TemplatedClass(
    {
        ("2",): NodeLocationWriter_2_2,
        ("2", "2"): NodeLocationWriter_2_2,
        ("3",): NodeLocationWriter_3_3,
        ("3", "3"): NodeLocationWriter_3_3,
    }
)

NodeVelocityWriter = TemplatedClass(
    {
        ("2",): NodeVelocityWriter_2_2,
        ("2", "2"): NodeVelocityWriter_2_2,
        ("3",): NodeVelocityWriter_3_3,
        ("3", "3"): NodeVelocityWriter_3_3,
    }
)

NormallyDistributedTargetAreaModifier = TemplatedClass(
    {
        ("2",): NormallyDistributedTargetAreaModifier_2,
        ("3",): NormallyDistributedTargetAreaModifier_3,
    }
)

OffLatticeSimulation = TemplatedClass(
    {
        ("2",): OffLatticeSimulation_2_2,
        ("2", "2"): OffLatticeSimulation_2_2,
        ("3",): OffLatticeSimulation_3_3,
        ("3", "3"): OffLatticeSimulation_3_3,
    }
)

OnLatticeSimulation = TemplatedClass(
    {
        ("2",): OnLatticeSimulation_2,
        ("3",): OnLatticeSimulation_3,
    }
)

PlanarPolarisedFarhadifarForce = TemplatedClass(
    {
        ("2",): PlanarPolarisedFarhadifarForce_2,
        ("3",): PlanarPolarisedFarhadifarForce_3,
    }
)

PlaneBasedCellKiller = TemplatedClass(
    {
        ("2",): PlaneBasedCellKiller_2,
        ("3",): PlaneBasedCellKiller_3,
    }
)

PlaneBoundaryCondition = TemplatedClass(
    {
        ("2",): PlaneBoundaryCondition_2_2,
        ("2", "2"): PlaneBoundaryCondition_2_2,
        ("3",): PlaneBoundaryCondition_3_3,
        ("3", "3"): PlaneBoundaryCondition_3_3,
    }
)

PottsBasedCellPopulation = TemplatedClass(
    {
        ("2",): PottsBasedCellPopulation_2,
        ("3",): PottsBasedCellPopulation_3,
    }
)

PythonSimulationModifier = TemplatedClass(
    {
        ("2",): PythonSimulationModifier_2,
        ("3",): PythonSimulationModifier_3,
    }
)

RadialCellDataDistributionWriter = TemplatedClass(
    {
        ("2",): RadialCellDataDistributionWriter_2_2,
        ("2", "2"): RadialCellDataDistributionWriter_2_2,
        ("3",): RadialCellDataDistributionWriter_3_3,
        ("3", "3"): RadialCellDataDistributionWriter_3_3,
    }
)

RandomCaSwitchingUpdateRule = TemplatedClass(
    {
        ("2",): RandomCaSwitchingUpdateRule_2,
        ("3",): RandomCaSwitchingUpdateRule_3,
    }
)

RandomCellKiller = TemplatedClass(
    {
        ("2",): RandomCellKiller_2,
        ("3",): RandomCellKiller_3,
    }
)

RandomDirectionCentreBasedDivisionRule = TemplatedClass(
    {
        ("2",): RandomDirectionCentreBasedDivisionRule_2_2,
        ("2", "2"): RandomDirectionCentreBasedDivisionRule_2_2,
        ("3",): RandomDirectionCentreBasedDivisionRule_3_3,
        ("3", "3"): RandomDirectionCentreBasedDivisionRule_3_3,
    }
)

RandomDirectionVertexBasedDivisionRule = TemplatedClass(
    {
        ("2",): RandomDirectionVertexBasedDivisionRule_2,
        ("3",): RandomDirectionVertexBasedDivisionRule_3,
    }
)

RepulsionForce = TemplatedClass(
    {
        ("2",): RepulsionForce_2,
        ("3",): RepulsionForce_3,
    }
)

ShortAxisImmersedBoundaryDivisionRule = TemplatedClass(
    {
        ("2",): ShortAxisImmersedBoundaryDivisionRule_2,
        ("3",): ShortAxisImmersedBoundaryDivisionRule_3,
    }
)

ShortAxisVertexBasedDivisionRule = TemplatedClass(
    {
        ("2",): ShortAxisVertexBasedDivisionRule_2,
        ("3",): ShortAxisVertexBasedDivisionRule_3,
    }
)

ShovingCaBasedDivisionRule = TemplatedClass(
    {
        ("2",): ShovingCaBasedDivisionRule_2,
        ("3",): ShovingCaBasedDivisionRule_3,
    }
)

SimpleTargetAreaModifier = TemplatedClass(
    {
        ("2",): SimpleTargetAreaModifier_2,
        ("3",): SimpleTargetAreaModifier_3,
    }
)

SlidingBoundaryCondition = TemplatedClass(
    {
        ("2",): SlidingBoundaryCondition_2,
        ("3",): SlidingBoundaryCondition_3,
    }
)

SphereGeometryBoundaryCondition = TemplatedClass(
    {
        ("2",): SphereGeometryBoundaryCondition_2,
        ("3",): SphereGeometryBoundaryCondition_3,
    }
)

SurfaceAreaConstraintPottsUpdateRule = TemplatedClass(
    {
        ("2",): SurfaceAreaConstraintPottsUpdateRule_2,
        ("3",): SurfaceAreaConstraintPottsUpdateRule_3,
    }
)

T2SwapCellKiller = TemplatedClass(
    {
        ("2",): T2SwapCellKiller_2,
        ("3",): T2SwapCellKiller_3,
    }
)

TargetAreaLinearGrowthModifier = TemplatedClass(
    {
        ("2",): TargetAreaLinearGrowthModifier_2,
        ("3",): TargetAreaLinearGrowthModifier_3,
    }
)

TargetedCellKiller = TemplatedClass(
    {
        ("2",): TargetedCellKiller_2,
        ("3",): TargetedCellKiller_3,
    }
)

VertexBasedCellPopulation = TemplatedClass(
    {
        ("2",): VertexBasedCellPopulation_2,
        ("3",): VertexBasedCellPopulation_3,
    }
)

VertexBasedPopulationSrn = TemplatedClass(
    {
        ("2",): VertexBasedPopulationSrn_2,
        ("3",): VertexBasedPopulationSrn_3,
    }
)

VertexIntersectionSwapLocationsWriter = TemplatedClass(
    {
        ("2",): VertexIntersectionSwapLocationsWriter_2_2,
        ("2", "2"): VertexIntersectionSwapLocationsWriter_2_2,
        ("3",): VertexIntersectionSwapLocationsWriter_3_3,
        ("3", "3"): VertexIntersectionSwapLocationsWriter_3_3,
    }
)

VertexT1SwapLocationsWriter = TemplatedClass(
    {
        ("2",): VertexT1SwapLocationsWriter_2_2,
        ("2", "2"): VertexT1SwapLocationsWriter_2_2,
        ("3",): VertexT1SwapLocationsWriter_3_3,
        ("3", "3"): VertexT1SwapLocationsWriter_3_3,
    }
)

VertexT2SwapLocationsWriter = TemplatedClass(
    {
        ("2",): VertexT2SwapLocationsWriter_2_2,
        ("2", "2"): VertexT2SwapLocationsWriter_2_2,
        ("3",): VertexT2SwapLocationsWriter_3_3,
        ("3", "3"): VertexT2SwapLocationsWriter_3_3,
    }
)

VertexT3SwapLocationsWriter = TemplatedClass(
    {
        ("2",): VertexT3SwapLocationsWriter_2_2,
        ("2", "2"): VertexT3SwapLocationsWriter_2_2,
        ("3",): VertexT3SwapLocationsWriter_3_3,
        ("3", "3"): VertexT3SwapLocationsWriter_3_3,
    }
)

VolumeConstraintPottsUpdateRule = TemplatedClass(
    {
        ("2",): VolumeConstraintPottsUpdateRule_2,
        ("3",): VolumeConstraintPottsUpdateRule_3,
    }
)

VolumeTrackingModifier = TemplatedClass(
    {
        ("2",): VolumeTrackingModifier_2,
        ("3",): VolumeTrackingModifier_3,
    }
)

VonMisesVertexBasedDivisionRule = TemplatedClass(
    {
        ("2",): VonMisesVertexBasedDivisionRule_2,
        ("3",): VonMisesVertexBasedDivisionRule_3,
    }
)

VoronoiDataWriter = TemplatedClass(
    {
        ("2",): VoronoiDataWriter_2_2,
        ("2", "2"): VoronoiDataWriter_2_2,
        ("3",): VoronoiDataWriter_3_3,
        ("3", "3"): VoronoiDataWriter_3_3,
    }
)

VtkSceneModifier = TemplatedClass(
    {
        ("2",): VtkSceneModifier_2,
        ("3",): VtkSceneModifier_3,
    }
)

WelikyOsterForce = TemplatedClass(
    {
        ("2",): WelikyOsterForce_2,
        ("3",): WelikyOsterForce_3,
    }
)
