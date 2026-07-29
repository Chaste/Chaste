"""Cell-Based Module"""

__copyright__ = """Copyright (c) 2005-2026, University of Oxford.
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
    AdhesionPottsUpdateRule_1,
    AdhesionPottsUpdateRule_2,
    AdhesionPottsUpdateRule_3,
    Alarcon2004OxygenBasedCellCycleModel,
    AlwaysDivideCellCycleModel,
    ApcOneHitCellMutationState,
    ApcTwoHitCellMutationState,
    ApoptoticCellKiller_1,
    ApoptoticCellKiller_2,
    ApoptoticCellKiller_3,
    ApoptoticCellProperty,
    AttractingPlaneBoundaryCondition_1_1,
    AttractingPlaneBoundaryCondition_1_2,
    AttractingPlaneBoundaryCondition_1_3,
    AttractingPlaneBoundaryCondition_2_2,
    AttractingPlaneBoundaryCondition_2_3,
    AttractingPlaneBoundaryCondition_3_3,
    BernoulliTrialCellCycleModel,
    BetaCateninOneHitCellMutationState,
    BiasedBernoulliTrialCellCycleModel,
    BoundaryNodeWriter_1_1,
    BoundaryNodeWriter_1_2,
    BoundaryNodeWriter_1_3,
    BoundaryNodeWriter_2_2,
    BoundaryNodeWriter_2_3,
    BoundaryNodeWriter_3_3,
    BuskeAdhesiveForce_1,
    BuskeAdhesiveForce_2,
    BuskeAdhesiveForce_3,
    BuskeCompressionForce_1,
    BuskeCompressionForce_2,
    BuskeCompressionForce_3,
    BuskeElasticForce_1,
    BuskeElasticForce_2,
    BuskeElasticForce_3,
    CaBasedCellPopulation_1,
    CaBasedCellPopulation_2,
    CaBasedCellPopulation_3,
    Cell,
    CellAgesWriter_1_1,
    CellAgesWriter_1_2,
    CellAgesWriter_1_3,
    CellAgesWriter_2_2,
    CellAgesWriter_2_3,
    CellAgesWriter_3_3,
    CellAncestor,
    CellAncestorWriter_1_1,
    CellAncestorWriter_1_2,
    CellAncestorWriter_1_3,
    CellAncestorWriter_2_2,
    CellAncestorWriter_2_3,
    CellAncestorWriter_3_3,
    CellAppliedForceWriter_1_1,
    CellAppliedForceWriter_1_2,
    CellAppliedForceWriter_1_3,
    CellAppliedForceWriter_2_2,
    CellAppliedForceWriter_2_3,
    CellAppliedForceWriter_3_3,
    CellCycleModelProteinConcentrationsWriter_1_1,
    CellCycleModelProteinConcentrationsWriter_1_2,
    CellCycleModelProteinConcentrationsWriter_1_3,
    CellCycleModelProteinConcentrationsWriter_2_2,
    CellCycleModelProteinConcentrationsWriter_2_3,
    CellCycleModelProteinConcentrationsWriter_3_3,
    CellData,
    CellDataItemWriter_1_1,
    CellDataItemWriter_1_2,
    CellDataItemWriter_1_3,
    CellDataItemWriter_2_2,
    CellDataItemWriter_2_3,
    CellDataItemWriter_3_3,
    CellDeltaNotchWriter_1_1,
    CellDeltaNotchWriter_1_2,
    CellDeltaNotchWriter_1_3,
    CellDeltaNotchWriter_2_2,
    CellDeltaNotchWriter_2_3,
    CellDeltaNotchWriter_3_3,
    CellDivisionLocationsWriter_1_1,
    CellDivisionLocationsWriter_1_2,
    CellDivisionLocationsWriter_1_3,
    CellDivisionLocationsWriter_2_2,
    CellDivisionLocationsWriter_2_3,
    CellDivisionLocationsWriter_3_3,
    CellEdgeData,
    CellId,
    CellIdWriter_1_1,
    CellIdWriter_1_2,
    CellIdWriter_1_3,
    CellIdWriter_2_2,
    CellIdWriter_2_3,
    CellIdWriter_3_3,
    CellLabel,
    CellLabelWriter_1_1,
    CellLabelWriter_1_2,
    CellLabelWriter_1_3,
    CellLabelWriter_2_2,
    CellLabelWriter_2_3,
    CellLabelWriter_3_3,
    CellLocationIndexWriter_1_1,
    CellLocationIndexWriter_1_2,
    CellLocationIndexWriter_1_3,
    CellLocationIndexWriter_2_2,
    CellLocationIndexWriter_2_3,
    CellLocationIndexWriter_3_3,
    CellMutationStatesCountWriter_1_1,
    CellMutationStatesCountWriter_1_2,
    CellMutationStatesCountWriter_1_3,
    CellMutationStatesCountWriter_2_2,
    CellMutationStatesCountWriter_2_3,
    CellMutationStatesCountWriter_3_3,
    CellMutationStatesWriter_1_1,
    CellMutationStatesWriter_1_2,
    CellMutationStatesWriter_1_3,
    CellMutationStatesWriter_2_2,
    CellMutationStatesWriter_2_3,
    CellMutationStatesWriter_3_3,
    CellPopulationAdjacencyMatrixWriter_1_1,
    CellPopulationAdjacencyMatrixWriter_1_2,
    CellPopulationAdjacencyMatrixWriter_1_3,
    CellPopulationAdjacencyMatrixWriter_2_2,
    CellPopulationAdjacencyMatrixWriter_2_3,
    CellPopulationAdjacencyMatrixWriter_3_3,
    CellPopulationAreaWriter_1_1,
    CellPopulationAreaWriter_1_2,
    CellPopulationAreaWriter_1_3,
    CellPopulationAreaWriter_2_2,
    CellPopulationAreaWriter_2_3,
    CellPopulationAreaWriter_3_3,
    CellPopulationElementWriter_1_1,
    CellPopulationElementWriter_1_2,
    CellPopulationElementWriter_1_3,
    CellPopulationElementWriter_2_2,
    CellPopulationElementWriter_2_3,
    CellPopulationElementWriter_3_3,
    CellProliferativePhasesCountWriter_1_1,
    CellProliferativePhasesCountWriter_1_2,
    CellProliferativePhasesCountWriter_1_3,
    CellProliferativePhasesCountWriter_2_2,
    CellProliferativePhasesCountWriter_2_3,
    CellProliferativePhasesCountWriter_3_3,
    CellProliferativePhasesWriter_1_1,
    CellProliferativePhasesWriter_1_2,
    CellProliferativePhasesWriter_1_3,
    CellProliferativePhasesWriter_2_2,
    CellProliferativePhasesWriter_2_3,
    CellProliferativePhasesWriter_3_3,
    CellProliferativeTypesCountWriter_1_1,
    CellProliferativeTypesCountWriter_1_2,
    CellProliferativeTypesCountWriter_1_3,
    CellProliferativeTypesCountWriter_2_2,
    CellProliferativeTypesCountWriter_2_3,
    CellProliferativeTypesCountWriter_3_3,
    CellProliferativeTypesWriter_1_1,
    CellProliferativeTypesWriter_1_2,
    CellProliferativeTypesWriter_1_3,
    CellProliferativeTypesWriter_2_2,
    CellProliferativeTypesWriter_2_3,
    CellProliferativeTypesWriter_3_3,
    CellPropertyCollection,
    CellPropertyRegistry,
    CellRadiusWriter_1_1,
    CellRadiusWriter_1_2,
    CellRadiusWriter_1_3,
    CellRadiusWriter_2_2,
    CellRadiusWriter_2_3,
    CellRadiusWriter_3_3,
    CellRemovalLocationsWriter_1_1,
    CellRemovalLocationsWriter_1_2,
    CellRemovalLocationsWriter_1_3,
    CellRemovalLocationsWriter_2_2,
    CellRemovalLocationsWriter_2_3,
    CellRemovalLocationsWriter_3_3,
    CellRosetteRankWriter_1_1,
    CellRosetteRankWriter_1_2,
    CellRosetteRankWriter_1_3,
    CellRosetteRankWriter_2_2,
    CellRosetteRankWriter_2_3,
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
    CellVolumesWriter_1_1,
    CellVolumesWriter_1_2,
    CellVolumesWriter_1_3,
    CellVolumesWriter_2_2,
    CellVolumesWriter_2_3,
    CellVolumesWriter_3_3,
    ChemotacticForce_1,
    ChemotacticForce_2,
    ChemotacticForce_3,
    ChemotaxisPottsUpdateRule_1,
    ChemotaxisPottsUpdateRule_2,
    ChemotaxisPottsUpdateRule_3,
    ContactInhibitionCellCycleModel,
    DefaultCellProliferativeType,
    DeltaNotchEdgeInteriorTrackingModifier_1,
    DeltaNotchEdgeInteriorTrackingModifier_2,
    DeltaNotchEdgeInteriorTrackingModifier_3,
    DeltaNotchEdgeSrnModel,
    DeltaNotchEdgeTrackingModifier_1,
    DeltaNotchEdgeTrackingModifier_2,
    DeltaNotchEdgeTrackingModifier_3,
    DeltaNotchInteriorSrnModel,
    DeltaNotchSrnModel,
    DeltaNotchTrackingModifier_1,
    DeltaNotchTrackingModifier_2,
    DeltaNotchTrackingModifier_3,
    DifferentialAdhesionGeneralisedLinearSpringForce_1_1,
    DifferentialAdhesionGeneralisedLinearSpringForce_1_2,
    DifferentialAdhesionGeneralisedLinearSpringForce_1_3,
    DifferentialAdhesionGeneralisedLinearSpringForce_2_2,
    DifferentialAdhesionGeneralisedLinearSpringForce_2_3,
    DifferentialAdhesionGeneralisedLinearSpringForce_3_3,
    DifferentialAdhesionPottsUpdateRule_1,
    DifferentialAdhesionPottsUpdateRule_2,
    DifferentialAdhesionPottsUpdateRule_3,
    DifferentiatedCellProliferativeType,
    DiffusionCaUpdateRule_1,
    DiffusionCaUpdateRule_2,
    DiffusionCaUpdateRule_3,
    DiffusionForce_1,
    DiffusionForce_2,
    DiffusionForce_3,
    DivisionBiasTrackingModifier_1,
    DivisionBiasTrackingModifier_2,
    DivisionBiasTrackingModifier_3,
    ExclusionCaBasedDivisionRule_1,
    ExclusionCaBasedDivisionRule_2,
    ExclusionCaBasedDivisionRule_3,
    ExponentialG1GenerationalCellCycleModel,
    ExtrinsicPullModifier_1,
    ExtrinsicPullModifier_2,
    ExtrinsicPullModifier_3,
    FarhadifarForce_1,
    FarhadifarForce_2,
    FarhadifarForce_3,
    FixedCentreBasedDivisionRule_1_1,
    FixedCentreBasedDivisionRule_1_2,
    FixedCentreBasedDivisionRule_1_3,
    FixedCentreBasedDivisionRule_2_2,
    FixedCentreBasedDivisionRule_2_3,
    FixedCentreBasedDivisionRule_3_3,
    FixedG1GenerationalCellCycleModel,
    FixedSequenceCellCycleModel,
    FixedVertexBasedDivisionRule_1,
    FixedVertexBasedDivisionRule_2,
    FixedVertexBasedDivisionRule_3,
    ForwardEulerNumericalMethod_1_1,
    ForwardEulerNumericalMethod_1_2,
    ForwardEulerNumericalMethod_1_3,
    ForwardEulerNumericalMethod_2_2,
    ForwardEulerNumericalMethod_2_3,
    ForwardEulerNumericalMethod_3_3,
    GammaG1CellCycleModel,
    GeneralisedLinearSpringForce_1_1,
    GeneralisedLinearSpringForce_1_2,
    GeneralisedLinearSpringForce_1_3,
    GeneralisedLinearSpringForce_2_2,
    GeneralisedLinearSpringForce_2_3,
    GeneralisedLinearSpringForce_3_3,
    Goldbeter1991SrnModel,
    HeterotypicBoundaryLengthWriter_1_1,
    HeterotypicBoundaryLengthWriter_1_2,
    HeterotypicBoundaryLengthWriter_1_3,
    HeterotypicBoundaryLengthWriter_2_2,
    HeterotypicBoundaryLengthWriter_2_3,
    HeterotypicBoundaryLengthWriter_3_3,
    ImmersedBoundaryBoundaryCellWriter_1_1,
    ImmersedBoundaryBoundaryCellWriter_1_2,
    ImmersedBoundaryBoundaryCellWriter_1_3,
    ImmersedBoundaryBoundaryCellWriter_2_2,
    ImmersedBoundaryBoundaryCellWriter_2_3,
    ImmersedBoundaryBoundaryCellWriter_3_3,
    ImmersedBoundaryCellPopulation_1,
    ImmersedBoundaryCellPopulation_2,
    ImmersedBoundaryCellPopulation_3,
    ImmersedBoundaryKinematicFeedbackForce_1,
    ImmersedBoundaryKinematicFeedbackForce_2,
    ImmersedBoundaryKinematicFeedbackForce_3,
    ImmersedBoundaryLinearDifferentialAdhesionForce_1,
    ImmersedBoundaryLinearDifferentialAdhesionForce_2,
    ImmersedBoundaryLinearDifferentialAdhesionForce_3,
    ImmersedBoundaryLinearInteractionForce_1,
    ImmersedBoundaryLinearInteractionForce_2,
    ImmersedBoundaryLinearInteractionForce_3,
    ImmersedBoundaryLinearMembraneForce_1,
    ImmersedBoundaryLinearMembraneForce_2,
    ImmersedBoundaryLinearMembraneForce_3,
    ImmersedBoundaryMorseInteractionForce_1,
    ImmersedBoundaryMorseInteractionForce_2,
    ImmersedBoundaryMorseInteractionForce_3,
    ImmersedBoundaryMorseMembraneForce_1,
    ImmersedBoundaryMorseMembraneForce_2,
    ImmersedBoundaryMorseMembraneForce_3,
    ImmersedBoundaryNeighbourNumberWriter_1_1,
    ImmersedBoundaryNeighbourNumberWriter_1_2,
    ImmersedBoundaryNeighbourNumberWriter_1_3,
    ImmersedBoundaryNeighbourNumberWriter_2_2,
    ImmersedBoundaryNeighbourNumberWriter_2_3,
    ImmersedBoundaryNeighbourNumberWriter_3_3,
    ImmersedBoundarySimulationModifier_1,
    ImmersedBoundarySimulationModifier_2,
    ImmersedBoundarySimulationModifier_3,
    ImmersedBoundarySvgWriter_1,
    ImmersedBoundarySvgWriter_2,
    ImmersedBoundarySvgWriter_3,
    IsolatedLabelledCellKiller_1,
    IsolatedLabelledCellKiller_2,
    IsolatedLabelledCellKiller_3,
    LabelDependentBernoulliTrialCellCycleModel,
    LegacyCellProliferativeTypesWriter_1_1,
    LegacyCellProliferativeTypesWriter_1_2,
    LegacyCellProliferativeTypesWriter_1_3,
    LegacyCellProliferativeTypesWriter_2_2,
    LegacyCellProliferativeTypesWriter_2_3,
    LegacyCellProliferativeTypesWriter_3_3,
    MeshBasedCellPopulation_1_1,
    MeshBasedCellPopulation_1_2,
    MeshBasedCellPopulation_1_3,
    MeshBasedCellPopulation_2_2,
    MeshBasedCellPopulation_2_3,
    MeshBasedCellPopulation_3_3,
    MeshBasedCellPopulationWithGhostNodes_1,
    MeshBasedCellPopulationWithGhostNodes_2,
    MeshBasedCellPopulationWithGhostNodes_3,
    NagaiHondaDifferentialAdhesionForce_1,
    NagaiHondaDifferentialAdhesionForce_2,
    NagaiHondaDifferentialAdhesionForce_3,
    NagaiHondaForce_1,
    NagaiHondaForce_2,
    NagaiHondaForce_3,
    NoCellCycleModel,
    NodeBasedCellPopulation_1,
    NodeBasedCellPopulation_2,
    NodeBasedCellPopulation_3,
    NodeBasedCellPopulationWithBuskeUpdate_1,
    NodeBasedCellPopulationWithBuskeUpdate_2,
    NodeBasedCellPopulationWithBuskeUpdate_3,
    NodeBasedCellPopulationWithParticles_1,
    NodeBasedCellPopulationWithParticles_2,
    NodeBasedCellPopulationWithParticles_3,
    NodeLocationWriter_1_1,
    NodeLocationWriter_1_2,
    NodeLocationWriter_1_3,
    NodeLocationWriter_2_2,
    NodeLocationWriter_2_3,
    NodeLocationWriter_3_3,
    NodeVelocityWriter_1_1,
    NodeVelocityWriter_1_2,
    NodeVelocityWriter_1_3,
    NodeVelocityWriter_2_2,
    NodeVelocityWriter_2_3,
    NodeVelocityWriter_3_3,
    NormallyDistributedTargetAreaModifier_1,
    NormallyDistributedTargetAreaModifier_2,
    NormallyDistributedTargetAreaModifier_3,
    NullSrnModel,
    OffLatticeSimulation_1_1,
    OffLatticeSimulation_1_2,
    OffLatticeSimulation_1_3,
    OffLatticeSimulation_2_2,
    OffLatticeSimulation_2_3,
    OffLatticeSimulation_3_3,
    OnLatticeSimulation_1,
    OnLatticeSimulation_2,
    OnLatticeSimulation_3,
    PlanarPolarisedFarhadifarForce_1,
    PlanarPolarisedFarhadifarForce_2,
    PlanarPolarisedFarhadifarForce_3,
    PlaneBasedCellKiller_1,
    PlaneBasedCellKiller_2,
    PlaneBasedCellKiller_3,
    PlaneBoundaryCondition_1_1,
    PlaneBoundaryCondition_1_2,
    PlaneBoundaryCondition_1_3,
    PlaneBoundaryCondition_2_2,
    PlaneBoundaryCondition_2_3,
    PlaneBoundaryCondition_3_3,
    PottsBasedCellPopulation_1,
    PottsBasedCellPopulation_2,
    PottsBasedCellPopulation_3,
    PythonSimulationModifier_2,
    PythonSimulationModifier_3,
    RadialCellDataDistributionWriter_1_1,
    RadialCellDataDistributionWriter_1_2,
    RadialCellDataDistributionWriter_1_3,
    RadialCellDataDistributionWriter_2_2,
    RadialCellDataDistributionWriter_2_3,
    RadialCellDataDistributionWriter_3_3,
    RandomCaSwitchingUpdateRule_1,
    RandomCaSwitchingUpdateRule_2,
    RandomCaSwitchingUpdateRule_3,
    RandomCellKiller_1,
    RandomCellKiller_2,
    RandomCellKiller_3,
    RandomDirectionCentreBasedDivisionRule_1_1,
    RandomDirectionCentreBasedDivisionRule_1_2,
    RandomDirectionCentreBasedDivisionRule_1_3,
    RandomDirectionCentreBasedDivisionRule_2_2,
    RandomDirectionCentreBasedDivisionRule_2_3,
    RandomDirectionCentreBasedDivisionRule_3_3,
    RandomDirectionVertexBasedDivisionRule_1,
    RandomDirectionVertexBasedDivisionRule_2,
    RandomDirectionVertexBasedDivisionRule_3,
    RepulsionForce_1,
    RepulsionForce_2,
    RepulsionForce_3,
    ShortAxisImmersedBoundaryDivisionRule_1,
    ShortAxisImmersedBoundaryDivisionRule_2,
    ShortAxisImmersedBoundaryDivisionRule_3,
    ShortAxisVertexBasedDivisionRule_1,
    ShortAxisVertexBasedDivisionRule_2,
    ShortAxisVertexBasedDivisionRule_3,
    ShovingCaBasedDivisionRule_1,
    ShovingCaBasedDivisionRule_2,
    ShovingCaBasedDivisionRule_3,
    SimpleOxygenBasedCellCycleModel,
    SimpleTargetAreaModifier_1,
    SimpleTargetAreaModifier_2,
    SimpleTargetAreaModifier_3,
    SimulationTime,
    SlidingBoundaryCondition_1,
    SlidingBoundaryCondition_2,
    SlidingBoundaryCondition_3,
    SphereGeometryBoundaryCondition_1,
    SphereGeometryBoundaryCondition_2,
    SphereGeometryBoundaryCondition_3,
    StemCellProliferativeType,
    StochasticOxygenBasedCellCycleModel,
    SurfaceAreaConstraintPottsUpdateRule_1,
    SurfaceAreaConstraintPottsUpdateRule_2,
    SurfaceAreaConstraintPottsUpdateRule_3,
    T2SwapCellKiller_1,
    T2SwapCellKiller_2,
    T2SwapCellKiller_3,
    TargetAreaLinearGrowthModifier_1,
    TargetAreaLinearGrowthModifier_2,
    TargetAreaLinearGrowthModifier_3,
    TargetedCellKiller_1,
    TargetedCellKiller_2,
    TargetedCellKiller_3,
    TransitCellProliferativeType,
    TysonNovakCellCycleModel,
    UniformCellCycleModel,
    UniformG1GenerationalCellCycleModel,
    VertexBasedCellPopulation_1,
    VertexBasedCellPopulation_2,
    VertexBasedCellPopulation_3,
    VertexBasedPopulationSrn_1,
    VertexBasedPopulationSrn_2,
    VertexBasedPopulationSrn_3,
    VertexIntersectionSwapLocationsWriter_1_1,
    VertexIntersectionSwapLocationsWriter_1_2,
    VertexIntersectionSwapLocationsWriter_1_3,
    VertexIntersectionSwapLocationsWriter_2_2,
    VertexIntersectionSwapLocationsWriter_2_3,
    VertexIntersectionSwapLocationsWriter_3_3,
    VertexT1SwapLocationsWriter_1_1,
    VertexT1SwapLocationsWriter_1_2,
    VertexT1SwapLocationsWriter_1_3,
    VertexT1SwapLocationsWriter_2_2,
    VertexT1SwapLocationsWriter_2_3,
    VertexT1SwapLocationsWriter_3_3,
    VertexT2SwapLocationsWriter_1_1,
    VertexT2SwapLocationsWriter_1_2,
    VertexT2SwapLocationsWriter_1_3,
    VertexT2SwapLocationsWriter_2_2,
    VertexT2SwapLocationsWriter_2_3,
    VertexT2SwapLocationsWriter_3_3,
    VertexT3SwapLocationsWriter_1_1,
    VertexT3SwapLocationsWriter_1_2,
    VertexT3SwapLocationsWriter_1_3,
    VertexT3SwapLocationsWriter_2_2,
    VertexT3SwapLocationsWriter_2_3,
    VertexT3SwapLocationsWriter_3_3,
    VolumeConstraintPottsUpdateRule_1,
    VolumeConstraintPottsUpdateRule_2,
    VolumeConstraintPottsUpdateRule_3,
    VolumeTrackingModifier_1,
    VolumeTrackingModifier_2,
    VolumeTrackingModifier_3,
    VonMisesVertexBasedDivisionRule_1,
    VonMisesVertexBasedDivisionRule_2,
    VonMisesVertexBasedDivisionRule_3,
    VoronoiDataWriter_1_1,
    VoronoiDataWriter_1_2,
    VoronoiDataWriter_1_3,
    VoronoiDataWriter_2_2,
    VoronoiDataWriter_2_3,
    VoronoiDataWriter_3_3,
    VtkSceneModifier_2,
    VtkSceneModifier_3,
    WelikyOsterForce_1,
    WelikyOsterForce_2,
    WelikyOsterForce_3,
    WildTypeCellMutationState,
)
from chaste._syntax import TemplateClass, TemplateMethod
from chaste.cell_based._fortests import (
    AbstractCellBasedTestSuite,
    AbstractCellBasedWithTimingsTestSuite,
    SetupNotebookTest,
    TearDownNotebookTest,
)


# Template Class Syntax
class AdhesionPottsUpdateRule(TemplateClass):
    _instantiations = {
        ("1",): AdhesionPottsUpdateRule_1,
        ("2",): AdhesionPottsUpdateRule_2,
        ("3",): AdhesionPottsUpdateRule_3,
    }


class ApoptoticCellKiller(TemplateClass):
    _instantiations = {
        ("1",): ApoptoticCellKiller_1,
        ("2",): ApoptoticCellKiller_2,
        ("3",): ApoptoticCellKiller_3,
    }


class AttractingPlaneBoundaryCondition(TemplateClass):
    _instantiations = {
        ("1",): AttractingPlaneBoundaryCondition_1_1,
        ("1", "1"): AttractingPlaneBoundaryCondition_1_1,
        ("1", "2"): AttractingPlaneBoundaryCondition_1_2,
        ("1", "3"): AttractingPlaneBoundaryCondition_1_3,
        ("2",): AttractingPlaneBoundaryCondition_2_2,
        ("2", "2"): AttractingPlaneBoundaryCondition_2_2,
        ("2", "3"): AttractingPlaneBoundaryCondition_2_3,
        ("3",): AttractingPlaneBoundaryCondition_3_3,
        ("3", "3"): AttractingPlaneBoundaryCondition_3_3,
    }


class BoundaryNodeWriter(TemplateClass):
    _instantiations = {
        ("1",): BoundaryNodeWriter_1_1,
        ("1", "1"): BoundaryNodeWriter_1_1,
        ("1", "2"): BoundaryNodeWriter_1_2,
        ("1", "3"): BoundaryNodeWriter_1_3,
        ("2",): BoundaryNodeWriter_2_2,
        ("2", "2"): BoundaryNodeWriter_2_2,
        ("2", "3"): BoundaryNodeWriter_2_3,
        ("3",): BoundaryNodeWriter_3_3,
        ("3", "3"): BoundaryNodeWriter_3_3,
    }


class BuskeAdhesiveForce(TemplateClass):
    _instantiations = {
        ("1",): BuskeAdhesiveForce_1,
        ("2",): BuskeAdhesiveForce_2,
        ("3",): BuskeAdhesiveForce_3,
    }


class BuskeCompressionForce(TemplateClass):
    _instantiations = {
        ("1",): BuskeCompressionForce_1,
        ("2",): BuskeCompressionForce_2,
        ("3",): BuskeCompressionForce_3,
    }


class BuskeElasticForce(TemplateClass):
    _instantiations = {
        ("1",): BuskeElasticForce_1,
        ("2",): BuskeElasticForce_2,
        ("3",): BuskeElasticForce_3,
    }


class CaBasedCellPopulation(TemplateClass):
    _instantiations = {
        ("1",): CaBasedCellPopulation_1,
        ("2",): CaBasedCellPopulation_2,
        ("3",): CaBasedCellPopulation_3,
    }


class CellAgesWriter(TemplateClass):
    _instantiations = {
        ("1",): CellAgesWriter_1_1,
        ("1", "1"): CellAgesWriter_1_1,
        ("1", "2"): CellAgesWriter_1_2,
        ("1", "3"): CellAgesWriter_1_3,
        ("2",): CellAgesWriter_2_2,
        ("2", "2"): CellAgesWriter_2_2,
        ("2", "3"): CellAgesWriter_2_3,
        ("3",): CellAgesWriter_3_3,
        ("3", "3"): CellAgesWriter_3_3,
    }


class CellAncestorWriter(TemplateClass):
    _instantiations = {
        ("1",): CellAncestorWriter_1_1,
        ("1", "1"): CellAncestorWriter_1_1,
        ("1", "2"): CellAncestorWriter_1_2,
        ("1", "3"): CellAncestorWriter_1_3,
        ("2",): CellAncestorWriter_2_2,
        ("2", "2"): CellAncestorWriter_2_2,
        ("2", "3"): CellAncestorWriter_2_3,
        ("3",): CellAncestorWriter_3_3,
        ("3", "3"): CellAncestorWriter_3_3,
    }


class CellAppliedForceWriter(TemplateClass):
    _instantiations = {
        ("1",): CellAppliedForceWriter_1_1,
        ("1", "1"): CellAppliedForceWriter_1_1,
        ("1", "2"): CellAppliedForceWriter_1_2,
        ("1", "3"): CellAppliedForceWriter_1_3,
        ("2",): CellAppliedForceWriter_2_2,
        ("2", "2"): CellAppliedForceWriter_2_2,
        ("2", "3"): CellAppliedForceWriter_2_3,
        ("3",): CellAppliedForceWriter_3_3,
        ("3", "3"): CellAppliedForceWriter_3_3,
    }


class CellCycleModelProteinConcentrationsWriter(TemplateClass):
    _instantiations = {
        ("1",): CellCycleModelProteinConcentrationsWriter_1_1,
        ("1", "1"): CellCycleModelProteinConcentrationsWriter_1_1,
        ("1", "2"): CellCycleModelProteinConcentrationsWriter_1_2,
        ("1", "3"): CellCycleModelProteinConcentrationsWriter_1_3,
        ("2",): CellCycleModelProteinConcentrationsWriter_2_2,
        ("2", "2"): CellCycleModelProteinConcentrationsWriter_2_2,
        ("2", "3"): CellCycleModelProteinConcentrationsWriter_2_3,
        ("3",): CellCycleModelProteinConcentrationsWriter_3_3,
        ("3", "3"): CellCycleModelProteinConcentrationsWriter_3_3,
    }


class CellDataItemWriter(TemplateClass):
    _instantiations = {
        ("1",): CellDataItemWriter_1_1,
        ("1", "1"): CellDataItemWriter_1_1,
        ("1", "2"): CellDataItemWriter_1_2,
        ("1", "3"): CellDataItemWriter_1_3,
        ("2",): CellDataItemWriter_2_2,
        ("2", "2"): CellDataItemWriter_2_2,
        ("2", "3"): CellDataItemWriter_2_3,
        ("3",): CellDataItemWriter_3_3,
        ("3", "3"): CellDataItemWriter_3_3,
    }


class CellDeltaNotchWriter(TemplateClass):
    _instantiations = {
        ("1",): CellDeltaNotchWriter_1_1,
        ("1", "1"): CellDeltaNotchWriter_1_1,
        ("1", "2"): CellDeltaNotchWriter_1_2,
        ("1", "3"): CellDeltaNotchWriter_1_3,
        ("2",): CellDeltaNotchWriter_2_2,
        ("2", "2"): CellDeltaNotchWriter_2_2,
        ("2", "3"): CellDeltaNotchWriter_2_3,
        ("3",): CellDeltaNotchWriter_3_3,
        ("3", "3"): CellDeltaNotchWriter_3_3,
    }


class CellDivisionLocationsWriter(TemplateClass):
    _instantiations = {
        ("1",): CellDivisionLocationsWriter_1_1,
        ("1", "1"): CellDivisionLocationsWriter_1_1,
        ("1", "2"): CellDivisionLocationsWriter_1_2,
        ("1", "3"): CellDivisionLocationsWriter_1_3,
        ("2",): CellDivisionLocationsWriter_2_2,
        ("2", "2"): CellDivisionLocationsWriter_2_2,
        ("2", "3"): CellDivisionLocationsWriter_2_3,
        ("3",): CellDivisionLocationsWriter_3_3,
        ("3", "3"): CellDivisionLocationsWriter_3_3,
    }


class CellIdWriter(TemplateClass):
    _instantiations = {
        ("1",): CellIdWriter_1_1,
        ("1", "1"): CellIdWriter_1_1,
        ("1", "2"): CellIdWriter_1_2,
        ("1", "3"): CellIdWriter_1_3,
        ("2",): CellIdWriter_2_2,
        ("2", "2"): CellIdWriter_2_2,
        ("2", "3"): CellIdWriter_2_3,
        ("3",): CellIdWriter_3_3,
        ("3", "3"): CellIdWriter_3_3,
    }


class CellLabelWriter(TemplateClass):
    _instantiations = {
        ("1",): CellLabelWriter_1_1,
        ("1", "1"): CellLabelWriter_1_1,
        ("1", "2"): CellLabelWriter_1_2,
        ("1", "3"): CellLabelWriter_1_3,
        ("2",): CellLabelWriter_2_2,
        ("2", "2"): CellLabelWriter_2_2,
        ("2", "3"): CellLabelWriter_2_3,
        ("3",): CellLabelWriter_3_3,
        ("3", "3"): CellLabelWriter_3_3,
    }


class CellLocationIndexWriter(TemplateClass):
    _instantiations = {
        ("1",): CellLocationIndexWriter_1_1,
        ("1", "1"): CellLocationIndexWriter_1_1,
        ("1", "2"): CellLocationIndexWriter_1_2,
        ("1", "3"): CellLocationIndexWriter_1_3,
        ("2",): CellLocationIndexWriter_2_2,
        ("2", "2"): CellLocationIndexWriter_2_2,
        ("2", "3"): CellLocationIndexWriter_2_3,
        ("3",): CellLocationIndexWriter_3_3,
        ("3", "3"): CellLocationIndexWriter_3_3,
    }


class CellMutationStatesCountWriter(TemplateClass):
    _instantiations = {
        ("1",): CellMutationStatesCountWriter_1_1,
        ("1", "1"): CellMutationStatesCountWriter_1_1,
        ("1", "2"): CellMutationStatesCountWriter_1_2,
        ("1", "3"): CellMutationStatesCountWriter_1_3,
        ("2",): CellMutationStatesCountWriter_2_2,
        ("2", "2"): CellMutationStatesCountWriter_2_2,
        ("2", "3"): CellMutationStatesCountWriter_2_3,
        ("3",): CellMutationStatesCountWriter_3_3,
        ("3", "3"): CellMutationStatesCountWriter_3_3,
    }


class CellMutationStatesWriter(TemplateClass):
    _instantiations = {
        ("1",): CellMutationStatesWriter_1_1,
        ("1", "1"): CellMutationStatesWriter_1_1,
        ("1", "2"): CellMutationStatesWriter_1_2,
        ("1", "3"): CellMutationStatesWriter_1_3,
        ("2",): CellMutationStatesWriter_2_2,
        ("2", "2"): CellMutationStatesWriter_2_2,
        ("2", "3"): CellMutationStatesWriter_2_3,
        ("3",): CellMutationStatesWriter_3_3,
        ("3", "3"): CellMutationStatesWriter_3_3,
    }


class CellPopulationAdjacencyMatrixWriter(TemplateClass):
    _instantiations = {
        ("1",): CellPopulationAdjacencyMatrixWriter_1_1,
        ("1", "1"): CellPopulationAdjacencyMatrixWriter_1_1,
        ("1", "2"): CellPopulationAdjacencyMatrixWriter_1_2,
        ("1", "3"): CellPopulationAdjacencyMatrixWriter_1_3,
        ("2",): CellPopulationAdjacencyMatrixWriter_2_2,
        ("2", "2"): CellPopulationAdjacencyMatrixWriter_2_2,
        ("2", "3"): CellPopulationAdjacencyMatrixWriter_2_3,
        ("3",): CellPopulationAdjacencyMatrixWriter_3_3,
        ("3", "3"): CellPopulationAdjacencyMatrixWriter_3_3,
    }


class CellPopulationAreaWriter(TemplateClass):
    _instantiations = {
        ("1",): CellPopulationAreaWriter_1_1,
        ("1", "1"): CellPopulationAreaWriter_1_1,
        ("1", "2"): CellPopulationAreaWriter_1_2,
        ("1", "3"): CellPopulationAreaWriter_1_3,
        ("2",): CellPopulationAreaWriter_2_2,
        ("2", "2"): CellPopulationAreaWriter_2_2,
        ("2", "3"): CellPopulationAreaWriter_2_3,
        ("3",): CellPopulationAreaWriter_3_3,
        ("3", "3"): CellPopulationAreaWriter_3_3,
    }


class CellPopulationElementWriter(TemplateClass):
    _instantiations = {
        ("1",): CellPopulationElementWriter_1_1,
        ("1", "1"): CellPopulationElementWriter_1_1,
        ("1", "2"): CellPopulationElementWriter_1_2,
        ("1", "3"): CellPopulationElementWriter_1_3,
        ("2",): CellPopulationElementWriter_2_2,
        ("2", "2"): CellPopulationElementWriter_2_2,
        ("2", "3"): CellPopulationElementWriter_2_3,
        ("3",): CellPopulationElementWriter_3_3,
        ("3", "3"): CellPopulationElementWriter_3_3,
    }


class CellProliferativePhasesCountWriter(TemplateClass):
    _instantiations = {
        ("1",): CellProliferativePhasesCountWriter_1_1,
        ("1", "1"): CellProliferativePhasesCountWriter_1_1,
        ("1", "2"): CellProliferativePhasesCountWriter_1_2,
        ("1", "3"): CellProliferativePhasesCountWriter_1_3,
        ("2",): CellProliferativePhasesCountWriter_2_2,
        ("2", "2"): CellProliferativePhasesCountWriter_2_2,
        ("2", "3"): CellProliferativePhasesCountWriter_2_3,
        ("3",): CellProliferativePhasesCountWriter_3_3,
        ("3", "3"): CellProliferativePhasesCountWriter_3_3,
    }


class CellProliferativePhasesWriter(TemplateClass):
    _instantiations = {
        ("1",): CellProliferativePhasesWriter_1_1,
        ("1", "1"): CellProliferativePhasesWriter_1_1,
        ("1", "2"): CellProliferativePhasesWriter_1_2,
        ("1", "3"): CellProliferativePhasesWriter_1_3,
        ("2",): CellProliferativePhasesWriter_2_2,
        ("2", "2"): CellProliferativePhasesWriter_2_2,
        ("2", "3"): CellProliferativePhasesWriter_2_3,
        ("3",): CellProliferativePhasesWriter_3_3,
        ("3", "3"): CellProliferativePhasesWriter_3_3,
    }


class CellProliferativeTypesCountWriter(TemplateClass):
    _instantiations = {
        ("1",): CellProliferativeTypesCountWriter_1_1,
        ("1", "1"): CellProliferativeTypesCountWriter_1_1,
        ("1", "2"): CellProliferativeTypesCountWriter_1_2,
        ("1", "3"): CellProliferativeTypesCountWriter_1_3,
        ("2",): CellProliferativeTypesCountWriter_2_2,
        ("2", "2"): CellProliferativeTypesCountWriter_2_2,
        ("2", "3"): CellProliferativeTypesCountWriter_2_3,
        ("3",): CellProliferativeTypesCountWriter_3_3,
        ("3", "3"): CellProliferativeTypesCountWriter_3_3,
    }


class CellProliferativeTypesWriter(TemplateClass):
    _instantiations = {
        ("1",): CellProliferativeTypesWriter_1_1,
        ("1", "1"): CellProliferativeTypesWriter_1_1,
        ("1", "2"): CellProliferativeTypesWriter_1_2,
        ("1", "3"): CellProliferativeTypesWriter_1_3,
        ("2",): CellProliferativeTypesWriter_2_2,
        ("2", "2"): CellProliferativeTypesWriter_2_2,
        ("2", "3"): CellProliferativeTypesWriter_2_3,
        ("3",): CellProliferativeTypesWriter_3_3,
        ("3", "3"): CellProliferativeTypesWriter_3_3,
    }


class CellRadiusWriter(TemplateClass):
    _instantiations = {
        ("1",): CellRadiusWriter_1_1,
        ("1", "1"): CellRadiusWriter_1_1,
        ("1", "2"): CellRadiusWriter_1_2,
        ("1", "3"): CellRadiusWriter_1_3,
        ("2",): CellRadiusWriter_2_2,
        ("2", "2"): CellRadiusWriter_2_2,
        ("2", "3"): CellRadiusWriter_2_3,
        ("3",): CellRadiusWriter_3_3,
        ("3", "3"): CellRadiusWriter_3_3,
    }


class CellRemovalLocationsWriter(TemplateClass):
    _instantiations = {
        ("1",): CellRemovalLocationsWriter_1_1,
        ("1", "1"): CellRemovalLocationsWriter_1_1,
        ("1", "2"): CellRemovalLocationsWriter_1_2,
        ("1", "3"): CellRemovalLocationsWriter_1_3,
        ("2",): CellRemovalLocationsWriter_2_2,
        ("2", "2"): CellRemovalLocationsWriter_2_2,
        ("2", "3"): CellRemovalLocationsWriter_2_3,
        ("3",): CellRemovalLocationsWriter_3_3,
        ("3", "3"): CellRemovalLocationsWriter_3_3,
    }


class CellRosetteRankWriter(TemplateClass):
    _instantiations = {
        ("1",): CellRosetteRankWriter_1_1,
        ("1", "1"): CellRosetteRankWriter_1_1,
        ("1", "2"): CellRosetteRankWriter_1_2,
        ("1", "3"): CellRosetteRankWriter_1_3,
        ("2",): CellRosetteRankWriter_2_2,
        ("2", "2"): CellRosetteRankWriter_2_2,
        ("2", "3"): CellRosetteRankWriter_2_3,
        ("3",): CellRosetteRankWriter_3_3,
        ("3", "3"): CellRosetteRankWriter_3_3,
    }


class CellsGenerator(TemplateClass):
    _instantiations = {
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


class CellVolumesWriter(TemplateClass):
    _instantiations = {
        ("1",): CellVolumesWriter_1_1,
        ("1", "1"): CellVolumesWriter_1_1,
        ("1", "2"): CellVolumesWriter_1_2,
        ("1", "3"): CellVolumesWriter_1_3,
        ("2",): CellVolumesWriter_2_2,
        ("2", "2"): CellVolumesWriter_2_2,
        ("2", "3"): CellVolumesWriter_2_3,
        ("3",): CellVolumesWriter_3_3,
        ("3", "3"): CellVolumesWriter_3_3,
    }


class ChemotacticForce(TemplateClass):
    _instantiations = {
        ("1",): ChemotacticForce_1,
        ("2",): ChemotacticForce_2,
        ("3",): ChemotacticForce_3,
    }


class ChemotaxisPottsUpdateRule(TemplateClass):
    _instantiations = {
        ("1",): ChemotaxisPottsUpdateRule_1,
        ("2",): ChemotaxisPottsUpdateRule_2,
        ("3",): ChemotaxisPottsUpdateRule_3,
    }


class DeltaNotchEdgeInteriorTrackingModifier(TemplateClass):
    _instantiations = {
        ("1",): DeltaNotchEdgeInteriorTrackingModifier_1,
        ("2",): DeltaNotchEdgeInteriorTrackingModifier_2,
        ("3",): DeltaNotchEdgeInteriorTrackingModifier_3,
    }


class DeltaNotchTrackingModifier(TemplateClass):
    _instantiations = {
        ("1",): DeltaNotchTrackingModifier_1,
        ("2",): DeltaNotchTrackingModifier_2,
        ("3",): DeltaNotchTrackingModifier_3,
    }


class DeltaNotchEdgeTrackingModifier(TemplateClass):
    _instantiations = {
        ("1",): DeltaNotchEdgeTrackingModifier_1,
        ("2",): DeltaNotchEdgeTrackingModifier_2,
        ("3",): DeltaNotchEdgeTrackingModifier_3,
    }


class DifferentialAdhesionGeneralisedLinearSpringForce(TemplateClass):
    _instantiations = {
        ("1",): DifferentialAdhesionGeneralisedLinearSpringForce_1_1,
        ("1", "1"): DifferentialAdhesionGeneralisedLinearSpringForce_1_1,
        ("1", "2"): DifferentialAdhesionGeneralisedLinearSpringForce_1_2,
        ("1", "3"): DifferentialAdhesionGeneralisedLinearSpringForce_1_3,
        ("2",): DifferentialAdhesionGeneralisedLinearSpringForce_2_2,
        ("2", "2"): DifferentialAdhesionGeneralisedLinearSpringForce_2_2,
        ("2", "3"): DifferentialAdhesionGeneralisedLinearSpringForce_2_3,
        ("3",): DifferentialAdhesionGeneralisedLinearSpringForce_3_3,
        ("3", "3"): DifferentialAdhesionGeneralisedLinearSpringForce_3_3,
    }


class DifferentialAdhesionPottsUpdateRule(TemplateClass):
    _instantiations = {
        ("1",): DifferentialAdhesionPottsUpdateRule_1,
        ("2",): DifferentialAdhesionPottsUpdateRule_2,
        ("3",): DifferentialAdhesionPottsUpdateRule_3,
    }


class DiffusionCaUpdateRule(TemplateClass):
    _instantiations = {
        ("1",): DiffusionCaUpdateRule_1,
        ("2",): DiffusionCaUpdateRule_2,
        ("3",): DiffusionCaUpdateRule_3,
    }


class DiffusionForce(TemplateClass):
    _instantiations = {
        ("1",): DiffusionForce_1,
        ("2",): DiffusionForce_2,
        ("3",): DiffusionForce_3,
    }


class DivisionBiasTrackingModifier(TemplateClass):
    _instantiations = {
        ("1",): DivisionBiasTrackingModifier_1,
        ("2",): DivisionBiasTrackingModifier_2,
        ("3",): DivisionBiasTrackingModifier_3,
    }


class ExclusionCaBasedDivisionRule(TemplateClass):
    _instantiations = {
        ("1",): ExclusionCaBasedDivisionRule_1,
        ("2",): ExclusionCaBasedDivisionRule_2,
        ("3",): ExclusionCaBasedDivisionRule_3,
    }


class ExtrinsicPullModifier(TemplateClass):
    _instantiations = {
        ("1",): ExtrinsicPullModifier_1,
        ("2",): ExtrinsicPullModifier_2,
        ("3",): ExtrinsicPullModifier_3,
    }


class FarhadifarForce(TemplateClass):
    _instantiations = {
        ("1",): FarhadifarForce_1,
        ("2",): FarhadifarForce_2,
        ("3",): FarhadifarForce_3,
    }


class FixedCentreBasedDivisionRule(TemplateClass):
    _instantiations = {
        ("1",): FixedCentreBasedDivisionRule_1_1,
        ("1", "1"): FixedCentreBasedDivisionRule_1_1,
        ("1", "2"): FixedCentreBasedDivisionRule_1_2,
        ("1", "3"): FixedCentreBasedDivisionRule_1_3,
        ("2",): FixedCentreBasedDivisionRule_2_2,
        ("2", "2"): FixedCentreBasedDivisionRule_2_2,
        ("2", "3"): FixedCentreBasedDivisionRule_2_3,
        ("3",): FixedCentreBasedDivisionRule_3_3,
        ("3", "3"): FixedCentreBasedDivisionRule_3_3,
    }


class FixedVertexBasedDivisionRule(TemplateClass):
    _instantiations = {
        ("1",): FixedVertexBasedDivisionRule_1,
        ("2",): FixedVertexBasedDivisionRule_2,
        ("3",): FixedVertexBasedDivisionRule_3,
    }


class ForwardEulerNumericalMethod(TemplateClass):
    _instantiations = {
        ("1",): ForwardEulerNumericalMethod_1_1,
        ("1", "1"): ForwardEulerNumericalMethod_1_1,
        ("1", "2"): ForwardEulerNumericalMethod_1_2,
        ("1", "3"): ForwardEulerNumericalMethod_1_3,
        ("2",): ForwardEulerNumericalMethod_2_2,
        ("2", "2"): ForwardEulerNumericalMethod_2_2,
        ("2", "3"): ForwardEulerNumericalMethod_2_3,
        ("3",): ForwardEulerNumericalMethod_3_3,
        ("3", "3"): ForwardEulerNumericalMethod_3_3,
    }


class GeneralisedLinearSpringForce(TemplateClass):
    _instantiations = {
        ("1",): GeneralisedLinearSpringForce_1_1,
        ("1", "1"): GeneralisedLinearSpringForce_1_1,
        ("1", "2"): GeneralisedLinearSpringForce_1_2,
        ("1", "3"): GeneralisedLinearSpringForce_1_3,
        ("2",): GeneralisedLinearSpringForce_2_2,
        ("2", "2"): GeneralisedLinearSpringForce_2_2,
        ("2", "3"): GeneralisedLinearSpringForce_2_3,
        ("3",): GeneralisedLinearSpringForce_3_3,
        ("3", "3"): GeneralisedLinearSpringForce_3_3,
    }


class HeterotypicBoundaryLengthWriter(TemplateClass):
    _instantiations = {
        ("1",): HeterotypicBoundaryLengthWriter_1_1,
        ("1", "1"): HeterotypicBoundaryLengthWriter_1_1,
        ("1", "2"): HeterotypicBoundaryLengthWriter_1_2,
        ("1", "3"): HeterotypicBoundaryLengthWriter_1_3,
        ("2",): HeterotypicBoundaryLengthWriter_2_2,
        ("2", "2"): HeterotypicBoundaryLengthWriter_2_2,
        ("2", "3"): HeterotypicBoundaryLengthWriter_2_3,
        ("3",): HeterotypicBoundaryLengthWriter_3_3,
        ("3", "3"): HeterotypicBoundaryLengthWriter_3_3,
    }


class ImmersedBoundaryBoundaryCellWriter(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundaryBoundaryCellWriter_1_1,
        ("1", "1"): ImmersedBoundaryBoundaryCellWriter_1_1,
        ("1", "2"): ImmersedBoundaryBoundaryCellWriter_1_2,
        ("1", "3"): ImmersedBoundaryBoundaryCellWriter_1_3,
        ("2",): ImmersedBoundaryBoundaryCellWriter_2_2,
        ("2", "2"): ImmersedBoundaryBoundaryCellWriter_2_2,
        ("2", "3"): ImmersedBoundaryBoundaryCellWriter_2_3,
        ("3",): ImmersedBoundaryBoundaryCellWriter_3_3,
        ("3", "3"): ImmersedBoundaryBoundaryCellWriter_3_3,
    }


class ImmersedBoundaryCellPopulation(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundaryCellPopulation_1,
        ("2",): ImmersedBoundaryCellPopulation_2,
        ("3",): ImmersedBoundaryCellPopulation_3,
    }


class ImmersedBoundaryKinematicFeedbackForce(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundaryKinematicFeedbackForce_1,
        ("2",): ImmersedBoundaryKinematicFeedbackForce_2,
        ("3",): ImmersedBoundaryKinematicFeedbackForce_3,
    }


class ImmersedBoundaryLinearDifferentialAdhesionForce(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundaryLinearDifferentialAdhesionForce_1,
        ("2",): ImmersedBoundaryLinearDifferentialAdhesionForce_2,
        ("3",): ImmersedBoundaryLinearDifferentialAdhesionForce_3,
    }


class ImmersedBoundaryLinearInteractionForce(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundaryLinearInteractionForce_1,
        ("2",): ImmersedBoundaryLinearInteractionForce_2,
        ("3",): ImmersedBoundaryLinearInteractionForce_3,
    }


class ImmersedBoundaryLinearMembraneForce(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundaryLinearMembraneForce_1,
        ("2",): ImmersedBoundaryLinearMembraneForce_2,
        ("3",): ImmersedBoundaryLinearMembraneForce_3,
    }


class ImmersedBoundaryMorseInteractionForce(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundaryMorseInteractionForce_1,
        ("2",): ImmersedBoundaryMorseInteractionForce_2,
        ("3",): ImmersedBoundaryMorseInteractionForce_3,
    }


class ImmersedBoundaryMorseMembraneForce(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundaryMorseMembraneForce_1,
        ("2",): ImmersedBoundaryMorseMembraneForce_2,
        ("3",): ImmersedBoundaryMorseMembraneForce_3,
    }


class ImmersedBoundaryNeighbourNumberWriter(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundaryNeighbourNumberWriter_1_1,
        ("1", "1"): ImmersedBoundaryNeighbourNumberWriter_1_1,
        ("1", "2"): ImmersedBoundaryNeighbourNumberWriter_1_2,
        ("1", "3"): ImmersedBoundaryNeighbourNumberWriter_1_3,
        ("2",): ImmersedBoundaryNeighbourNumberWriter_2_2,
        ("2", "2"): ImmersedBoundaryNeighbourNumberWriter_2_2,
        ("2", "3"): ImmersedBoundaryNeighbourNumberWriter_2_3,
        ("3",): ImmersedBoundaryNeighbourNumberWriter_3_3,
        ("3", "3"): ImmersedBoundaryNeighbourNumberWriter_3_3,
    }


class ImmersedBoundarySimulationModifier(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundarySimulationModifier_1,
        ("2",): ImmersedBoundarySimulationModifier_2,
        ("3",): ImmersedBoundarySimulationModifier_3,
    }


class ImmersedBoundarySvgWriter(TemplateClass):
    _instantiations = {
        ("1",): ImmersedBoundarySvgWriter_1,
        ("2",): ImmersedBoundarySvgWriter_2,
        ("3",): ImmersedBoundarySvgWriter_3,
    }


class IsolatedLabelledCellKiller(TemplateClass):
    _instantiations = {
        ("1",): IsolatedLabelledCellKiller_1,
        ("2",): IsolatedLabelledCellKiller_2,
        ("3",): IsolatedLabelledCellKiller_3,
    }


class LegacyCellProliferativeTypesWriter(TemplateClass):
    _instantiations = {
        ("1",): LegacyCellProliferativeTypesWriter_1_1,
        ("1", "1"): LegacyCellProliferativeTypesWriter_1_1,
        ("1", "2"): LegacyCellProliferativeTypesWriter_1_2,
        ("1", "3"): LegacyCellProliferativeTypesWriter_1_3,
        ("2",): LegacyCellProliferativeTypesWriter_2_2,
        ("2", "2"): LegacyCellProliferativeTypesWriter_2_2,
        ("2", "3"): LegacyCellProliferativeTypesWriter_2_3,
        ("3",): LegacyCellProliferativeTypesWriter_3_3,
        ("3", "3"): LegacyCellProliferativeTypesWriter_3_3,
    }


class MeshBasedCellPopulation(TemplateClass):
    _instantiations = {
        ("1",): MeshBasedCellPopulation_1_1,
        ("1", "1"): MeshBasedCellPopulation_1_1,
        ("1", "2"): MeshBasedCellPopulation_1_2,
        ("1", "3"): MeshBasedCellPopulation_1_3,
        ("2",): MeshBasedCellPopulation_2_2,
        ("2", "2"): MeshBasedCellPopulation_2_2,
        ("2", "3"): MeshBasedCellPopulation_2_3,
        ("3",): MeshBasedCellPopulation_3_3,
        ("3", "3"): MeshBasedCellPopulation_3_3,
    }


class MeshBasedCellPopulationWithGhostNodes(TemplateClass):
    _instantiations = {
        ("1",): MeshBasedCellPopulationWithGhostNodes_1,
        ("2",): MeshBasedCellPopulationWithGhostNodes_2,
        ("3",): MeshBasedCellPopulationWithGhostNodes_3,
    }


class NagaiHondaDifferentialAdhesionForce(TemplateClass):
    _instantiations = {
        ("1",): NagaiHondaDifferentialAdhesionForce_1,
        ("2",): NagaiHondaDifferentialAdhesionForce_2,
        ("3",): NagaiHondaDifferentialAdhesionForce_3,
    }


class NagaiHondaForce(TemplateClass):
    _instantiations = {
        ("1",): NagaiHondaForce_1,
        ("2",): NagaiHondaForce_2,
        ("3",): NagaiHondaForce_3,
    }


class NodeBasedCellPopulation(TemplateClass):
    _instantiations = {
        ("1",): NodeBasedCellPopulation_1,
        ("2",): NodeBasedCellPopulation_2,
        ("3",): NodeBasedCellPopulation_3,
    }


class NodeBasedCellPopulationWithBuskeUpdate(TemplateClass):
    _instantiations = {
        ("1",): NodeBasedCellPopulationWithBuskeUpdate_1,
        ("2",): NodeBasedCellPopulationWithBuskeUpdate_2,
        ("3",): NodeBasedCellPopulationWithBuskeUpdate_3,
    }


class NodeBasedCellPopulationWithParticles(TemplateClass):
    _instantiations = {
        ("1",): NodeBasedCellPopulationWithParticles_1,
        ("2",): NodeBasedCellPopulationWithParticles_2,
        ("3",): NodeBasedCellPopulationWithParticles_3,
    }


class NodeLocationWriter(TemplateClass):
    _instantiations = {
        ("1",): NodeLocationWriter_1_1,
        ("1", "1"): NodeLocationWriter_1_1,
        ("1", "2"): NodeLocationWriter_1_2,
        ("1", "3"): NodeLocationWriter_1_3,
        ("2",): NodeLocationWriter_2_2,
        ("2", "2"): NodeLocationWriter_2_2,
        ("2", "3"): NodeLocationWriter_2_3,
        ("3",): NodeLocationWriter_3_3,
        ("3", "3"): NodeLocationWriter_3_3,
    }


class NodeVelocityWriter(TemplateClass):
    _instantiations = {
        ("1",): NodeVelocityWriter_1_1,
        ("1", "1"): NodeVelocityWriter_1_1,
        ("1", "2"): NodeVelocityWriter_1_2,
        ("1", "3"): NodeVelocityWriter_1_3,
        ("2",): NodeVelocityWriter_2_2,
        ("2", "2"): NodeVelocityWriter_2_2,
        ("2", "3"): NodeVelocityWriter_2_3,
        ("3",): NodeVelocityWriter_3_3,
        ("3", "3"): NodeVelocityWriter_3_3,
    }


class NormallyDistributedTargetAreaModifier(TemplateClass):
    _instantiations = {
        ("1",): NormallyDistributedTargetAreaModifier_1,
        ("2",): NormallyDistributedTargetAreaModifier_2,
        ("3",): NormallyDistributedTargetAreaModifier_3,
    }


class OffLatticeSimulation(TemplateClass):
    _instantiations = {
        ("1",): OffLatticeSimulation_1_1,
        ("1", "1"): OffLatticeSimulation_1_1,
        ("1", "2"): OffLatticeSimulation_1_2,
        ("1", "3"): OffLatticeSimulation_1_3,
        ("2",): OffLatticeSimulation_2_2,
        ("2", "2"): OffLatticeSimulation_2_2,
        ("2", "3"): OffLatticeSimulation_2_3,
        ("3",): OffLatticeSimulation_3_3,
        ("3", "3"): OffLatticeSimulation_3_3,
    }


class OnLatticeSimulation(TemplateClass):
    _instantiations = {
        ("1",): OnLatticeSimulation_1,
        ("2",): OnLatticeSimulation_2,
        ("3",): OnLatticeSimulation_3,
    }


class PlanarPolarisedFarhadifarForce(TemplateClass):
    _instantiations = {
        ("1",): PlanarPolarisedFarhadifarForce_1,
        ("2",): PlanarPolarisedFarhadifarForce_2,
        ("3",): PlanarPolarisedFarhadifarForce_3,
    }


class PlaneBasedCellKiller(TemplateClass):
    _instantiations = {
        ("1",): PlaneBasedCellKiller_1,
        ("2",): PlaneBasedCellKiller_2,
        ("3",): PlaneBasedCellKiller_3,
    }


class PlaneBoundaryCondition(TemplateClass):
    _instantiations = {
        ("1",): PlaneBoundaryCondition_1_1,
        ("1", "1"): PlaneBoundaryCondition_1_1,
        ("1", "2"): PlaneBoundaryCondition_1_2,
        ("1", "3"): PlaneBoundaryCondition_1_3,
        ("2",): PlaneBoundaryCondition_2_2,
        ("2", "2"): PlaneBoundaryCondition_2_2,
        ("2", "3"): PlaneBoundaryCondition_2_3,
        ("3",): PlaneBoundaryCondition_3_3,
        ("3", "3"): PlaneBoundaryCondition_3_3,
    }


class PottsBasedCellPopulation(TemplateClass):
    _instantiations = {
        ("1",): PottsBasedCellPopulation_1,
        ("2",): PottsBasedCellPopulation_2,
        ("3",): PottsBasedCellPopulation_3,
    }


class PythonSimulationModifier(TemplateClass):
    _instantiations = {
        ("2",): PythonSimulationModifier_2,
        ("3",): PythonSimulationModifier_3,
    }


class RadialCellDataDistributionWriter(TemplateClass):
    _instantiations = {
        ("1",): RadialCellDataDistributionWriter_1_1,
        ("1", "1"): RadialCellDataDistributionWriter_1_1,
        ("1", "2"): RadialCellDataDistributionWriter_1_2,
        ("1", "3"): RadialCellDataDistributionWriter_1_3,
        ("2",): RadialCellDataDistributionWriter_2_2,
        ("2", "2"): RadialCellDataDistributionWriter_2_2,
        ("2", "3"): RadialCellDataDistributionWriter_2_3,
        ("3",): RadialCellDataDistributionWriter_3_3,
        ("3", "3"): RadialCellDataDistributionWriter_3_3,
    }


class RandomCaSwitchingUpdateRule(TemplateClass):
    _instantiations = {
        ("1",): RandomCaSwitchingUpdateRule_1,
        ("2",): RandomCaSwitchingUpdateRule_2,
        ("3",): RandomCaSwitchingUpdateRule_3,
    }


class RandomCellKiller(TemplateClass):
    _instantiations = {
        ("1",): RandomCellKiller_1,
        ("2",): RandomCellKiller_2,
        ("3",): RandomCellKiller_3,
    }


class RandomDirectionCentreBasedDivisionRule(TemplateClass):
    _instantiations = {
        ("1",): RandomDirectionCentreBasedDivisionRule_1_1,
        ("1", "1"): RandomDirectionCentreBasedDivisionRule_1_1,
        ("1", "2"): RandomDirectionCentreBasedDivisionRule_1_2,
        ("1", "3"): RandomDirectionCentreBasedDivisionRule_1_3,
        ("2",): RandomDirectionCentreBasedDivisionRule_2_2,
        ("2", "2"): RandomDirectionCentreBasedDivisionRule_2_2,
        ("2", "3"): RandomDirectionCentreBasedDivisionRule_2_3,
        ("3",): RandomDirectionCentreBasedDivisionRule_3_3,
        ("3", "3"): RandomDirectionCentreBasedDivisionRule_3_3,
    }


class RandomDirectionVertexBasedDivisionRule(TemplateClass):
    _instantiations = {
        ("1",): RandomDirectionVertexBasedDivisionRule_1,
        ("2",): RandomDirectionVertexBasedDivisionRule_2,
        ("3",): RandomDirectionVertexBasedDivisionRule_3,
    }


class RepulsionForce(TemplateClass):
    _instantiations = {
        ("1",): RepulsionForce_1,
        ("2",): RepulsionForce_2,
        ("3",): RepulsionForce_3,
    }


class ShortAxisImmersedBoundaryDivisionRule(TemplateClass):
    _instantiations = {
        ("1",): ShortAxisImmersedBoundaryDivisionRule_1,
        ("2",): ShortAxisImmersedBoundaryDivisionRule_2,
        ("3",): ShortAxisImmersedBoundaryDivisionRule_3,
    }


class ShortAxisVertexBasedDivisionRule(TemplateClass):
    _instantiations = {
        ("1",): ShortAxisVertexBasedDivisionRule_1,
        ("2",): ShortAxisVertexBasedDivisionRule_2,
        ("3",): ShortAxisVertexBasedDivisionRule_3,
    }


class ShovingCaBasedDivisionRule(TemplateClass):
    _instantiations = {
        ("1",): ShovingCaBasedDivisionRule_1,
        ("2",): ShovingCaBasedDivisionRule_2,
        ("3",): ShovingCaBasedDivisionRule_3,
    }


class SimpleTargetAreaModifier(TemplateClass):
    _instantiations = {
        ("1",): SimpleTargetAreaModifier_1,
        ("2",): SimpleTargetAreaModifier_2,
        ("3",): SimpleTargetAreaModifier_3,
    }


class SlidingBoundaryCondition(TemplateClass):
    _instantiations = {
        ("1",): SlidingBoundaryCondition_1,
        ("2",): SlidingBoundaryCondition_2,
        ("3",): SlidingBoundaryCondition_3,
    }


class SphereGeometryBoundaryCondition(TemplateClass):
    _instantiations = {
        ("1",): SphereGeometryBoundaryCondition_1,
        ("2",): SphereGeometryBoundaryCondition_2,
        ("3",): SphereGeometryBoundaryCondition_3,
    }


class SurfaceAreaConstraintPottsUpdateRule(TemplateClass):
    _instantiations = {
        ("1",): SurfaceAreaConstraintPottsUpdateRule_1,
        ("2",): SurfaceAreaConstraintPottsUpdateRule_2,
        ("3",): SurfaceAreaConstraintPottsUpdateRule_3,
    }


class T2SwapCellKiller(TemplateClass):
    _instantiations = {
        ("1",): T2SwapCellKiller_1,
        ("2",): T2SwapCellKiller_2,
        ("3",): T2SwapCellKiller_3,
    }


class TargetAreaLinearGrowthModifier(TemplateClass):
    _instantiations = {
        ("1",): TargetAreaLinearGrowthModifier_1,
        ("2",): TargetAreaLinearGrowthModifier_2,
        ("3",): TargetAreaLinearGrowthModifier_3,
    }


class TargetedCellKiller(TemplateClass):
    _instantiations = {
        ("1",): TargetedCellKiller_1,
        ("2",): TargetedCellKiller_2,
        ("3",): TargetedCellKiller_3,
    }


class VertexBasedCellPopulation(TemplateClass):
    _instantiations = {
        ("1",): VertexBasedCellPopulation_1,
        ("2",): VertexBasedCellPopulation_2,
        ("3",): VertexBasedCellPopulation_3,
    }


class VertexBasedPopulationSrn(TemplateClass):
    _instantiations = {
        ("1",): VertexBasedPopulationSrn_1,
        ("2",): VertexBasedPopulationSrn_2,
        ("3",): VertexBasedPopulationSrn_3,
    }


class VertexIntersectionSwapLocationsWriter(TemplateClass):
    _instantiations = {
        ("1",): VertexIntersectionSwapLocationsWriter_1_1,
        ("1", "1"): VertexIntersectionSwapLocationsWriter_1_1,
        ("1", "2"): VertexIntersectionSwapLocationsWriter_1_2,
        ("1", "3"): VertexIntersectionSwapLocationsWriter_1_3,
        ("2",): VertexIntersectionSwapLocationsWriter_2_2,
        ("2", "2"): VertexIntersectionSwapLocationsWriter_2_2,
        ("2", "3"): VertexIntersectionSwapLocationsWriter_2_3,
        ("3",): VertexIntersectionSwapLocationsWriter_3_3,
        ("3", "3"): VertexIntersectionSwapLocationsWriter_3_3,
    }


class VertexT1SwapLocationsWriter(TemplateClass):
    _instantiations = {
        ("1",): VertexT1SwapLocationsWriter_1_1,
        ("1", "1"): VertexT1SwapLocationsWriter_1_1,
        ("1", "2"): VertexT1SwapLocationsWriter_1_2,
        ("1", "3"): VertexT1SwapLocationsWriter_1_3,
        ("2",): VertexT1SwapLocationsWriter_2_2,
        ("2", "2"): VertexT1SwapLocationsWriter_2_2,
        ("2", "3"): VertexT1SwapLocationsWriter_2_3,
        ("3",): VertexT1SwapLocationsWriter_3_3,
        ("3", "3"): VertexT1SwapLocationsWriter_3_3,
    }


class VertexT2SwapLocationsWriter(TemplateClass):
    _instantiations = {
        ("1",): VertexT2SwapLocationsWriter_1_1,
        ("1", "1"): VertexT2SwapLocationsWriter_1_1,
        ("1", "2"): VertexT2SwapLocationsWriter_1_2,
        ("1", "3"): VertexT2SwapLocationsWriter_1_3,
        ("2",): VertexT2SwapLocationsWriter_2_2,
        ("2", "2"): VertexT2SwapLocationsWriter_2_2,
        ("2", "3"): VertexT2SwapLocationsWriter_2_3,
        ("3",): VertexT2SwapLocationsWriter_3_3,
        ("3", "3"): VertexT2SwapLocationsWriter_3_3,
    }


class VertexT3SwapLocationsWriter(TemplateClass):
    _instantiations = {
        ("1",): VertexT3SwapLocationsWriter_1_1,
        ("1", "1"): VertexT3SwapLocationsWriter_1_1,
        ("1", "2"): VertexT3SwapLocationsWriter_1_2,
        ("1", "3"): VertexT3SwapLocationsWriter_1_3,
        ("2",): VertexT3SwapLocationsWriter_2_2,
        ("2", "2"): VertexT3SwapLocationsWriter_2_2,
        ("2", "3"): VertexT3SwapLocationsWriter_2_3,
        ("3",): VertexT3SwapLocationsWriter_3_3,
        ("3", "3"): VertexT3SwapLocationsWriter_3_3,
    }


class VolumeConstraintPottsUpdateRule(TemplateClass):
    _instantiations = {
        ("1",): VolumeConstraintPottsUpdateRule_1,
        ("2",): VolumeConstraintPottsUpdateRule_2,
        ("3",): VolumeConstraintPottsUpdateRule_3,
    }


class VolumeTrackingModifier(TemplateClass):
    _instantiations = {
        ("1",): VolumeTrackingModifier_1,
        ("2",): VolumeTrackingModifier_2,
        ("3",): VolumeTrackingModifier_3,
    }


class VonMisesVertexBasedDivisionRule(TemplateClass):
    _instantiations = {
        ("1",): VonMisesVertexBasedDivisionRule_1,
        ("2",): VonMisesVertexBasedDivisionRule_2,
        ("3",): VonMisesVertexBasedDivisionRule_3,
    }


class VoronoiDataWriter(TemplateClass):
    _instantiations = {
        ("1",): VoronoiDataWriter_1_1,
        ("1", "1"): VoronoiDataWriter_1_1,
        ("1", "2"): VoronoiDataWriter_1_2,
        ("1", "3"): VoronoiDataWriter_1_3,
        ("2",): VoronoiDataWriter_2_2,
        ("2", "2"): VoronoiDataWriter_2_2,
        ("2", "3"): VoronoiDataWriter_2_3,
        ("3",): VoronoiDataWriter_3_3,
        ("3", "3"): VoronoiDataWriter_3_3,
    }


class VtkSceneModifier(TemplateClass):
    _instantiations = {
        ("2",): VtkSceneModifier_2,
        ("3",): VtkSceneModifier_3,
    }


class WelikyOsterForce(TemplateClass):
    _instantiations = {
        ("1",): WelikyOsterForce_1,
        ("2",): WelikyOsterForce_2,
        ("3",): WelikyOsterForce_3,
    }


# Cell populations expose their templated writer-adders (AddCellWriter<Writer>,
# AddPopulationWriter<Writer>, ...) as per-writer mangled bindings, generated by
# the PopulationWriter custom generator. Attach a TemplateMethod descriptor for
# each family so the C++-style subscript form works, e.g.
# population.AddCellWriter[CellVolumesWriter]() in place of the mangled
# population.AddCellWriterCellVolumesWriter(). The plain non-templated overload
# inherited from AbstractCellPopulation is kept as a fallback, so the instance
# form population.AddCellWriter(writer) still works. A class is recognised as a
# population by carrying at least one such mangled binding.
_WRITER_METHOD_BASES = (
    "AddCellWriter",
    "AddCellPopulationCountWriter",
    "AddCellPopulationEventWriter",
    "AddPopulationWriter",
)


def _enable_writer_subscript_syntax(namespace):
    for obj in list(namespace.values()):
        if not isinstance(obj, type):
            continue
        for base in _WRITER_METHOD_BASES:
            if isinstance(obj.__dict__.get(base), TemplateMethod):
                continue  # already wired
            if any(name.startswith(base + "_") for name in dir(obj)):
                setattr(obj, base, TemplateMethod(base, getattr(obj, base, None)))


_enable_writer_subscript_syntax(globals())
