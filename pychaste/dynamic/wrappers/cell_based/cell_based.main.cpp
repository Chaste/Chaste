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
#include "AbstractCellProperty.cppwg.hpp"
#include "ApoptoticCellProperty.cppwg.hpp"
#include "CellPropertyCollection.cppwg.hpp"
#include "CellData.cppwg.hpp"
#include "CellLabel.cppwg.hpp"
#include "CellAncestor.cppwg.hpp"
#include "CellId.cppwg.hpp"
#include "CellEdgeData.cppwg.hpp"
#include "CellPropertyRegistry.cppwg.hpp"
#include "AbstractCellProliferativeType.cppwg.hpp"
#include "StemCellProliferativeType.cppwg.hpp"
#include "DefaultCellProliferativeType.cppwg.hpp"
#include "TransitCellProliferativeType.cppwg.hpp"
#include "DifferentiatedCellProliferativeType.cppwg.hpp"
#include "AbstractCellMutationState.cppwg.hpp"
#include "ApcOneHitCellMutationState.cppwg.hpp"
#include "ApcTwoHitCellMutationState.cppwg.hpp"
#include "BetaCateninOneHitCellMutationState.cppwg.hpp"
#include "WildTypeCellMutationState.cppwg.hpp"
#include "AbstractCellCycleModel.cppwg.hpp"
#include "AbstractPhaseBasedCellCycleModel.cppwg.hpp"
#include "AbstractSimpleCellCycleModel.cppwg.hpp"
#include "AbstractSimplePhaseBasedCellCycleModel.cppwg.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.cppwg.hpp"
#include "AbstractCellCycleModelOdeSolver.cppwg.hpp"
#include "AbstractOdeBasedCellCycleModel.cppwg.hpp"
#include "AbstractOdeBasedPhaseBasedCellCycleModel.cppwg.hpp"
#include "NoCellCycleModel.cppwg.hpp"
#include "UniformCellCycleModel.cppwg.hpp"
#include "UniformG1GenerationalCellCycleModel.cppwg.hpp"
#include "StochasticOxygenBasedCellCycleModel.cppwg.hpp"
#include "SimpleOxygenBasedCellCycleModel.cppwg.hpp"
#include "BiasedBernoulliTrialCellCycleModel.cppwg.hpp"
#include "LabelDependentBernoulliTrialCellCycleModel.cppwg.hpp"
#include "AlwaysDivideCellCycleModel.cppwg.hpp"
#include "ContactInhibitionCellCycleModel.cppwg.hpp"
#include "GammaG1CellCycleModel.cppwg.hpp"
#include "ExponentialG1GenerationalCellCycleModel.cppwg.hpp"
#include "TysonNovakCellCycleModel.cppwg.hpp"
#include "Alarcon2004OxygenBasedCellCycleModel.cppwg.hpp"
#include "FixedSequenceCellCycleModel.cppwg.hpp"
#include "BernoulliTrialCellCycleModel.cppwg.hpp"
#include "FixedG1GenerationalCellCycleModel.cppwg.hpp"
#include "AbstractSrnModel.cppwg.hpp"
#include "AbstractOdeSrnModel.cppwg.hpp"
#include "NullSrnModel.cppwg.hpp"
#include "CellSrnModel.cppwg.hpp"
#include "DeltaNotchSrnModel.cppwg.hpp"
#include "DeltaNotchEdgeSrnModel.cppwg.hpp"
#include "DeltaNotchInteriorSrnModel.cppwg.hpp"
#include "Goldbeter1991SrnModel.cppwg.hpp"
#include "Cell.cppwg.hpp"
#include "CellsGeneratorNoCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorNoCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorUniformCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorUniformCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorSimpleOxygenBasedCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorSimpleOxygenBasedCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorUniformG1GenerationalCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorUniformG1GenerationalCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorBiasedBernoulliTrialCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorBiasedBernoulliTrialCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorAlwaysDivideCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorAlwaysDivideCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorContactInhibitionCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorContactInhibitionCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorStochasticOxygenBasedCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorStochasticOxygenBasedCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorGammaG1CellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorGammaG1CellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorExponentialG1GenerationalCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorExponentialG1GenerationalCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorTysonNovakCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorTysonNovakCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorAlarcon2004OxygenBasedCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorAlarcon2004OxygenBasedCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorFixedSequenceCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorFixedSequenceCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorBernoulliTrialCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorBernoulliTrialCellCycleModel_3.cppwg.hpp"
#include "CellsGeneratorFixedG1GenerationalCellCycleModel_2.cppwg.hpp"
#include "CellsGeneratorFixedG1GenerationalCellCycleModel_3.cppwg.hpp"
#include "VertexBasedPopulationSrn2.cppwg.hpp"
#include "VertexBasedPopulationSrn3.cppwg.hpp"
#include "AbstractCellBasedWriter2_2.cppwg.hpp"
#include "AbstractCellBasedWriter3_3.cppwg.hpp"
#include "AbstractCellWriter2_2.cppwg.hpp"
#include "AbstractCellWriter3_3.cppwg.hpp"
#include "CellAgesWriter2_2.cppwg.hpp"
#include "CellAgesWriter3_3.cppwg.hpp"
#include "CellAncestorWriter2_2.cppwg.hpp"
#include "CellAncestorWriter3_3.cppwg.hpp"
#include "CellAppliedForceWriter2_2.cppwg.hpp"
#include "CellAppliedForceWriter3_3.cppwg.hpp"
#include "CellCycleModelProteinConcentrationsWriter2_2.cppwg.hpp"
#include "CellCycleModelProteinConcentrationsWriter3_3.cppwg.hpp"
#include "CellDataItemWriter2_2.cppwg.hpp"
#include "CellDataItemWriter3_3.cppwg.hpp"
#include "CellDeltaNotchWriter2_2.cppwg.hpp"
#include "CellDeltaNotchWriter3_3.cppwg.hpp"
#include "CellIdWriter2_2.cppwg.hpp"
#include "CellIdWriter3_3.cppwg.hpp"
#include "CellLabelWriter2_2.cppwg.hpp"
#include "CellLabelWriter3_3.cppwg.hpp"
#include "CellLocationIndexWriter2_2.cppwg.hpp"
#include "CellLocationIndexWriter3_3.cppwg.hpp"
#include "CellMutationStatesWriter2_2.cppwg.hpp"
#include "CellMutationStatesWriter3_3.cppwg.hpp"
#include "CellProliferativePhasesWriter2_2.cppwg.hpp"
#include "CellProliferativePhasesWriter3_3.cppwg.hpp"
#include "CellProliferativeTypesWriter2_2.cppwg.hpp"
#include "CellProliferativeTypesWriter3_3.cppwg.hpp"
#include "CellRadiusWriter2_2.cppwg.hpp"
#include "CellRadiusWriter3_3.cppwg.hpp"
#include "CellRosetteRankWriter2_2.cppwg.hpp"
#include "CellRosetteRankWriter3_3.cppwg.hpp"
#include "CellVolumesWriter2_2.cppwg.hpp"
#include "CellVolumesWriter3_3.cppwg.hpp"
#include "ImmersedBoundaryBoundaryCellWriter2_2.cppwg.hpp"
#include "ImmersedBoundaryBoundaryCellWriter3_3.cppwg.hpp"
#include "ImmersedBoundaryNeighbourNumberWriter2_2.cppwg.hpp"
#include "ImmersedBoundaryNeighbourNumberWriter3_3.cppwg.hpp"
#include "LegacyCellProliferativeTypesWriter2_2.cppwg.hpp"
#include "LegacyCellProliferativeTypesWriter3_3.cppwg.hpp"
#include "AbstractCellPopulationWriter2_2.cppwg.hpp"
#include "AbstractCellPopulationWriter3_3.cppwg.hpp"
#include "BoundaryNodeWriter2_2.cppwg.hpp"
#include "BoundaryNodeWriter3_3.cppwg.hpp"
#include "CellPopulationAdjacencyMatrixWriter2_2.cppwg.hpp"
#include "CellPopulationAdjacencyMatrixWriter3_3.cppwg.hpp"
#include "CellPopulationAreaWriter2_2.cppwg.hpp"
#include "CellPopulationAreaWriter3_3.cppwg.hpp"
#include "CellPopulationElementWriter2_2.cppwg.hpp"
#include "CellPopulationElementWriter3_3.cppwg.hpp"
#include "HeterotypicBoundaryLengthWriter2_2.cppwg.hpp"
#include "HeterotypicBoundaryLengthWriter3_3.cppwg.hpp"
#include "NodeLocationWriter2_2.cppwg.hpp"
#include "NodeLocationWriter3_3.cppwg.hpp"
#include "NodeVelocityWriter2_2.cppwg.hpp"
#include "NodeVelocityWriter3_3.cppwg.hpp"
#include "RadialCellDataDistributionWriter2_2.cppwg.hpp"
#include "RadialCellDataDistributionWriter3_3.cppwg.hpp"
#include "VertexIntersectionSwapLocationsWriter2_2.cppwg.hpp"
#include "VertexIntersectionSwapLocationsWriter3_3.cppwg.hpp"
#include "VertexT1SwapLocationsWriter2_2.cppwg.hpp"
#include "VertexT1SwapLocationsWriter3_3.cppwg.hpp"
#include "VertexT2SwapLocationsWriter2_2.cppwg.hpp"
#include "VertexT2SwapLocationsWriter3_3.cppwg.hpp"
#include "VertexT3SwapLocationsWriter2_2.cppwg.hpp"
#include "VertexT3SwapLocationsWriter3_3.cppwg.hpp"
#include "VoronoiDataWriter2_2.cppwg.hpp"
#include "VoronoiDataWriter3_3.cppwg.hpp"
#include "AbstractCellPopulationCountWriter2_2.cppwg.hpp"
#include "AbstractCellPopulationCountWriter3_3.cppwg.hpp"
#include "CellMutationStatesCountWriter2_2.cppwg.hpp"
#include "CellMutationStatesCountWriter3_3.cppwg.hpp"
#include "CellProliferativePhasesCountWriter2_2.cppwg.hpp"
#include "CellProliferativePhasesCountWriter3_3.cppwg.hpp"
#include "CellProliferativeTypesCountWriter2_2.cppwg.hpp"
#include "CellProliferativeTypesCountWriter3_3.cppwg.hpp"
#include "AbstractCellPopulationEventWriter2_2.cppwg.hpp"
#include "AbstractCellPopulationEventWriter3_3.cppwg.hpp"
#include "CellDivisionLocationsWriter2_2.cppwg.hpp"
#include "CellDivisionLocationsWriter3_3.cppwg.hpp"
#include "CellRemovalLocationsWriter2_2.cppwg.hpp"
#include "CellRemovalLocationsWriter3_3.cppwg.hpp"
#include "AbstractUpdateRule2.cppwg.hpp"
#include "AbstractUpdateRule3.cppwg.hpp"
#include "AbstractCaUpdateRule2.cppwg.hpp"
#include "AbstractCaUpdateRule3.cppwg.hpp"
#include "AbstractPottsUpdateRule2.cppwg.hpp"
#include "AbstractPottsUpdateRule3.cppwg.hpp"
#include "AbstractCaSwitchingUpdateRule2.cppwg.hpp"
#include "AbstractCaSwitchingUpdateRule3.cppwg.hpp"
#include "AdhesionPottsUpdateRule2.cppwg.hpp"
#include "AdhesionPottsUpdateRule3.cppwg.hpp"
#include "ChemotaxisPottsUpdateRule2.cppwg.hpp"
#include "ChemotaxisPottsUpdateRule3.cppwg.hpp"
#include "DifferentialAdhesionPottsUpdateRule2.cppwg.hpp"
#include "DifferentialAdhesionPottsUpdateRule3.cppwg.hpp"
#include "DiffusionCaUpdateRule2.cppwg.hpp"
#include "DiffusionCaUpdateRule3.cppwg.hpp"
#include "RandomCaSwitchingUpdateRule2.cppwg.hpp"
#include "RandomCaSwitchingUpdateRule3.cppwg.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule2.cppwg.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule3.cppwg.hpp"
#include "VolumeConstraintPottsUpdateRule2.cppwg.hpp"
#include "VolumeConstraintPottsUpdateRule3.cppwg.hpp"
#include "AbstractCentreBasedDivisionRule2_2.cppwg.hpp"
#include "AbstractCentreBasedDivisionRule3_3.cppwg.hpp"
#include "AbstractCaBasedDivisionRule2.cppwg.hpp"
#include "AbstractCaBasedDivisionRule3.cppwg.hpp"
#include "AbstractImmersedBoundaryDivisionRule2.cppwg.hpp"
#include "AbstractImmersedBoundaryDivisionRule3.cppwg.hpp"
#include "AbstractVertexBasedDivisionRule2.cppwg.hpp"
#include "AbstractVertexBasedDivisionRule3.cppwg.hpp"
#include "ExclusionCaBasedDivisionRule2.cppwg.hpp"
#include "ExclusionCaBasedDivisionRule3.cppwg.hpp"
#include "FixedCentreBasedDivisionRule2_2.cppwg.hpp"
#include "FixedCentreBasedDivisionRule3_3.cppwg.hpp"
#include "FixedVertexBasedDivisionRule2.cppwg.hpp"
#include "FixedVertexBasedDivisionRule3.cppwg.hpp"
#include "RandomDirectionCentreBasedDivisionRule2_2.cppwg.hpp"
#include "RandomDirectionCentreBasedDivisionRule3_3.cppwg.hpp"
#include "RandomDirectionVertexBasedDivisionRule2.cppwg.hpp"
#include "RandomDirectionVertexBasedDivisionRule3.cppwg.hpp"
#include "ShortAxisImmersedBoundaryDivisionRule2.cppwg.hpp"
#include "ShortAxisImmersedBoundaryDivisionRule3.cppwg.hpp"
#include "ShortAxisVertexBasedDivisionRule2.cppwg.hpp"
#include "ShortAxisVertexBasedDivisionRule3.cppwg.hpp"
#include "ShovingCaBasedDivisionRule2.cppwg.hpp"
#include "ShovingCaBasedDivisionRule3.cppwg.hpp"
#include "VonMisesVertexBasedDivisionRule2.cppwg.hpp"
#include "VonMisesVertexBasedDivisionRule3.cppwg.hpp"
#include "AbstractForce2_2.cppwg.hpp"
#include "AbstractForce3_3.cppwg.hpp"
#include "AbstractTwoBodyInteractionForce2_2.cppwg.hpp"
#include "AbstractTwoBodyInteractionForce3_3.cppwg.hpp"
#include "AbstractImmersedBoundaryForce2.cppwg.hpp"
#include "AbstractImmersedBoundaryForce3.cppwg.hpp"
#include "BuskeAdhesiveForce2.cppwg.hpp"
#include "BuskeAdhesiveForce3.cppwg.hpp"
#include "BuskeCompressionForce2.cppwg.hpp"
#include "BuskeCompressionForce3.cppwg.hpp"
#include "BuskeElasticForce2.cppwg.hpp"
#include "BuskeElasticForce3.cppwg.hpp"
#include "ChemotacticForce2.cppwg.hpp"
#include "ChemotacticForce3.cppwg.hpp"
#include "DiffusionForce2.cppwg.hpp"
#include "DiffusionForce3.cppwg.hpp"
#include "FarhadifarForce2.cppwg.hpp"
#include "FarhadifarForce3.cppwg.hpp"
#include "GeneralisedLinearSpringForce2_2.cppwg.hpp"
#include "GeneralisedLinearSpringForce3_3.cppwg.hpp"
#include "ImmersedBoundaryKinematicFeedbackForce2.cppwg.hpp"
#include "ImmersedBoundaryKinematicFeedbackForce3.cppwg.hpp"
#include "ImmersedBoundaryLinearDifferentialAdhesionForce2.cppwg.hpp"
#include "ImmersedBoundaryLinearDifferentialAdhesionForce3.cppwg.hpp"
#include "ImmersedBoundaryLinearInteractionForce2.cppwg.hpp"
#include "ImmersedBoundaryLinearInteractionForce3.cppwg.hpp"
#include "ImmersedBoundaryLinearMembraneForce2.cppwg.hpp"
#include "ImmersedBoundaryLinearMembraneForce3.cppwg.hpp"
#include "ImmersedBoundaryMorseInteractionForce2.cppwg.hpp"
#include "ImmersedBoundaryMorseInteractionForce3.cppwg.hpp"
#include "ImmersedBoundaryMorseMembraneForce2.cppwg.hpp"
#include "ImmersedBoundaryMorseMembraneForce3.cppwg.hpp"
#include "NagaiHondaForce2.cppwg.hpp"
#include "NagaiHondaForce3.cppwg.hpp"
#include "RepulsionForce2.cppwg.hpp"
#include "RepulsionForce3.cppwg.hpp"
#include "WelikyOsterForce2.cppwg.hpp"
#include "WelikyOsterForce3.cppwg.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce2_2.cppwg.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce3_3.cppwg.hpp"
#include "NagaiHondaDifferentialAdhesionForce2.cppwg.hpp"
#include "NagaiHondaDifferentialAdhesionForce3.cppwg.hpp"
#include "PlanarPolarisedFarhadifarForce2.cppwg.hpp"
#include "PlanarPolarisedFarhadifarForce3.cppwg.hpp"
#include "AbstractCellKiller2.cppwg.hpp"
#include "AbstractCellKiller3.cppwg.hpp"
#include "ApoptoticCellKiller2.cppwg.hpp"
#include "ApoptoticCellKiller3.cppwg.hpp"
#include "IsolatedLabelledCellKiller2.cppwg.hpp"
#include "IsolatedLabelledCellKiller3.cppwg.hpp"
#include "PlaneBasedCellKiller2.cppwg.hpp"
#include "PlaneBasedCellKiller3.cppwg.hpp"
#include "RandomCellKiller2.cppwg.hpp"
#include "RandomCellKiller3.cppwg.hpp"
#include "T2SwapCellKiller2.cppwg.hpp"
#include "T2SwapCellKiller3.cppwg.hpp"
#include "TargetedCellKiller2.cppwg.hpp"
#include "TargetedCellKiller3.cppwg.hpp"
#include "AbstractCellPopulationBoundaryCondition2_2.cppwg.hpp"
#include "AbstractCellPopulationBoundaryCondition3_3.cppwg.hpp"
#include "PlaneBoundaryCondition2_2.cppwg.hpp"
#include "PlaneBoundaryCondition3_3.cppwg.hpp"
#include "SlidingBoundaryCondition2.cppwg.hpp"
#include "SlidingBoundaryCondition3.cppwg.hpp"
#include "SphereGeometryBoundaryCondition2.cppwg.hpp"
#include "SphereGeometryBoundaryCondition3.cppwg.hpp"
#include "AbstractCellPopulation2_2.cppwg.hpp"
#include "AbstractCellPopulation3_3.cppwg.hpp"
#include "AbstractOffLatticeCellPopulation2_2.cppwg.hpp"
#include "AbstractOffLatticeCellPopulation3_3.cppwg.hpp"
#include "AbstractCentreBasedCellPopulation2_2.cppwg.hpp"
#include "AbstractCentreBasedCellPopulation3_3.cppwg.hpp"
#include "AbstractOnLatticeCellPopulation2.cppwg.hpp"
#include "AbstractOnLatticeCellPopulation3.cppwg.hpp"
#include "CaBasedCellPopulation2.cppwg.hpp"
#include "CaBasedCellPopulation3.cppwg.hpp"
#include "ImmersedBoundaryCellPopulation2.cppwg.hpp"
#include "ImmersedBoundaryCellPopulation3.cppwg.hpp"
#include "MeshBasedCellPopulation2_2.cppwg.hpp"
#include "MeshBasedCellPopulation3_3.cppwg.hpp"
#include "MeshBasedCellPopulationWithGhostNodes2.cppwg.hpp"
#include "MeshBasedCellPopulationWithGhostNodes3.cppwg.hpp"
#include "NodeBasedCellPopulation2.cppwg.hpp"
#include "NodeBasedCellPopulation3.cppwg.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate2.cppwg.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate3.cppwg.hpp"
#include "NodeBasedCellPopulationWithParticles2.cppwg.hpp"
#include "NodeBasedCellPopulationWithParticles3.cppwg.hpp"
#include "VertexBasedCellPopulation2.cppwg.hpp"
#include "VertexBasedCellPopulation3.cppwg.hpp"
#include "PottsBasedCellPopulation2.cppwg.hpp"
#include "PottsBasedCellPopulation3.cppwg.hpp"
#include "AbstractCellBasedSimulationModifier2_2.cppwg.hpp"
#include "AbstractCellBasedSimulationModifier3_3.cppwg.hpp"
#include "AbstractTargetAreaModifier2.cppwg.hpp"
#include "AbstractTargetAreaModifier3.cppwg.hpp"
#include "DeltaNotchTrackingModifier2.cppwg.hpp"
#include "DeltaNotchTrackingModifier3.cppwg.hpp"
#include "DeltaNotchEdgeTrackingModifier2.cppwg.hpp"
#include "DeltaNotchEdgeTrackingModifier3.cppwg.hpp"
#include "DeltaNotchEdgeInteriorTrackingModifier2.cppwg.hpp"
#include "DeltaNotchEdgeInteriorTrackingModifier3.cppwg.hpp"
#include "DivisionBiasTrackingModifier2.cppwg.hpp"
#include "DivisionBiasTrackingModifier3.cppwg.hpp"
#include "ExtrinsicPullModifier2.cppwg.hpp"
#include "ExtrinsicPullModifier3.cppwg.hpp"
#include "ImmersedBoundarySimulationModifier2.cppwg.hpp"
#include "ImmersedBoundarySimulationModifier3.cppwg.hpp"
#include "ImmersedBoundarySvgWriter2.cppwg.hpp"
#include "ImmersedBoundarySvgWriter3.cppwg.hpp"
#include "NormallyDistributedTargetAreaModifier2.cppwg.hpp"
#include "NormallyDistributedTargetAreaModifier3.cppwg.hpp"
#include "SimpleTargetAreaModifier2.cppwg.hpp"
#include "SimpleTargetAreaModifier3.cppwg.hpp"
#include "TargetAreaLinearGrowthModifier2.cppwg.hpp"
#include "TargetAreaLinearGrowthModifier3.cppwg.hpp"
#include "VolumeTrackingModifier2.cppwg.hpp"
#include "VolumeTrackingModifier3.cppwg.hpp"
#include "AbstractNumericalMethod2_2.cppwg.hpp"
#include "AbstractNumericalMethod3_3.cppwg.hpp"
#include "ForwardEulerNumericalMethod2_2.cppwg.hpp"
#include "ForwardEulerNumericalMethod3_3.cppwg.hpp"
#include "SimulationTime.cppwg.hpp"
#include "AbstractCellBasedSimulation2_2.cppwg.hpp"
#include "AbstractCellBasedSimulation3_3.cppwg.hpp"
#include "OffLatticeSimulation2_2.cppwg.hpp"
#include "OffLatticeSimulation3_3.cppwg.hpp"
#include "OnLatticeSimulation2.cppwg.hpp"
#include "OnLatticeSimulation3.cppwg.hpp"
#include "AttractingPlaneBoundaryCondition2_2.cppwg.hpp"
#include "AttractingPlaneBoundaryCondition3_3.cppwg.hpp"
#include "VtkSceneModifier2.cppwg.hpp"
#include "VtkSceneModifier3.cppwg.hpp"
#include "PythonSimulationModifier2.cppwg.hpp"
#include "PythonSimulationModifier3.cppwg.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_pychaste_cell_based, m)
{
    register_AbstractCellProperty_class(m);
    register_ApoptoticCellProperty_class(m);
    register_CellPropertyCollection_class(m);
    register_CellData_class(m);
    register_CellLabel_class(m);
    register_CellAncestor_class(m);
    register_CellId_class(m);
    register_CellEdgeData_class(m);
    register_CellPropertyRegistry_class(m);
    register_AbstractCellProliferativeType_class(m);
    register_StemCellProliferativeType_class(m);
    register_DefaultCellProliferativeType_class(m);
    register_TransitCellProliferativeType_class(m);
    register_DifferentiatedCellProliferativeType_class(m);
    register_AbstractCellMutationState_class(m);
    register_ApcOneHitCellMutationState_class(m);
    register_ApcTwoHitCellMutationState_class(m);
    register_BetaCateninOneHitCellMutationState_class(m);
    register_WildTypeCellMutationState_class(m);
    register_AbstractCellCycleModel_class(m);
    register_AbstractPhaseBasedCellCycleModel_class(m);
    register_AbstractSimpleCellCycleModel_class(m);
    register_AbstractSimplePhaseBasedCellCycleModel_class(m);
    register_AbstractSimpleGenerationalCellCycleModel_class(m);
    register_AbstractCellCycleModelOdeSolver_class(m);
    register_AbstractOdeBasedCellCycleModel_class(m);
    register_AbstractOdeBasedPhaseBasedCellCycleModel_class(m);
    register_NoCellCycleModel_class(m);
    register_UniformCellCycleModel_class(m);
    register_UniformG1GenerationalCellCycleModel_class(m);
    register_StochasticOxygenBasedCellCycleModel_class(m);
    register_SimpleOxygenBasedCellCycleModel_class(m);
    register_BiasedBernoulliTrialCellCycleModel_class(m);
    register_LabelDependentBernoulliTrialCellCycleModel_class(m);
    register_AlwaysDivideCellCycleModel_class(m);
    register_ContactInhibitionCellCycleModel_class(m);
    register_GammaG1CellCycleModel_class(m);
    register_ExponentialG1GenerationalCellCycleModel_class(m);
    register_TysonNovakCellCycleModel_class(m);
    register_Alarcon2004OxygenBasedCellCycleModel_class(m);
    register_FixedSequenceCellCycleModel_class(m);
    register_BernoulliTrialCellCycleModel_class(m);
    register_FixedG1GenerationalCellCycleModel_class(m);
    register_AbstractSrnModel_class(m);
    register_AbstractOdeSrnModel_class(m);
    register_NullSrnModel_class(m);
    register_CellSrnModel_class(m);
    register_DeltaNotchSrnModel_class(m);
    register_DeltaNotchEdgeSrnModel_class(m);
    register_DeltaNotchInteriorSrnModel_class(m);
    register_Goldbeter1991SrnModel_class(m);
    register_Cell_class(m);
    register_CellsGeneratorNoCellCycleModel_2_class(m);
    register_CellsGeneratorNoCellCycleModel_3_class(m);
    register_CellsGeneratorUniformCellCycleModel_2_class(m);
    register_CellsGeneratorUniformCellCycleModel_3_class(m);
    register_CellsGeneratorSimpleOxygenBasedCellCycleModel_2_class(m);
    register_CellsGeneratorSimpleOxygenBasedCellCycleModel_3_class(m);
    register_CellsGeneratorUniformG1GenerationalCellCycleModel_2_class(m);
    register_CellsGeneratorUniformG1GenerationalCellCycleModel_3_class(m);
    register_CellsGeneratorBiasedBernoulliTrialCellCycleModel_2_class(m);
    register_CellsGeneratorBiasedBernoulliTrialCellCycleModel_3_class(m);
    register_CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2_class(m);
    register_CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_3_class(m);
    register_CellsGeneratorAlwaysDivideCellCycleModel_2_class(m);
    register_CellsGeneratorAlwaysDivideCellCycleModel_3_class(m);
    register_CellsGeneratorContactInhibitionCellCycleModel_2_class(m);
    register_CellsGeneratorContactInhibitionCellCycleModel_3_class(m);
    register_CellsGeneratorStochasticOxygenBasedCellCycleModel_2_class(m);
    register_CellsGeneratorStochasticOxygenBasedCellCycleModel_3_class(m);
    register_CellsGeneratorGammaG1CellCycleModel_2_class(m);
    register_CellsGeneratorGammaG1CellCycleModel_3_class(m);
    register_CellsGeneratorExponentialG1GenerationalCellCycleModel_2_class(m);
    register_CellsGeneratorExponentialG1GenerationalCellCycleModel_3_class(m);
    register_CellsGeneratorTysonNovakCellCycleModel_2_class(m);
    register_CellsGeneratorTysonNovakCellCycleModel_3_class(m);
    register_CellsGeneratorAlarcon2004OxygenBasedCellCycleModel_2_class(m);
    register_CellsGeneratorAlarcon2004OxygenBasedCellCycleModel_3_class(m);
    register_CellsGeneratorFixedSequenceCellCycleModel_2_class(m);
    register_CellsGeneratorFixedSequenceCellCycleModel_3_class(m);
    register_CellsGeneratorBernoulliTrialCellCycleModel_2_class(m);
    register_CellsGeneratorBernoulliTrialCellCycleModel_3_class(m);
    register_CellsGeneratorFixedG1GenerationalCellCycleModel_2_class(m);
    register_CellsGeneratorFixedG1GenerationalCellCycleModel_3_class(m);
    register_VertexBasedPopulationSrn2_class(m);
    register_VertexBasedPopulationSrn3_class(m);
    register_AbstractCellBasedWriter2_2_class(m);
    register_AbstractCellBasedWriter3_3_class(m);
    register_AbstractCellWriter2_2_class(m);
    register_AbstractCellWriter3_3_class(m);
    register_CellAgesWriter2_2_class(m);
    register_CellAgesWriter3_3_class(m);
    register_CellAncestorWriter2_2_class(m);
    register_CellAncestorWriter3_3_class(m);
    register_CellAppliedForceWriter2_2_class(m);
    register_CellAppliedForceWriter3_3_class(m);
    register_CellCycleModelProteinConcentrationsWriter2_2_class(m);
    register_CellCycleModelProteinConcentrationsWriter3_3_class(m);
    register_CellDataItemWriter2_2_class(m);
    register_CellDataItemWriter3_3_class(m);
    register_CellDeltaNotchWriter2_2_class(m);
    register_CellDeltaNotchWriter3_3_class(m);
    register_CellIdWriter2_2_class(m);
    register_CellIdWriter3_3_class(m);
    register_CellLabelWriter2_2_class(m);
    register_CellLabelWriter3_3_class(m);
    register_CellLocationIndexWriter2_2_class(m);
    register_CellLocationIndexWriter3_3_class(m);
    register_CellMutationStatesWriter2_2_class(m);
    register_CellMutationStatesWriter3_3_class(m);
    register_CellProliferativePhasesWriter2_2_class(m);
    register_CellProliferativePhasesWriter3_3_class(m);
    register_CellProliferativeTypesWriter2_2_class(m);
    register_CellProliferativeTypesWriter3_3_class(m);
    register_CellRadiusWriter2_2_class(m);
    register_CellRadiusWriter3_3_class(m);
    register_CellRosetteRankWriter2_2_class(m);
    register_CellRosetteRankWriter3_3_class(m);
    register_CellVolumesWriter2_2_class(m);
    register_CellVolumesWriter3_3_class(m);
    register_ImmersedBoundaryBoundaryCellWriter2_2_class(m);
    register_ImmersedBoundaryBoundaryCellWriter3_3_class(m);
    register_ImmersedBoundaryNeighbourNumberWriter2_2_class(m);
    register_ImmersedBoundaryNeighbourNumberWriter3_3_class(m);
    register_LegacyCellProliferativeTypesWriter2_2_class(m);
    register_LegacyCellProliferativeTypesWriter3_3_class(m);
    register_AbstractCellPopulationWriter2_2_class(m);
    register_AbstractCellPopulationWriter3_3_class(m);
    register_BoundaryNodeWriter2_2_class(m);
    register_BoundaryNodeWriter3_3_class(m);
    register_CellPopulationAdjacencyMatrixWriter2_2_class(m);
    register_CellPopulationAdjacencyMatrixWriter3_3_class(m);
    register_CellPopulationAreaWriter2_2_class(m);
    register_CellPopulationAreaWriter3_3_class(m);
    register_CellPopulationElementWriter2_2_class(m);
    register_CellPopulationElementWriter3_3_class(m);
    register_HeterotypicBoundaryLengthWriter2_2_class(m);
    register_HeterotypicBoundaryLengthWriter3_3_class(m);
    register_NodeLocationWriter2_2_class(m);
    register_NodeLocationWriter3_3_class(m);
    register_NodeVelocityWriter2_2_class(m);
    register_NodeVelocityWriter3_3_class(m);
    register_RadialCellDataDistributionWriter2_2_class(m);
    register_RadialCellDataDistributionWriter3_3_class(m);
    register_VertexIntersectionSwapLocationsWriter2_2_class(m);
    register_VertexIntersectionSwapLocationsWriter3_3_class(m);
    register_VertexT1SwapLocationsWriter2_2_class(m);
    register_VertexT1SwapLocationsWriter3_3_class(m);
    register_VertexT2SwapLocationsWriter2_2_class(m);
    register_VertexT2SwapLocationsWriter3_3_class(m);
    register_VertexT3SwapLocationsWriter2_2_class(m);
    register_VertexT3SwapLocationsWriter3_3_class(m);
    register_VoronoiDataWriter2_2_class(m);
    register_VoronoiDataWriter3_3_class(m);
    register_AbstractCellPopulationCountWriter2_2_class(m);
    register_AbstractCellPopulationCountWriter3_3_class(m);
    register_CellMutationStatesCountWriter2_2_class(m);
    register_CellMutationStatesCountWriter3_3_class(m);
    register_CellProliferativePhasesCountWriter2_2_class(m);
    register_CellProliferativePhasesCountWriter3_3_class(m);
    register_CellProliferativeTypesCountWriter2_2_class(m);
    register_CellProliferativeTypesCountWriter3_3_class(m);
    register_AbstractCellPopulationEventWriter2_2_class(m);
    register_AbstractCellPopulationEventWriter3_3_class(m);
    register_CellDivisionLocationsWriter2_2_class(m);
    register_CellDivisionLocationsWriter3_3_class(m);
    register_CellRemovalLocationsWriter2_2_class(m);
    register_CellRemovalLocationsWriter3_3_class(m);
    register_AbstractUpdateRule2_class(m);
    register_AbstractUpdateRule3_class(m);
    register_AbstractCaUpdateRule2_class(m);
    register_AbstractCaUpdateRule3_class(m);
    register_AbstractPottsUpdateRule2_class(m);
    register_AbstractPottsUpdateRule3_class(m);
    register_AbstractCaSwitchingUpdateRule2_class(m);
    register_AbstractCaSwitchingUpdateRule3_class(m);
    register_AdhesionPottsUpdateRule2_class(m);
    register_AdhesionPottsUpdateRule3_class(m);
    register_ChemotaxisPottsUpdateRule2_class(m);
    register_ChemotaxisPottsUpdateRule3_class(m);
    register_DifferentialAdhesionPottsUpdateRule2_class(m);
    register_DifferentialAdhesionPottsUpdateRule3_class(m);
    register_DiffusionCaUpdateRule2_class(m);
    register_DiffusionCaUpdateRule3_class(m);
    register_RandomCaSwitchingUpdateRule2_class(m);
    register_RandomCaSwitchingUpdateRule3_class(m);
    register_SurfaceAreaConstraintPottsUpdateRule2_class(m);
    register_SurfaceAreaConstraintPottsUpdateRule3_class(m);
    register_VolumeConstraintPottsUpdateRule2_class(m);
    register_VolumeConstraintPottsUpdateRule3_class(m);
    register_AbstractCentreBasedDivisionRule2_2_class(m);
    register_AbstractCentreBasedDivisionRule3_3_class(m);
    register_AbstractCaBasedDivisionRule2_class(m);
    register_AbstractCaBasedDivisionRule3_class(m);
    register_AbstractImmersedBoundaryDivisionRule2_class(m);
    register_AbstractImmersedBoundaryDivisionRule3_class(m);
    register_AbstractVertexBasedDivisionRule2_class(m);
    register_AbstractVertexBasedDivisionRule3_class(m);
    register_ExclusionCaBasedDivisionRule2_class(m);
    register_ExclusionCaBasedDivisionRule3_class(m);
    register_FixedCentreBasedDivisionRule2_2_class(m);
    register_FixedCentreBasedDivisionRule3_3_class(m);
    register_FixedVertexBasedDivisionRule2_class(m);
    register_FixedVertexBasedDivisionRule3_class(m);
    register_RandomDirectionCentreBasedDivisionRule2_2_class(m);
    register_RandomDirectionCentreBasedDivisionRule3_3_class(m);
    register_RandomDirectionVertexBasedDivisionRule2_class(m);
    register_RandomDirectionVertexBasedDivisionRule3_class(m);
    register_ShortAxisImmersedBoundaryDivisionRule2_class(m);
    register_ShortAxisImmersedBoundaryDivisionRule3_class(m);
    register_ShortAxisVertexBasedDivisionRule2_class(m);
    register_ShortAxisVertexBasedDivisionRule3_class(m);
    register_ShovingCaBasedDivisionRule2_class(m);
    register_ShovingCaBasedDivisionRule3_class(m);
    register_VonMisesVertexBasedDivisionRule2_class(m);
    register_VonMisesVertexBasedDivisionRule3_class(m);
    register_AbstractForce2_2_class(m);
    register_AbstractForce3_3_class(m);
    register_AbstractTwoBodyInteractionForce2_2_class(m);
    register_AbstractTwoBodyInteractionForce3_3_class(m);
    register_AbstractImmersedBoundaryForce2_class(m);
    register_AbstractImmersedBoundaryForce3_class(m);
    register_BuskeAdhesiveForce2_class(m);
    register_BuskeAdhesiveForce3_class(m);
    register_BuskeCompressionForce2_class(m);
    register_BuskeCompressionForce3_class(m);
    register_BuskeElasticForce2_class(m);
    register_BuskeElasticForce3_class(m);
    register_ChemotacticForce2_class(m);
    register_ChemotacticForce3_class(m);
    register_DiffusionForce2_class(m);
    register_DiffusionForce3_class(m);
    register_FarhadifarForce2_class(m);
    register_FarhadifarForce3_class(m);
    register_GeneralisedLinearSpringForce2_2_class(m);
    register_GeneralisedLinearSpringForce3_3_class(m);
    register_ImmersedBoundaryKinematicFeedbackForce2_class(m);
    register_ImmersedBoundaryKinematicFeedbackForce3_class(m);
    register_ImmersedBoundaryLinearDifferentialAdhesionForce2_class(m);
    register_ImmersedBoundaryLinearDifferentialAdhesionForce3_class(m);
    register_ImmersedBoundaryLinearInteractionForce2_class(m);
    register_ImmersedBoundaryLinearInteractionForce3_class(m);
    register_ImmersedBoundaryLinearMembraneForce2_class(m);
    register_ImmersedBoundaryLinearMembraneForce3_class(m);
    register_ImmersedBoundaryMorseInteractionForce2_class(m);
    register_ImmersedBoundaryMorseInteractionForce3_class(m);
    register_ImmersedBoundaryMorseMembraneForce2_class(m);
    register_ImmersedBoundaryMorseMembraneForce3_class(m);
    register_NagaiHondaForce2_class(m);
    register_NagaiHondaForce3_class(m);
    register_RepulsionForce2_class(m);
    register_RepulsionForce3_class(m);
    register_WelikyOsterForce2_class(m);
    register_WelikyOsterForce3_class(m);
    register_DifferentialAdhesionGeneralisedLinearSpringForce2_2_class(m);
    register_DifferentialAdhesionGeneralisedLinearSpringForce3_3_class(m);
    register_NagaiHondaDifferentialAdhesionForce2_class(m);
    register_NagaiHondaDifferentialAdhesionForce3_class(m);
    register_PlanarPolarisedFarhadifarForce2_class(m);
    register_PlanarPolarisedFarhadifarForce3_class(m);
    register_AbstractCellKiller2_class(m);
    register_AbstractCellKiller3_class(m);
    register_ApoptoticCellKiller2_class(m);
    register_ApoptoticCellKiller3_class(m);
    register_IsolatedLabelledCellKiller2_class(m);
    register_IsolatedLabelledCellKiller3_class(m);
    register_PlaneBasedCellKiller2_class(m);
    register_PlaneBasedCellKiller3_class(m);
    register_RandomCellKiller2_class(m);
    register_RandomCellKiller3_class(m);
    register_T2SwapCellKiller2_class(m);
    register_T2SwapCellKiller3_class(m);
    register_TargetedCellKiller2_class(m);
    register_TargetedCellKiller3_class(m);
    register_AbstractCellPopulationBoundaryCondition2_2_class(m);
    register_AbstractCellPopulationBoundaryCondition3_3_class(m);
    register_PlaneBoundaryCondition2_2_class(m);
    register_PlaneBoundaryCondition3_3_class(m);
    register_SlidingBoundaryCondition2_class(m);
    register_SlidingBoundaryCondition3_class(m);
    register_SphereGeometryBoundaryCondition2_class(m);
    register_SphereGeometryBoundaryCondition3_class(m);
    register_AbstractCellPopulation2_2_class(m);
    register_AbstractCellPopulation3_3_class(m);
    register_AbstractOffLatticeCellPopulation2_2_class(m);
    register_AbstractOffLatticeCellPopulation3_3_class(m);
    register_AbstractCentreBasedCellPopulation2_2_class(m);
    register_AbstractCentreBasedCellPopulation3_3_class(m);
    register_AbstractOnLatticeCellPopulation2_class(m);
    register_AbstractOnLatticeCellPopulation3_class(m);
    register_CaBasedCellPopulation2_class(m);
    register_CaBasedCellPopulation3_class(m);
    register_ImmersedBoundaryCellPopulation2_class(m);
    register_ImmersedBoundaryCellPopulation3_class(m);
    register_MeshBasedCellPopulation2_2_class(m);
    register_MeshBasedCellPopulation3_3_class(m);
    register_MeshBasedCellPopulationWithGhostNodes2_class(m);
    register_MeshBasedCellPopulationWithGhostNodes3_class(m);
    register_NodeBasedCellPopulation2_class(m);
    register_NodeBasedCellPopulation3_class(m);
    register_NodeBasedCellPopulationWithBuskeUpdate2_class(m);
    register_NodeBasedCellPopulationWithBuskeUpdate3_class(m);
    register_NodeBasedCellPopulationWithParticles2_class(m);
    register_NodeBasedCellPopulationWithParticles3_class(m);
    register_VertexBasedCellPopulation2_class(m);
    register_VertexBasedCellPopulation3_class(m);
    register_PottsBasedCellPopulation2_class(m);
    register_PottsBasedCellPopulation3_class(m);
    register_AbstractCellBasedSimulationModifier2_2_class(m);
    register_AbstractCellBasedSimulationModifier3_3_class(m);
    register_AbstractTargetAreaModifier2_class(m);
    register_AbstractTargetAreaModifier3_class(m);
    register_DeltaNotchTrackingModifier2_class(m);
    register_DeltaNotchTrackingModifier3_class(m);
    register_DeltaNotchEdgeTrackingModifier2_class(m);
    register_DeltaNotchEdgeTrackingModifier3_class(m);
    register_DeltaNotchEdgeInteriorTrackingModifier2_class(m);
    register_DeltaNotchEdgeInteriorTrackingModifier3_class(m);
    register_DivisionBiasTrackingModifier2_class(m);
    register_DivisionBiasTrackingModifier3_class(m);
    register_ExtrinsicPullModifier2_class(m);
    register_ExtrinsicPullModifier3_class(m);
    register_ImmersedBoundarySimulationModifier2_class(m);
    register_ImmersedBoundarySimulationModifier3_class(m);
    register_ImmersedBoundarySvgWriter2_class(m);
    register_ImmersedBoundarySvgWriter3_class(m);
    register_NormallyDistributedTargetAreaModifier2_class(m);
    register_NormallyDistributedTargetAreaModifier3_class(m);
    register_SimpleTargetAreaModifier2_class(m);
    register_SimpleTargetAreaModifier3_class(m);
    register_TargetAreaLinearGrowthModifier2_class(m);
    register_TargetAreaLinearGrowthModifier3_class(m);
    register_VolumeTrackingModifier2_class(m);
    register_VolumeTrackingModifier3_class(m);
    register_AbstractNumericalMethod2_2_class(m);
    register_AbstractNumericalMethod3_3_class(m);
    register_ForwardEulerNumericalMethod2_2_class(m);
    register_ForwardEulerNumericalMethod3_3_class(m);
    register_SimulationTime_class(m);
    register_AbstractCellBasedSimulation2_2_class(m);
    register_AbstractCellBasedSimulation3_3_class(m);
    register_OffLatticeSimulation2_2_class(m);
    register_OffLatticeSimulation3_3_class(m);
    register_OnLatticeSimulation2_class(m);
    register_OnLatticeSimulation3_class(m);
    register_AttractingPlaneBoundaryCondition2_2_class(m);
    register_AttractingPlaneBoundaryCondition3_3_class(m);
    register_VtkSceneModifier2_class(m);
    register_VtkSceneModifier3_class(m);
    register_PythonSimulationModifier2_class(m);
    register_PythonSimulationModifier3_class(m);
}
