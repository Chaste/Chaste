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

// This file is auto-generated; manual changes will be overwritten.
// To make enduring changes, see pychaste/dynamic/config.yaml.

#ifndef pychaste_HEADERS_HPP_
#define pychaste_HEADERS_HPP_

// Includes
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "ProgressReporter.hpp"
#include "RandomNumberGenerator.hpp"
#include "TimeStepper.hpp"
#include "Version.hpp"
#include "Identifiable.hpp"
#include "PetscSetupUtils.hpp"
#include "PetscTools.hpp"
#include "ReplicatableVector.hpp"
#include "Timer.hpp"
#include "AbstractOdeSystemInformation.hpp"
#include "AbstractOdeSystem.hpp"
#include "DeltaNotchOdeSystem.hpp"
#include "DeltaNotchEdgeOdeSystem.hpp"
#include "DeltaNotchInteriorOdeSystem.hpp"
#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"
#include "Goldbeter1991OdeSystem.hpp"
#include "TysonNovak2001OdeSystem.hpp"
#include "AbstractPythonOdeSystemInformation.hpp"
#include "AbstractBoundaryCondition.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PdeSimulationTime.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "AbstractNonlinearEllipticPde.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "VolumeDependentAveragedSourceEllipticPde.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "CellBasedEllipticPdeSolver.hpp"
#include "CellBasedParabolicPdeSolver.hpp"
#include "AbstractPdeModifier.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"
#include "AbstractGrowingDomainPdeModifier.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ChastePoint.hpp"
#include "AbstractChasteRegion.hpp"
#include "ChasteCuboid.hpp"
#include "ChasteEllipsoid.hpp"
#include "AbstractElement.hpp"
#include "AbstractMesh.hpp"
#include "NodeAttributes.hpp"
#include "Node.hpp"
#include "Edge.hpp"
#include "EdgeHelper.hpp"
#include "EdgeOperation.hpp"
#include "Element.hpp"
#include "MutableElement.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "MutableMesh.hpp"
#include "NodesOnlyMesh.hpp"
#include "FluidSource.hpp"
#include "ImmersedBoundaryElement.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryHoneycombMeshGenerator.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "Cylindrical2dMesh.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"
#include "PeriodicNodesOnlyMesh.hpp"
#include "Toroidal2dMesh.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "ToroidalHoneycombMeshGenerator.hpp"
#include "VertexMesh.hpp"
#include "MutableVertexMesh.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "VoronoiVertexMeshGenerator.hpp"
#include "PottsElement.hpp"
#include "PottsMesh.hpp"
#include "PottsMeshGenerator.hpp"
#include "PottsMeshWriter.hpp"
#include "Cell.hpp"
#include "CellsGenerator.hpp"
#include "AbstractCellCycleModel.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"
#include "AbstractSimpleCellCycleModel.hpp"
#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "AbstractCellCycleModelOdeSolver.hpp"
#include "AbstractOdeBasedCellCycleModel.hpp"
#include "AbstractOdeBasedPhaseBasedCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "BiasedBernoulliTrialCellCycleModel.hpp"
#include "LabelDependentBernoulliTrialCellCycleModel.hpp"
#include "AlwaysDivideCellCycleModel.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "GammaG1CellCycleModel.hpp"
#include "ExponentialG1GenerationalCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "Alarcon2004OxygenBasedCellCycleModel.hpp"
#include "FixedSequenceCellCycleModel.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellProperty.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellPropertyCollection.hpp"
#include "CellData.hpp"
#include "CellLabel.hpp"
#include "CellAncestor.hpp"
#include "CellId.hpp"
#include "CellEdgeData.hpp"
#include "CellPropertyRegistry.hpp"
#include "SimulationTime.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractOnLatticeCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "SlidingBoundaryCondition.hpp"
#include "SphereGeometryBoundaryCondition.hpp"
#include "AbstractCentreBasedDivisionRule.hpp"
#include "AbstractCaBasedDivisionRule.hpp"
#include "AbstractImmersedBoundaryDivisionRule.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"
#include "ExclusionCaBasedDivisionRule.hpp"
#include "FixedCentreBasedDivisionRule.hpp"
#include "FixedVertexBasedDivisionRule.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"
#include "RandomDirectionVertexBasedDivisionRule.hpp"
#include "ShortAxisImmersedBoundaryDivisionRule.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"
#include "ShovingCaBasedDivisionRule.hpp"
#include "VonMisesVertexBasedDivisionRule.hpp"
#include "AbstractForce.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"
#include "AbstractImmersedBoundaryForce.hpp"
#include "BuskeAdhesiveForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "BuskeElasticForce.hpp"
#include "ChemotacticForce.hpp"
#include "DiffusionForce.hpp"
#include "FarhadifarForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "ImmersedBoundaryKinematicFeedbackForce.hpp"
#include "ImmersedBoundaryLinearDifferentialAdhesionForce.hpp"
#include "ImmersedBoundaryLinearInteractionForce.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"
#include "ImmersedBoundaryMorseInteractionForce.hpp"
#include "ImmersedBoundaryMorseMembraneForce.hpp"
#include "NagaiHondaForce.hpp"
#include "RepulsionForce.hpp"
#include "WelikyOsterForce.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "PlanarPolarisedFarhadifarForce.hpp"
#include "AbstractCellKiller.hpp"
#include "ApoptoticCellKiller.hpp"
#include "IsolatedLabelledCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "RandomCellKiller.hpp"
#include "T2SwapCellKiller.hpp"
#include "TargetedCellKiller.hpp"
#include "VertexBasedPopulationSrn.hpp"
#include "AbstractUpdateRule.hpp"
#include "AbstractCaUpdateRule.hpp"
#include "AbstractPottsUpdateRule.hpp"
#include "AbstractCaSwitchingUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "ChemotaxisPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "RandomCaSwitchingUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AbstractCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellBasedSimulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractTargetAreaModifier.hpp"
#include "DeltaNotchTrackingModifier.hpp"
#include "DeltaNotchEdgeTrackingModifier.hpp"
#include "DeltaNotchEdgeInteriorTrackingModifier.hpp"
#include "DivisionBiasTrackingModifier.hpp"
#include "ExtrinsicPullModifier.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundarySvgWriter.hpp"
#include "NormallyDistributedTargetAreaModifier.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"
#include "VolumeTrackingModifier.hpp"
#include "AbstractNumericalMethod.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "AbstractSrnModel.hpp"
#include "AbstractOdeSrnModel.hpp"
#include "NullSrnModel.hpp"
#include "CellSrnModel.hpp"
#include "DeltaNotchSrnModel.hpp"
#include "DeltaNotchEdgeSrnModel.hpp"
#include "DeltaNotchInteriorSrnModel.hpp"
#include "Goldbeter1991SrnModel.hpp"
#include "AbstractCellBasedWriter.hpp"
#include "AbstractCellWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAppliedForceWriter.hpp"
#include "CellCycleModelProteinConcentrationsWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "CellDeltaNotchWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellRadiusWriter.hpp"
#include "CellRosetteRankWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "ImmersedBoundaryBoundaryCellWriter.hpp"
#include "ImmersedBoundaryNeighbourNumberWriter.hpp"
#include "LegacyCellProliferativeTypesWriter.hpp"
#include "AbstractCellPopulationWriter.hpp"
#include "BoundaryNodeWriter.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"
#include "CellPopulationAreaWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "HeterotypicBoundaryLengthWriter.hpp"
#include "NodeLocationWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "RadialCellDataDistributionWriter.hpp"
#include "VertexIntersectionSwapLocationsWriter.hpp"
#include "VertexT1SwapLocationsWriter.hpp"
#include "VertexT2SwapLocationsWriter.hpp"
#include "VertexT3SwapLocationsWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "AbstractCellPopulationCountWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "AbstractCellPopulationEventWriter.hpp"
#include "CellDivisionLocationsWriter.hpp"
#include "CellRemovalLocationsWriter.hpp"
#include "AttractingPlaneBoundaryCondition.hpp"
#include "VtkSceneModifier.hpp"
#include "PythonSimulationModifier.hpp"
#include "VtkScene.hpp"
#include "AbstractPyChasteActorGenerator.hpp"
#include "CellPopulationPyChasteActorGenerator.hpp"

// Instantiate Template Classes
template class AbstractBoundaryCondition<2>;
template class AbstractBoundaryCondition<3>;
template class ConstBoundaryCondition<2>;
template class ConstBoundaryCondition<3>;
template class AbstractLinearPde<2, 2>;
template class AbstractLinearPde<3, 3>;
template class AbstractLinearParabolicPde<2, 2>;
template class AbstractLinearParabolicPde<3, 3>;
template class AbstractLinearEllipticPde<2, 2>;
template class AbstractLinearEllipticPde<3, 3>;
template class AbstractLinearParabolicPdeSystemForCoupledOdeSystem<2, 2, 1>;
template class AbstractLinearParabolicPdeSystemForCoupledOdeSystem<3, 3, 1>;
template class AbstractNonlinearEllipticPde<2>;
template class AbstractNonlinearEllipticPde<3>;
template class CellwiseSourceEllipticPde<2>;
template class CellwiseSourceEllipticPde<3>;
template class AveragedSourceEllipticPde<2>;
template class AveragedSourceEllipticPde<3>;
template class VolumeDependentAveragedSourceEllipticPde<2>;
template class VolumeDependentAveragedSourceEllipticPde<3>;
template class UniformSourceEllipticPde<2>;
template class UniformSourceEllipticPde<3>;
template class CellwiseSourceParabolicPde<2>;
template class CellwiseSourceParabolicPde<3>;
template class UniformSourceParabolicPde<2>;
template class UniformSourceParabolicPde<3>;
template class AveragedSourceParabolicPde<2>;
template class AveragedSourceParabolicPde<3>;
template class CellBasedEllipticPdeSolver<2>;
template class CellBasedEllipticPdeSolver<3>;
template class CellBasedParabolicPdeSolver<2>;
template class CellBasedParabolicPdeSolver<3>;
template class AbstractPdeModifier<2>;
template class AbstractPdeModifier<3>;
template class AbstractBoxDomainPdeModifier<2>;
template class AbstractBoxDomainPdeModifier<3>;
template class AbstractGrowingDomainPdeModifier<2>;
template class AbstractGrowingDomainPdeModifier<3>;
template class EllipticGrowingDomainPdeModifier<2>;
template class EllipticGrowingDomainPdeModifier<3>;
template class ParabolicGrowingDomainPdeModifier<2>;
template class ParabolicGrowingDomainPdeModifier<3>;
template class EllipticBoxDomainPdeModifier<2>;
template class EllipticBoxDomainPdeModifier<3>;
template class ParabolicBoxDomainPdeModifier<2>;
template class ParabolicBoxDomainPdeModifier<3>;
template class ChastePoint<2>;
template class ChastePoint<3>;
template class AbstractChasteRegion<2>;
template class AbstractChasteRegion<3>;
template class ChasteCuboid<2>;
template class ChasteCuboid<3>;
template class ChasteEllipsoid<2>;
template class ChasteEllipsoid<3>;
template class AbstractElement<1, 2>;
template class AbstractElement<2, 2>;
template class AbstractElement<2, 3>;
template class AbstractElement<3, 3>;
template class AbstractMesh<2, 2>;
template class AbstractMesh<3, 3>;
template class NodeAttributes<2>;
template class NodeAttributes<3>;
template class Node<2>;
template class Node<3>;
template class Edge<2>;
template class Edge<3>;
template class EdgeHelper<2>;
template class EdgeHelper<3>;
template class Element<2, 2>;
template class Element<3, 3>;
template class MutableElement<1, 2>;
template class MutableElement<2, 2>;
template class MutableElement<2, 3>;
template class MutableElement<3, 3>;
template class AbstractTetrahedralMesh<2, 2>;
template class AbstractTetrahedralMesh<3, 3>;
template class TetrahedralMesh<2, 2>;
template class TetrahedralMesh<3, 3>;
template class MutableMesh<2, 2>;
template class MutableMesh<3, 3>;
template class NodesOnlyMesh<2>;
template class NodesOnlyMesh<3>;
template class FluidSource<2>;
template class FluidSource<3>;
template class ImmersedBoundaryElement<1, 2>;
template class ImmersedBoundaryElement<2, 2>;
template class ImmersedBoundaryElement<2, 3>;
template class ImmersedBoundaryElement<3, 3>;
template class ImmersedBoundaryMesh<2, 2>;
template class ImmersedBoundaryMesh<3, 3>;
template class PeriodicNodesOnlyMesh<2>;
template class PeriodicNodesOnlyMesh<3>;
template class VertexMesh<2, 2>;
template class VertexMesh<3, 3>;
template class MutableVertexMesh<2, 2>;
template class MutableVertexMesh<3, 3>;
template class PottsElement<2>;
template class PottsElement<3>;
template class PottsMesh<2>;
template class PottsMesh<3>;
template class PottsMeshGenerator<2>;
template class PottsMeshGenerator<3>;
template class PottsMeshWriter<2>;
template class PottsMeshWriter<3>;
template class CellsGenerator<Alarcon2004OxygenBasedCellCycleModel, 2>;
template class CellsGenerator<Alarcon2004OxygenBasedCellCycleModel, 3>;
template class CellsGenerator<AlwaysDivideCellCycleModel, 2>;
template class CellsGenerator<AlwaysDivideCellCycleModel, 3>;
template class CellsGenerator<BernoulliTrialCellCycleModel, 2>;
template class CellsGenerator<BernoulliTrialCellCycleModel, 3>;
template class CellsGenerator<BiasedBernoulliTrialCellCycleModel, 2>;
template class CellsGenerator<BiasedBernoulliTrialCellCycleModel, 3>;
template class CellsGenerator<ContactInhibitionCellCycleModel, 2>;
template class CellsGenerator<ContactInhibitionCellCycleModel, 3>;
template class CellsGenerator<ExponentialG1GenerationalCellCycleModel, 2>;
template class CellsGenerator<ExponentialG1GenerationalCellCycleModel, 3>;
template class CellsGenerator<FixedG1GenerationalCellCycleModel, 2>;
template class CellsGenerator<FixedG1GenerationalCellCycleModel, 3>;
template class CellsGenerator<FixedSequenceCellCycleModel, 2>;
template class CellsGenerator<FixedSequenceCellCycleModel, 3>;
template class CellsGenerator<GammaG1CellCycleModel, 2>;
template class CellsGenerator<GammaG1CellCycleModel, 3>;
template class CellsGenerator<LabelDependentBernoulliTrialCellCycleModel, 2>;
template class CellsGenerator<LabelDependentBernoulliTrialCellCycleModel, 3>;
template class CellsGenerator<NoCellCycleModel, 2>;
template class CellsGenerator<NoCellCycleModel, 3>;
template class CellsGenerator<SimpleOxygenBasedCellCycleModel, 2>;
template class CellsGenerator<SimpleOxygenBasedCellCycleModel, 3>;
template class CellsGenerator<StochasticOxygenBasedCellCycleModel, 2>;
template class CellsGenerator<StochasticOxygenBasedCellCycleModel, 3>;
template class CellsGenerator<TysonNovakCellCycleModel, 2>;
template class CellsGenerator<TysonNovakCellCycleModel, 3>;
template class CellsGenerator<UniformCellCycleModel, 2>;
template class CellsGenerator<UniformCellCycleModel, 3>;
template class CellsGenerator<UniformG1GenerationalCellCycleModel, 2>;
template class CellsGenerator<UniformG1GenerationalCellCycleModel, 3>;
template class AbstractCellPopulation<2, 2>;
template class AbstractCellPopulation<3, 3>;
template class AbstractOffLatticeCellPopulation<2, 2>;
template class AbstractOffLatticeCellPopulation<3, 3>;
template class AbstractCentreBasedCellPopulation<2, 2>;
template class AbstractCentreBasedCellPopulation<3, 3>;
template class AbstractOnLatticeCellPopulation<2>;
template class AbstractOnLatticeCellPopulation<3>;
template class CaBasedCellPopulation<2>;
template class CaBasedCellPopulation<3>;
template class ImmersedBoundaryCellPopulation<2>;
template class ImmersedBoundaryCellPopulation<3>;
template class MeshBasedCellPopulation<2, 2>;
template class MeshBasedCellPopulation<3, 3>;
template class MeshBasedCellPopulationWithGhostNodes<2>;
template class MeshBasedCellPopulationWithGhostNodes<3>;
template class NodeBasedCellPopulation<2>;
template class NodeBasedCellPopulation<3>;
template class NodeBasedCellPopulationWithBuskeUpdate<2>;
template class NodeBasedCellPopulationWithBuskeUpdate<3>;
template class NodeBasedCellPopulationWithParticles<2>;
template class NodeBasedCellPopulationWithParticles<3>;
template class VertexBasedCellPopulation<2>;
template class VertexBasedCellPopulation<3>;
template class PottsBasedCellPopulation<2>;
template class PottsBasedCellPopulation<3>;
template class AbstractCellPopulationBoundaryCondition<2, 2>;
template class AbstractCellPopulationBoundaryCondition<3, 3>;
template class PlaneBoundaryCondition<2, 2>;
template class PlaneBoundaryCondition<3, 3>;
template class SlidingBoundaryCondition<2>;
template class SlidingBoundaryCondition<3>;
template class SphereGeometryBoundaryCondition<2>;
template class SphereGeometryBoundaryCondition<3>;
template class AbstractCentreBasedDivisionRule<2, 2>;
template class AbstractCentreBasedDivisionRule<3, 3>;
template class AbstractCaBasedDivisionRule<2>;
template class AbstractCaBasedDivisionRule<3>;
template class AbstractImmersedBoundaryDivisionRule<2>;
template class AbstractImmersedBoundaryDivisionRule<3>;
template class AbstractVertexBasedDivisionRule<2>;
template class AbstractVertexBasedDivisionRule<3>;
template class ExclusionCaBasedDivisionRule<2>;
template class ExclusionCaBasedDivisionRule<3>;
template class FixedCentreBasedDivisionRule<2, 2>;
template class FixedCentreBasedDivisionRule<3, 3>;
template class FixedVertexBasedDivisionRule<2>;
template class FixedVertexBasedDivisionRule<3>;
template class RandomDirectionCentreBasedDivisionRule<2, 2>;
template class RandomDirectionCentreBasedDivisionRule<3, 3>;
template class RandomDirectionVertexBasedDivisionRule<2>;
template class RandomDirectionVertexBasedDivisionRule<3>;
template class ShortAxisImmersedBoundaryDivisionRule<2>;
template class ShortAxisImmersedBoundaryDivisionRule<3>;
template class ShortAxisVertexBasedDivisionRule<2>;
template class ShortAxisVertexBasedDivisionRule<3>;
template class ShovingCaBasedDivisionRule<2>;
template class ShovingCaBasedDivisionRule<3>;
template class VonMisesVertexBasedDivisionRule<2>;
template class VonMisesVertexBasedDivisionRule<3>;
template class AbstractForce<2, 2>;
template class AbstractForce<3, 3>;
template class AbstractTwoBodyInteractionForce<2, 2>;
template class AbstractTwoBodyInteractionForce<3, 3>;
template class AbstractImmersedBoundaryForce<2>;
template class AbstractImmersedBoundaryForce<3>;
template class BuskeAdhesiveForce<2>;
template class BuskeAdhesiveForce<3>;
template class BuskeCompressionForce<2>;
template class BuskeCompressionForce<3>;
template class BuskeElasticForce<2>;
template class BuskeElasticForce<3>;
template class ChemotacticForce<2>;
template class ChemotacticForce<3>;
template class DiffusionForce<2>;
template class DiffusionForce<3>;
template class FarhadifarForce<2>;
template class FarhadifarForce<3>;
template class GeneralisedLinearSpringForce<2, 2>;
template class GeneralisedLinearSpringForce<3, 3>;
template class ImmersedBoundaryKinematicFeedbackForce<2>;
template class ImmersedBoundaryKinematicFeedbackForce<3>;
template class ImmersedBoundaryLinearDifferentialAdhesionForce<2>;
template class ImmersedBoundaryLinearDifferentialAdhesionForce<3>;
template class ImmersedBoundaryLinearInteractionForce<2>;
template class ImmersedBoundaryLinearInteractionForce<3>;
template class ImmersedBoundaryLinearMembraneForce<2>;
template class ImmersedBoundaryLinearMembraneForce<3>;
template class ImmersedBoundaryMorseInteractionForce<2>;
template class ImmersedBoundaryMorseInteractionForce<3>;
template class ImmersedBoundaryMorseMembraneForce<2>;
template class ImmersedBoundaryMorseMembraneForce<3>;
template class NagaiHondaForce<2>;
template class NagaiHondaForce<3>;
template class RepulsionForce<2>;
template class RepulsionForce<3>;
template class WelikyOsterForce<2>;
template class WelikyOsterForce<3>;
template class DifferentialAdhesionGeneralisedLinearSpringForce<2, 2>;
template class DifferentialAdhesionGeneralisedLinearSpringForce<3, 3>;
template class NagaiHondaDifferentialAdhesionForce<2>;
template class NagaiHondaDifferentialAdhesionForce<3>;
template class PlanarPolarisedFarhadifarForce<2>;
template class PlanarPolarisedFarhadifarForce<3>;
template class AbstractCellKiller<2>;
template class AbstractCellKiller<3>;
template class ApoptoticCellKiller<2>;
template class ApoptoticCellKiller<3>;
template class IsolatedLabelledCellKiller<2>;
template class IsolatedLabelledCellKiller<3>;
template class PlaneBasedCellKiller<2>;
template class PlaneBasedCellKiller<3>;
template class RandomCellKiller<2>;
template class RandomCellKiller<3>;
template class T2SwapCellKiller<2>;
template class T2SwapCellKiller<3>;
template class TargetedCellKiller<2>;
template class TargetedCellKiller<3>;
template class VertexBasedPopulationSrn<2>;
template class VertexBasedPopulationSrn<3>;
template class AbstractUpdateRule<2>;
template class AbstractUpdateRule<3>;
template class AbstractCaUpdateRule<2>;
template class AbstractCaUpdateRule<3>;
template class AbstractPottsUpdateRule<2>;
template class AbstractPottsUpdateRule<3>;
template class AbstractCaSwitchingUpdateRule<2>;
template class AbstractCaSwitchingUpdateRule<3>;
template class AdhesionPottsUpdateRule<2>;
template class AdhesionPottsUpdateRule<3>;
template class ChemotaxisPottsUpdateRule<2>;
template class ChemotaxisPottsUpdateRule<3>;
template class DifferentialAdhesionPottsUpdateRule<2>;
template class DifferentialAdhesionPottsUpdateRule<3>;
template class DiffusionCaUpdateRule<2>;
template class DiffusionCaUpdateRule<3>;
template class RandomCaSwitchingUpdateRule<2>;
template class RandomCaSwitchingUpdateRule<3>;
template class SurfaceAreaConstraintPottsUpdateRule<2>;
template class SurfaceAreaConstraintPottsUpdateRule<3>;
template class VolumeConstraintPottsUpdateRule<2>;
template class VolumeConstraintPottsUpdateRule<3>;
template class AbstractCellBasedSimulation<2, 2>;
template class AbstractCellBasedSimulation<3, 3>;
template class OffLatticeSimulation<2, 2>;
template class OffLatticeSimulation<3, 3>;
template class OnLatticeSimulation<2>;
template class OnLatticeSimulation<3>;
template class AbstractCellBasedSimulationModifier<2, 2>;
template class AbstractCellBasedSimulationModifier<3, 3>;
template class AbstractTargetAreaModifier<2>;
template class AbstractTargetAreaModifier<3>;
template class DeltaNotchTrackingModifier<2>;
template class DeltaNotchTrackingModifier<3>;
template class DeltaNotchEdgeTrackingModifier<2>;
template class DeltaNotchEdgeTrackingModifier<3>;
template class DeltaNotchEdgeInteriorTrackingModifier<2>;
template class DeltaNotchEdgeInteriorTrackingModifier<3>;
template class DivisionBiasTrackingModifier<2>;
template class DivisionBiasTrackingModifier<3>;
template class ExtrinsicPullModifier<2>;
template class ExtrinsicPullModifier<3>;
template class ImmersedBoundarySimulationModifier<2>;
template class ImmersedBoundarySimulationModifier<3>;
template class ImmersedBoundarySvgWriter<2>;
template class ImmersedBoundarySvgWriter<3>;
template class NormallyDistributedTargetAreaModifier<2>;
template class NormallyDistributedTargetAreaModifier<3>;
template class SimpleTargetAreaModifier<2>;
template class SimpleTargetAreaModifier<3>;
template class TargetAreaLinearGrowthModifier<2>;
template class TargetAreaLinearGrowthModifier<3>;
template class VolumeTrackingModifier<2>;
template class VolumeTrackingModifier<3>;
template class AbstractNumericalMethod<2, 2>;
template class AbstractNumericalMethod<3, 3>;
template class ForwardEulerNumericalMethod<2, 2>;
template class ForwardEulerNumericalMethod<3, 3>;
template class AbstractCellBasedWriter<2, 2>;
template class AbstractCellBasedWriter<3, 3>;
template class AbstractCellWriter<2, 2>;
template class AbstractCellWriter<3, 3>;
template class CellAgesWriter<2, 2>;
template class CellAgesWriter<3, 3>;
template class CellAncestorWriter<2, 2>;
template class CellAncestorWriter<3, 3>;
template class CellAppliedForceWriter<2, 2>;
template class CellAppliedForceWriter<3, 3>;
template class CellCycleModelProteinConcentrationsWriter<2, 2>;
template class CellCycleModelProteinConcentrationsWriter<3, 3>;
template class CellDataItemWriter<2, 2>;
template class CellDataItemWriter<3, 3>;
template class CellDeltaNotchWriter<2, 2>;
template class CellDeltaNotchWriter<3, 3>;
template class CellIdWriter<2, 2>;
template class CellIdWriter<3, 3>;
template class CellLabelWriter<2, 2>;
template class CellLabelWriter<3, 3>;
template class CellLocationIndexWriter<2, 2>;
template class CellLocationIndexWriter<3, 3>;
template class CellMutationStatesWriter<2, 2>;
template class CellMutationStatesWriter<3, 3>;
template class CellProliferativePhasesWriter<2, 2>;
template class CellProliferativePhasesWriter<3, 3>;
template class CellProliferativeTypesWriter<2, 2>;
template class CellProliferativeTypesWriter<3, 3>;
template class CellRadiusWriter<2, 2>;
template class CellRadiusWriter<3, 3>;
template class CellRosetteRankWriter<2, 2>;
template class CellRosetteRankWriter<3, 3>;
template class CellVolumesWriter<2, 2>;
template class CellVolumesWriter<3, 3>;
template class ImmersedBoundaryBoundaryCellWriter<2, 2>;
template class ImmersedBoundaryBoundaryCellWriter<3, 3>;
template class ImmersedBoundaryNeighbourNumberWriter<2, 2>;
template class ImmersedBoundaryNeighbourNumberWriter<3, 3>;
template class LegacyCellProliferativeTypesWriter<2, 2>;
template class LegacyCellProliferativeTypesWriter<3, 3>;
template class AbstractCellPopulationWriter<2, 2>;
template class AbstractCellPopulationWriter<3, 3>;
template class BoundaryNodeWriter<2, 2>;
template class BoundaryNodeWriter<3, 3>;
template class CellPopulationAdjacencyMatrixWriter<2, 2>;
template class CellPopulationAdjacencyMatrixWriter<3, 3>;
template class CellPopulationAreaWriter<2, 2>;
template class CellPopulationAreaWriter<3, 3>;
template class CellPopulationElementWriter<2, 2>;
template class CellPopulationElementWriter<3, 3>;
template class HeterotypicBoundaryLengthWriter<2, 2>;
template class HeterotypicBoundaryLengthWriter<3, 3>;
template class NodeLocationWriter<2, 2>;
template class NodeLocationWriter<3, 3>;
template class NodeVelocityWriter<2, 2>;
template class NodeVelocityWriter<3, 3>;
template class RadialCellDataDistributionWriter<2, 2>;
template class RadialCellDataDistributionWriter<3, 3>;
template class VertexIntersectionSwapLocationsWriter<2, 2>;
template class VertexIntersectionSwapLocationsWriter<3, 3>;
template class VertexT1SwapLocationsWriter<2, 2>;
template class VertexT1SwapLocationsWriter<3, 3>;
template class VertexT2SwapLocationsWriter<2, 2>;
template class VertexT2SwapLocationsWriter<3, 3>;
template class VertexT3SwapLocationsWriter<2, 2>;
template class VertexT3SwapLocationsWriter<3, 3>;
template class VoronoiDataWriter<2, 2>;
template class VoronoiDataWriter<3, 3>;
template class AbstractCellPopulationCountWriter<2, 2>;
template class AbstractCellPopulationCountWriter<3, 3>;
template class CellMutationStatesCountWriter<2, 2>;
template class CellMutationStatesCountWriter<3, 3>;
template class CellProliferativePhasesCountWriter<2, 2>;
template class CellProliferativePhasesCountWriter<3, 3>;
template class CellProliferativeTypesCountWriter<2, 2>;
template class CellProliferativeTypesCountWriter<3, 3>;
template class AbstractCellPopulationEventWriter<2, 2>;
template class AbstractCellPopulationEventWriter<3, 3>;
template class CellDivisionLocationsWriter<2, 2>;
template class CellDivisionLocationsWriter<3, 3>;
template class CellRemovalLocationsWriter<2, 2>;
template class CellRemovalLocationsWriter<3, 3>;
template class AttractingPlaneBoundaryCondition<2, 2>;
template class AttractingPlaneBoundaryCondition<3, 3>;
template class VtkSceneModifier<2>;
template class VtkSceneModifier<3>;
template class PythonSimulationModifier<2>;
template class PythonSimulationModifier<3>;
template class VtkScene<2>;
template class VtkScene<3>;
template class AbstractPyChasteActorGenerator<2>;
template class AbstractPyChasteActorGenerator<3>;
template class CellPopulationPyChasteActorGenerator<2>;
template class CellPopulationPyChasteActorGenerator<3>;

// Typedefs for nicer naming
namespace cppwg
{
    typedef AbstractBoundaryCondition<2> AbstractBoundaryCondition_2;
    typedef AbstractBoundaryCondition<3> AbstractBoundaryCondition_3;
    typedef ConstBoundaryCondition<2> ConstBoundaryCondition_2;
    typedef ConstBoundaryCondition<3> ConstBoundaryCondition_3;
    typedef AbstractLinearPde<2, 2> AbstractLinearPde_2_2;
    typedef AbstractLinearPde<3, 3> AbstractLinearPde_3_3;
    typedef AbstractLinearParabolicPde<2, 2> AbstractLinearParabolicPde_2_2;
    typedef AbstractLinearParabolicPde<3, 3> AbstractLinearParabolicPde_3_3;
    typedef AbstractLinearEllipticPde<2, 2> AbstractLinearEllipticPde_2_2;
    typedef AbstractLinearEllipticPde<3, 3> AbstractLinearEllipticPde_3_3;
    typedef AbstractLinearParabolicPdeSystemForCoupledOdeSystem<2, 2, 1> AbstractLinearParabolicPdeSystemForCoupledOdeSystem_2_2_1;
    typedef AbstractLinearParabolicPdeSystemForCoupledOdeSystem<3, 3, 1> AbstractLinearParabolicPdeSystemForCoupledOdeSystem_3_3_1;
    typedef AbstractNonlinearEllipticPde<2> AbstractNonlinearEllipticPde_2;
    typedef AbstractNonlinearEllipticPde<3> AbstractNonlinearEllipticPde_3;
    typedef CellwiseSourceEllipticPde<2> CellwiseSourceEllipticPde_2;
    typedef CellwiseSourceEllipticPde<3> CellwiseSourceEllipticPde_3;
    typedef AveragedSourceEllipticPde<2> AveragedSourceEllipticPde_2;
    typedef AveragedSourceEllipticPde<3> AveragedSourceEllipticPde_3;
    typedef VolumeDependentAveragedSourceEllipticPde<2> VolumeDependentAveragedSourceEllipticPde_2;
    typedef VolumeDependentAveragedSourceEllipticPde<3> VolumeDependentAveragedSourceEllipticPde_3;
    typedef UniformSourceEllipticPde<2> UniformSourceEllipticPde_2;
    typedef UniformSourceEllipticPde<3> UniformSourceEllipticPde_3;
    typedef CellwiseSourceParabolicPde<2> CellwiseSourceParabolicPde_2;
    typedef CellwiseSourceParabolicPde<3> CellwiseSourceParabolicPde_3;
    typedef UniformSourceParabolicPde<2> UniformSourceParabolicPde_2;
    typedef UniformSourceParabolicPde<3> UniformSourceParabolicPde_3;
    typedef AveragedSourceParabolicPde<2> AveragedSourceParabolicPde_2;
    typedef AveragedSourceParabolicPde<3> AveragedSourceParabolicPde_3;
    typedef CellBasedEllipticPdeSolver<2> CellBasedEllipticPdeSolver_2;
    typedef CellBasedEllipticPdeSolver<3> CellBasedEllipticPdeSolver_3;
    typedef CellBasedParabolicPdeSolver<2> CellBasedParabolicPdeSolver_2;
    typedef CellBasedParabolicPdeSolver<3> CellBasedParabolicPdeSolver_3;
    typedef AbstractPdeModifier<2> AbstractPdeModifier_2;
    typedef AbstractPdeModifier<3> AbstractPdeModifier_3;
    typedef AbstractBoxDomainPdeModifier<2> AbstractBoxDomainPdeModifier_2;
    typedef AbstractBoxDomainPdeModifier<3> AbstractBoxDomainPdeModifier_3;
    typedef AbstractGrowingDomainPdeModifier<2> AbstractGrowingDomainPdeModifier_2;
    typedef AbstractGrowingDomainPdeModifier<3> AbstractGrowingDomainPdeModifier_3;
    typedef EllipticGrowingDomainPdeModifier<2> EllipticGrowingDomainPdeModifier_2;
    typedef EllipticGrowingDomainPdeModifier<3> EllipticGrowingDomainPdeModifier_3;
    typedef ParabolicGrowingDomainPdeModifier<2> ParabolicGrowingDomainPdeModifier_2;
    typedef ParabolicGrowingDomainPdeModifier<3> ParabolicGrowingDomainPdeModifier_3;
    typedef EllipticBoxDomainPdeModifier<2> EllipticBoxDomainPdeModifier_2;
    typedef EllipticBoxDomainPdeModifier<3> EllipticBoxDomainPdeModifier_3;
    typedef ParabolicBoxDomainPdeModifier<2> ParabolicBoxDomainPdeModifier_2;
    typedef ParabolicBoxDomainPdeModifier<3> ParabolicBoxDomainPdeModifier_3;
    typedef ChastePoint<2> ChastePoint_2;
    typedef ChastePoint<3> ChastePoint_3;
    typedef AbstractChasteRegion<2> AbstractChasteRegion_2;
    typedef AbstractChasteRegion<3> AbstractChasteRegion_3;
    typedef ChasteCuboid<2> ChasteCuboid_2;
    typedef ChasteCuboid<3> ChasteCuboid_3;
    typedef ChasteEllipsoid<2> ChasteEllipsoid_2;
    typedef ChasteEllipsoid<3> ChasteEllipsoid_3;
    typedef AbstractElement<1, 2> AbstractElement_1_2;
    typedef AbstractElement<2, 2> AbstractElement_2_2;
    typedef AbstractElement<2, 3> AbstractElement_2_3;
    typedef AbstractElement<3, 3> AbstractElement_3_3;
    typedef AbstractMesh<2, 2> AbstractMesh_2_2;
    typedef AbstractMesh<3, 3> AbstractMesh_3_3;
    typedef NodeAttributes<2> NodeAttributes_2;
    typedef NodeAttributes<3> NodeAttributes_3;
    typedef Node<2> Node_2;
    typedef Node<3> Node_3;
    typedef Edge<2> Edge_2;
    typedef Edge<3> Edge_3;
    typedef EdgeHelper<2> EdgeHelper_2;
    typedef EdgeHelper<3> EdgeHelper_3;
    typedef Element<2, 2> Element_2_2;
    typedef Element<3, 3> Element_3_3;
    typedef MutableElement<1, 2> MutableElement_1_2;
    typedef MutableElement<2, 2> MutableElement_2_2;
    typedef MutableElement<2, 3> MutableElement_2_3;
    typedef MutableElement<3, 3> MutableElement_3_3;
    typedef AbstractTetrahedralMesh<2, 2> AbstractTetrahedralMesh_2_2;
    typedef AbstractTetrahedralMesh<3, 3> AbstractTetrahedralMesh_3_3;
    typedef TetrahedralMesh<2, 2> TetrahedralMesh_2_2;
    typedef TetrahedralMesh<3, 3> TetrahedralMesh_3_3;
    typedef MutableMesh<2, 2> MutableMesh_2_2;
    typedef MutableMesh<3, 3> MutableMesh_3_3;
    typedef NodesOnlyMesh<2> NodesOnlyMesh_2;
    typedef NodesOnlyMesh<3> NodesOnlyMesh_3;
    typedef FluidSource<2> FluidSource_2;
    typedef FluidSource<3> FluidSource_3;
    typedef ImmersedBoundaryElement<1, 2> ImmersedBoundaryElement_1_2;
    typedef ImmersedBoundaryElement<2, 2> ImmersedBoundaryElement_2_2;
    typedef ImmersedBoundaryElement<2, 3> ImmersedBoundaryElement_2_3;
    typedef ImmersedBoundaryElement<3, 3> ImmersedBoundaryElement_3_3;
    typedef ImmersedBoundaryMesh<2, 2> ImmersedBoundaryMesh_2_2;
    typedef ImmersedBoundaryMesh<3, 3> ImmersedBoundaryMesh_3_3;
    typedef PeriodicNodesOnlyMesh<2> PeriodicNodesOnlyMesh_2;
    typedef PeriodicNodesOnlyMesh<3> PeriodicNodesOnlyMesh_3;
    typedef VertexMesh<2, 2> VertexMesh_2_2;
    typedef VertexMesh<3, 3> VertexMesh_3_3;
    typedef MutableVertexMesh<2, 2> MutableVertexMesh_2_2;
    typedef MutableVertexMesh<3, 3> MutableVertexMesh_3_3;
    typedef PottsElement<2> PottsElement_2;
    typedef PottsElement<3> PottsElement_3;
    typedef PottsMesh<2> PottsMesh_2;
    typedef PottsMesh<3> PottsMesh_3;
    typedef PottsMeshGenerator<2> PottsMeshGenerator_2;
    typedef PottsMeshGenerator<3> PottsMeshGenerator_3;
    typedef PottsMeshWriter<2> PottsMeshWriter_2;
    typedef PottsMeshWriter<3> PottsMeshWriter_3;
    typedef CellsGenerator<Alarcon2004OxygenBasedCellCycleModel, 2> CellsGenerator_Alarcon2004OxygenBasedCellCycleModel_2;
    typedef CellsGenerator<Alarcon2004OxygenBasedCellCycleModel, 3> CellsGenerator_Alarcon2004OxygenBasedCellCycleModel_3;
    typedef CellsGenerator<AlwaysDivideCellCycleModel, 2> CellsGenerator_AlwaysDivideCellCycleModel_2;
    typedef CellsGenerator<AlwaysDivideCellCycleModel, 3> CellsGenerator_AlwaysDivideCellCycleModel_3;
    typedef CellsGenerator<BernoulliTrialCellCycleModel, 2> CellsGenerator_BernoulliTrialCellCycleModel_2;
    typedef CellsGenerator<BernoulliTrialCellCycleModel, 3> CellsGenerator_BernoulliTrialCellCycleModel_3;
    typedef CellsGenerator<BiasedBernoulliTrialCellCycleModel, 2> CellsGenerator_BiasedBernoulliTrialCellCycleModel_2;
    typedef CellsGenerator<BiasedBernoulliTrialCellCycleModel, 3> CellsGenerator_BiasedBernoulliTrialCellCycleModel_3;
    typedef CellsGenerator<ContactInhibitionCellCycleModel, 2> CellsGenerator_ContactInhibitionCellCycleModel_2;
    typedef CellsGenerator<ContactInhibitionCellCycleModel, 3> CellsGenerator_ContactInhibitionCellCycleModel_3;
    typedef CellsGenerator<ExponentialG1GenerationalCellCycleModel, 2> CellsGenerator_ExponentialG1GenerationalCellCycleModel_2;
    typedef CellsGenerator<ExponentialG1GenerationalCellCycleModel, 3> CellsGenerator_ExponentialG1GenerationalCellCycleModel_3;
    typedef CellsGenerator<FixedG1GenerationalCellCycleModel, 2> CellsGenerator_FixedG1GenerationalCellCycleModel_2;
    typedef CellsGenerator<FixedG1GenerationalCellCycleModel, 3> CellsGenerator_FixedG1GenerationalCellCycleModel_3;
    typedef CellsGenerator<FixedSequenceCellCycleModel, 2> CellsGenerator_FixedSequenceCellCycleModel_2;
    typedef CellsGenerator<FixedSequenceCellCycleModel, 3> CellsGenerator_FixedSequenceCellCycleModel_3;
    typedef CellsGenerator<GammaG1CellCycleModel, 2> CellsGenerator_GammaG1CellCycleModel_2;
    typedef CellsGenerator<GammaG1CellCycleModel, 3> CellsGenerator_GammaG1CellCycleModel_3;
    typedef CellsGenerator<LabelDependentBernoulliTrialCellCycleModel, 2> CellsGenerator_LabelDependentBernoulliTrialCellCycleModel_2;
    typedef CellsGenerator<LabelDependentBernoulliTrialCellCycleModel, 3> CellsGenerator_LabelDependentBernoulliTrialCellCycleModel_3;
    typedef CellsGenerator<NoCellCycleModel, 2> CellsGenerator_NoCellCycleModel_2;
    typedef CellsGenerator<NoCellCycleModel, 3> CellsGenerator_NoCellCycleModel_3;
    typedef CellsGenerator<SimpleOxygenBasedCellCycleModel, 2> CellsGenerator_SimpleOxygenBasedCellCycleModel_2;
    typedef CellsGenerator<SimpleOxygenBasedCellCycleModel, 3> CellsGenerator_SimpleOxygenBasedCellCycleModel_3;
    typedef CellsGenerator<StochasticOxygenBasedCellCycleModel, 2> CellsGenerator_StochasticOxygenBasedCellCycleModel_2;
    typedef CellsGenerator<StochasticOxygenBasedCellCycleModel, 3> CellsGenerator_StochasticOxygenBasedCellCycleModel_3;
    typedef CellsGenerator<TysonNovakCellCycleModel, 2> CellsGenerator_TysonNovakCellCycleModel_2;
    typedef CellsGenerator<TysonNovakCellCycleModel, 3> CellsGenerator_TysonNovakCellCycleModel_3;
    typedef CellsGenerator<UniformCellCycleModel, 2> CellsGenerator_UniformCellCycleModel_2;
    typedef CellsGenerator<UniformCellCycleModel, 3> CellsGenerator_UniformCellCycleModel_3;
    typedef CellsGenerator<UniformG1GenerationalCellCycleModel, 2> CellsGenerator_UniformG1GenerationalCellCycleModel_2;
    typedef CellsGenerator<UniformG1GenerationalCellCycleModel, 3> CellsGenerator_UniformG1GenerationalCellCycleModel_3;
    typedef AbstractCellPopulation<2, 2> AbstractCellPopulation_2_2;
    typedef AbstractCellPopulation<3, 3> AbstractCellPopulation_3_3;
    typedef AbstractOffLatticeCellPopulation<2, 2> AbstractOffLatticeCellPopulation_2_2;
    typedef AbstractOffLatticeCellPopulation<3, 3> AbstractOffLatticeCellPopulation_3_3;
    typedef AbstractCentreBasedCellPopulation<2, 2> AbstractCentreBasedCellPopulation_2_2;
    typedef AbstractCentreBasedCellPopulation<3, 3> AbstractCentreBasedCellPopulation_3_3;
    typedef AbstractOnLatticeCellPopulation<2> AbstractOnLatticeCellPopulation_2;
    typedef AbstractOnLatticeCellPopulation<3> AbstractOnLatticeCellPopulation_3;
    typedef CaBasedCellPopulation<2> CaBasedCellPopulation_2;
    typedef CaBasedCellPopulation<3> CaBasedCellPopulation_3;
    typedef ImmersedBoundaryCellPopulation<2> ImmersedBoundaryCellPopulation_2;
    typedef ImmersedBoundaryCellPopulation<3> ImmersedBoundaryCellPopulation_3;
    typedef MeshBasedCellPopulation<2, 2> MeshBasedCellPopulation_2_2;
    typedef MeshBasedCellPopulation<3, 3> MeshBasedCellPopulation_3_3;
    typedef MeshBasedCellPopulationWithGhostNodes<2> MeshBasedCellPopulationWithGhostNodes_2;
    typedef MeshBasedCellPopulationWithGhostNodes<3> MeshBasedCellPopulationWithGhostNodes_3;
    typedef NodeBasedCellPopulation<2> NodeBasedCellPopulation_2;
    typedef NodeBasedCellPopulation<3> NodeBasedCellPopulation_3;
    typedef NodeBasedCellPopulationWithBuskeUpdate<2> NodeBasedCellPopulationWithBuskeUpdate_2;
    typedef NodeBasedCellPopulationWithBuskeUpdate<3> NodeBasedCellPopulationWithBuskeUpdate_3;
    typedef NodeBasedCellPopulationWithParticles<2> NodeBasedCellPopulationWithParticles_2;
    typedef NodeBasedCellPopulationWithParticles<3> NodeBasedCellPopulationWithParticles_3;
    typedef VertexBasedCellPopulation<2> VertexBasedCellPopulation_2;
    typedef VertexBasedCellPopulation<3> VertexBasedCellPopulation_3;
    typedef PottsBasedCellPopulation<2> PottsBasedCellPopulation_2;
    typedef PottsBasedCellPopulation<3> PottsBasedCellPopulation_3;
    typedef AbstractCellPopulationBoundaryCondition<2, 2> AbstractCellPopulationBoundaryCondition_2_2;
    typedef AbstractCellPopulationBoundaryCondition<3, 3> AbstractCellPopulationBoundaryCondition_3_3;
    typedef PlaneBoundaryCondition<2, 2> PlaneBoundaryCondition_2_2;
    typedef PlaneBoundaryCondition<3, 3> PlaneBoundaryCondition_3_3;
    typedef SlidingBoundaryCondition<2> SlidingBoundaryCondition_2;
    typedef SlidingBoundaryCondition<3> SlidingBoundaryCondition_3;
    typedef SphereGeometryBoundaryCondition<2> SphereGeometryBoundaryCondition_2;
    typedef SphereGeometryBoundaryCondition<3> SphereGeometryBoundaryCondition_3;
    typedef AbstractCentreBasedDivisionRule<2, 2> AbstractCentreBasedDivisionRule_2_2;
    typedef AbstractCentreBasedDivisionRule<3, 3> AbstractCentreBasedDivisionRule_3_3;
    typedef AbstractCaBasedDivisionRule<2> AbstractCaBasedDivisionRule_2;
    typedef AbstractCaBasedDivisionRule<3> AbstractCaBasedDivisionRule_3;
    typedef AbstractImmersedBoundaryDivisionRule<2> AbstractImmersedBoundaryDivisionRule_2;
    typedef AbstractImmersedBoundaryDivisionRule<3> AbstractImmersedBoundaryDivisionRule_3;
    typedef AbstractVertexBasedDivisionRule<2> AbstractVertexBasedDivisionRule_2;
    typedef AbstractVertexBasedDivisionRule<3> AbstractVertexBasedDivisionRule_3;
    typedef ExclusionCaBasedDivisionRule<2> ExclusionCaBasedDivisionRule_2;
    typedef ExclusionCaBasedDivisionRule<3> ExclusionCaBasedDivisionRule_3;
    typedef FixedCentreBasedDivisionRule<2, 2> FixedCentreBasedDivisionRule_2_2;
    typedef FixedCentreBasedDivisionRule<3, 3> FixedCentreBasedDivisionRule_3_3;
    typedef FixedVertexBasedDivisionRule<2> FixedVertexBasedDivisionRule_2;
    typedef FixedVertexBasedDivisionRule<3> FixedVertexBasedDivisionRule_3;
    typedef RandomDirectionCentreBasedDivisionRule<2, 2> RandomDirectionCentreBasedDivisionRule_2_2;
    typedef RandomDirectionCentreBasedDivisionRule<3, 3> RandomDirectionCentreBasedDivisionRule_3_3;
    typedef RandomDirectionVertexBasedDivisionRule<2> RandomDirectionVertexBasedDivisionRule_2;
    typedef RandomDirectionVertexBasedDivisionRule<3> RandomDirectionVertexBasedDivisionRule_3;
    typedef ShortAxisImmersedBoundaryDivisionRule<2> ShortAxisImmersedBoundaryDivisionRule_2;
    typedef ShortAxisImmersedBoundaryDivisionRule<3> ShortAxisImmersedBoundaryDivisionRule_3;
    typedef ShortAxisVertexBasedDivisionRule<2> ShortAxisVertexBasedDivisionRule_2;
    typedef ShortAxisVertexBasedDivisionRule<3> ShortAxisVertexBasedDivisionRule_3;
    typedef ShovingCaBasedDivisionRule<2> ShovingCaBasedDivisionRule_2;
    typedef ShovingCaBasedDivisionRule<3> ShovingCaBasedDivisionRule_3;
    typedef VonMisesVertexBasedDivisionRule<2> VonMisesVertexBasedDivisionRule_2;
    typedef VonMisesVertexBasedDivisionRule<3> VonMisesVertexBasedDivisionRule_3;
    typedef AbstractForce<2, 2> AbstractForce_2_2;
    typedef AbstractForce<3, 3> AbstractForce_3_3;
    typedef AbstractTwoBodyInteractionForce<2, 2> AbstractTwoBodyInteractionForce_2_2;
    typedef AbstractTwoBodyInteractionForce<3, 3> AbstractTwoBodyInteractionForce_3_3;
    typedef AbstractImmersedBoundaryForce<2> AbstractImmersedBoundaryForce_2;
    typedef AbstractImmersedBoundaryForce<3> AbstractImmersedBoundaryForce_3;
    typedef BuskeAdhesiveForce<2> BuskeAdhesiveForce_2;
    typedef BuskeAdhesiveForce<3> BuskeAdhesiveForce_3;
    typedef BuskeCompressionForce<2> BuskeCompressionForce_2;
    typedef BuskeCompressionForce<3> BuskeCompressionForce_3;
    typedef BuskeElasticForce<2> BuskeElasticForce_2;
    typedef BuskeElasticForce<3> BuskeElasticForce_3;
    typedef ChemotacticForce<2> ChemotacticForce_2;
    typedef ChemotacticForce<3> ChemotacticForce_3;
    typedef DiffusionForce<2> DiffusionForce_2;
    typedef DiffusionForce<3> DiffusionForce_3;
    typedef FarhadifarForce<2> FarhadifarForce_2;
    typedef FarhadifarForce<3> FarhadifarForce_3;
    typedef GeneralisedLinearSpringForce<2, 2> GeneralisedLinearSpringForce_2_2;
    typedef GeneralisedLinearSpringForce<3, 3> GeneralisedLinearSpringForce_3_3;
    typedef ImmersedBoundaryKinematicFeedbackForce<2> ImmersedBoundaryKinematicFeedbackForce_2;
    typedef ImmersedBoundaryKinematicFeedbackForce<3> ImmersedBoundaryKinematicFeedbackForce_3;
    typedef ImmersedBoundaryLinearDifferentialAdhesionForce<2> ImmersedBoundaryLinearDifferentialAdhesionForce_2;
    typedef ImmersedBoundaryLinearDifferentialAdhesionForce<3> ImmersedBoundaryLinearDifferentialAdhesionForce_3;
    typedef ImmersedBoundaryLinearInteractionForce<2> ImmersedBoundaryLinearInteractionForce_2;
    typedef ImmersedBoundaryLinearInteractionForce<3> ImmersedBoundaryLinearInteractionForce_3;
    typedef ImmersedBoundaryLinearMembraneForce<2> ImmersedBoundaryLinearMembraneForce_2;
    typedef ImmersedBoundaryLinearMembraneForce<3> ImmersedBoundaryLinearMembraneForce_3;
    typedef ImmersedBoundaryMorseInteractionForce<2> ImmersedBoundaryMorseInteractionForce_2;
    typedef ImmersedBoundaryMorseInteractionForce<3> ImmersedBoundaryMorseInteractionForce_3;
    typedef ImmersedBoundaryMorseMembraneForce<2> ImmersedBoundaryMorseMembraneForce_2;
    typedef ImmersedBoundaryMorseMembraneForce<3> ImmersedBoundaryMorseMembraneForce_3;
    typedef NagaiHondaForce<2> NagaiHondaForce_2;
    typedef NagaiHondaForce<3> NagaiHondaForce_3;
    typedef RepulsionForce<2> RepulsionForce_2;
    typedef RepulsionForce<3> RepulsionForce_3;
    typedef WelikyOsterForce<2> WelikyOsterForce_2;
    typedef WelikyOsterForce<3> WelikyOsterForce_3;
    typedef DifferentialAdhesionGeneralisedLinearSpringForce<2, 2> DifferentialAdhesionGeneralisedLinearSpringForce_2_2;
    typedef DifferentialAdhesionGeneralisedLinearSpringForce<3, 3> DifferentialAdhesionGeneralisedLinearSpringForce_3_3;
    typedef NagaiHondaDifferentialAdhesionForce<2> NagaiHondaDifferentialAdhesionForce_2;
    typedef NagaiHondaDifferentialAdhesionForce<3> NagaiHondaDifferentialAdhesionForce_3;
    typedef PlanarPolarisedFarhadifarForce<2> PlanarPolarisedFarhadifarForce_2;
    typedef PlanarPolarisedFarhadifarForce<3> PlanarPolarisedFarhadifarForce_3;
    typedef AbstractCellKiller<2> AbstractCellKiller_2;
    typedef AbstractCellKiller<3> AbstractCellKiller_3;
    typedef ApoptoticCellKiller<2> ApoptoticCellKiller_2;
    typedef ApoptoticCellKiller<3> ApoptoticCellKiller_3;
    typedef IsolatedLabelledCellKiller<2> IsolatedLabelledCellKiller_2;
    typedef IsolatedLabelledCellKiller<3> IsolatedLabelledCellKiller_3;
    typedef PlaneBasedCellKiller<2> PlaneBasedCellKiller_2;
    typedef PlaneBasedCellKiller<3> PlaneBasedCellKiller_3;
    typedef RandomCellKiller<2> RandomCellKiller_2;
    typedef RandomCellKiller<3> RandomCellKiller_3;
    typedef T2SwapCellKiller<2> T2SwapCellKiller_2;
    typedef T2SwapCellKiller<3> T2SwapCellKiller_3;
    typedef TargetedCellKiller<2> TargetedCellKiller_2;
    typedef TargetedCellKiller<3> TargetedCellKiller_3;
    typedef VertexBasedPopulationSrn<2> VertexBasedPopulationSrn_2;
    typedef VertexBasedPopulationSrn<3> VertexBasedPopulationSrn_3;
    typedef AbstractUpdateRule<2> AbstractUpdateRule_2;
    typedef AbstractUpdateRule<3> AbstractUpdateRule_3;
    typedef AbstractCaUpdateRule<2> AbstractCaUpdateRule_2;
    typedef AbstractCaUpdateRule<3> AbstractCaUpdateRule_3;
    typedef AbstractPottsUpdateRule<2> AbstractPottsUpdateRule_2;
    typedef AbstractPottsUpdateRule<3> AbstractPottsUpdateRule_3;
    typedef AbstractCaSwitchingUpdateRule<2> AbstractCaSwitchingUpdateRule_2;
    typedef AbstractCaSwitchingUpdateRule<3> AbstractCaSwitchingUpdateRule_3;
    typedef AdhesionPottsUpdateRule<2> AdhesionPottsUpdateRule_2;
    typedef AdhesionPottsUpdateRule<3> AdhesionPottsUpdateRule_3;
    typedef ChemotaxisPottsUpdateRule<2> ChemotaxisPottsUpdateRule_2;
    typedef ChemotaxisPottsUpdateRule<3> ChemotaxisPottsUpdateRule_3;
    typedef DifferentialAdhesionPottsUpdateRule<2> DifferentialAdhesionPottsUpdateRule_2;
    typedef DifferentialAdhesionPottsUpdateRule<3> DifferentialAdhesionPottsUpdateRule_3;
    typedef DiffusionCaUpdateRule<2> DiffusionCaUpdateRule_2;
    typedef DiffusionCaUpdateRule<3> DiffusionCaUpdateRule_3;
    typedef RandomCaSwitchingUpdateRule<2> RandomCaSwitchingUpdateRule_2;
    typedef RandomCaSwitchingUpdateRule<3> RandomCaSwitchingUpdateRule_3;
    typedef SurfaceAreaConstraintPottsUpdateRule<2> SurfaceAreaConstraintPottsUpdateRule_2;
    typedef SurfaceAreaConstraintPottsUpdateRule<3> SurfaceAreaConstraintPottsUpdateRule_3;
    typedef VolumeConstraintPottsUpdateRule<2> VolumeConstraintPottsUpdateRule_2;
    typedef VolumeConstraintPottsUpdateRule<3> VolumeConstraintPottsUpdateRule_3;
    typedef AbstractCellBasedSimulation<2, 2> AbstractCellBasedSimulation_2_2;
    typedef AbstractCellBasedSimulation<3, 3> AbstractCellBasedSimulation_3_3;
    typedef OffLatticeSimulation<2, 2> OffLatticeSimulation_2_2;
    typedef OffLatticeSimulation<3, 3> OffLatticeSimulation_3_3;
    typedef OnLatticeSimulation<2> OnLatticeSimulation_2;
    typedef OnLatticeSimulation<3> OnLatticeSimulation_3;
    typedef AbstractCellBasedSimulationModifier<2, 2> AbstractCellBasedSimulationModifier_2_2;
    typedef AbstractCellBasedSimulationModifier<3, 3> AbstractCellBasedSimulationModifier_3_3;
    typedef AbstractTargetAreaModifier<2> AbstractTargetAreaModifier_2;
    typedef AbstractTargetAreaModifier<3> AbstractTargetAreaModifier_3;
    typedef DeltaNotchTrackingModifier<2> DeltaNotchTrackingModifier_2;
    typedef DeltaNotchTrackingModifier<3> DeltaNotchTrackingModifier_3;
    typedef DeltaNotchEdgeTrackingModifier<2> DeltaNotchEdgeTrackingModifier_2;
    typedef DeltaNotchEdgeTrackingModifier<3> DeltaNotchEdgeTrackingModifier_3;
    typedef DeltaNotchEdgeInteriorTrackingModifier<2> DeltaNotchEdgeInteriorTrackingModifier_2;
    typedef DeltaNotchEdgeInteriorTrackingModifier<3> DeltaNotchEdgeInteriorTrackingModifier_3;
    typedef DivisionBiasTrackingModifier<2> DivisionBiasTrackingModifier_2;
    typedef DivisionBiasTrackingModifier<3> DivisionBiasTrackingModifier_3;
    typedef ExtrinsicPullModifier<2> ExtrinsicPullModifier_2;
    typedef ExtrinsicPullModifier<3> ExtrinsicPullModifier_3;
    typedef ImmersedBoundarySimulationModifier<2> ImmersedBoundarySimulationModifier_2;
    typedef ImmersedBoundarySimulationModifier<3> ImmersedBoundarySimulationModifier_3;
    typedef ImmersedBoundarySvgWriter<2> ImmersedBoundarySvgWriter_2;
    typedef ImmersedBoundarySvgWriter<3> ImmersedBoundarySvgWriter_3;
    typedef NormallyDistributedTargetAreaModifier<2> NormallyDistributedTargetAreaModifier_2;
    typedef NormallyDistributedTargetAreaModifier<3> NormallyDistributedTargetAreaModifier_3;
    typedef SimpleTargetAreaModifier<2> SimpleTargetAreaModifier_2;
    typedef SimpleTargetAreaModifier<3> SimpleTargetAreaModifier_3;
    typedef TargetAreaLinearGrowthModifier<2> TargetAreaLinearGrowthModifier_2;
    typedef TargetAreaLinearGrowthModifier<3> TargetAreaLinearGrowthModifier_3;
    typedef VolumeTrackingModifier<2> VolumeTrackingModifier_2;
    typedef VolumeTrackingModifier<3> VolumeTrackingModifier_3;
    typedef AbstractNumericalMethod<2, 2> AbstractNumericalMethod_2_2;
    typedef AbstractNumericalMethod<3, 3> AbstractNumericalMethod_3_3;
    typedef ForwardEulerNumericalMethod<2, 2> ForwardEulerNumericalMethod_2_2;
    typedef ForwardEulerNumericalMethod<3, 3> ForwardEulerNumericalMethod_3_3;
    typedef AbstractCellBasedWriter<2, 2> AbstractCellBasedWriter_2_2;
    typedef AbstractCellBasedWriter<3, 3> AbstractCellBasedWriter_3_3;
    typedef AbstractCellWriter<2, 2> AbstractCellWriter_2_2;
    typedef AbstractCellWriter<3, 3> AbstractCellWriter_3_3;
    typedef CellAgesWriter<2, 2> CellAgesWriter_2_2;
    typedef CellAgesWriter<3, 3> CellAgesWriter_3_3;
    typedef CellAncestorWriter<2, 2> CellAncestorWriter_2_2;
    typedef CellAncestorWriter<3, 3> CellAncestorWriter_3_3;
    typedef CellAppliedForceWriter<2, 2> CellAppliedForceWriter_2_2;
    typedef CellAppliedForceWriter<3, 3> CellAppliedForceWriter_3_3;
    typedef CellCycleModelProteinConcentrationsWriter<2, 2> CellCycleModelProteinConcentrationsWriter_2_2;
    typedef CellCycleModelProteinConcentrationsWriter<3, 3> CellCycleModelProteinConcentrationsWriter_3_3;
    typedef CellDataItemWriter<2, 2> CellDataItemWriter_2_2;
    typedef CellDataItemWriter<3, 3> CellDataItemWriter_3_3;
    typedef CellDeltaNotchWriter<2, 2> CellDeltaNotchWriter_2_2;
    typedef CellDeltaNotchWriter<3, 3> CellDeltaNotchWriter_3_3;
    typedef CellIdWriter<2, 2> CellIdWriter_2_2;
    typedef CellIdWriter<3, 3> CellIdWriter_3_3;
    typedef CellLabelWriter<2, 2> CellLabelWriter_2_2;
    typedef CellLabelWriter<3, 3> CellLabelWriter_3_3;
    typedef CellLocationIndexWriter<2, 2> CellLocationIndexWriter_2_2;
    typedef CellLocationIndexWriter<3, 3> CellLocationIndexWriter_3_3;
    typedef CellMutationStatesWriter<2, 2> CellMutationStatesWriter_2_2;
    typedef CellMutationStatesWriter<3, 3> CellMutationStatesWriter_3_3;
    typedef CellProliferativePhasesWriter<2, 2> CellProliferativePhasesWriter_2_2;
    typedef CellProliferativePhasesWriter<3, 3> CellProliferativePhasesWriter_3_3;
    typedef CellProliferativeTypesWriter<2, 2> CellProliferativeTypesWriter_2_2;
    typedef CellProliferativeTypesWriter<3, 3> CellProliferativeTypesWriter_3_3;
    typedef CellRadiusWriter<2, 2> CellRadiusWriter_2_2;
    typedef CellRadiusWriter<3, 3> CellRadiusWriter_3_3;
    typedef CellRosetteRankWriter<2, 2> CellRosetteRankWriter_2_2;
    typedef CellRosetteRankWriter<3, 3> CellRosetteRankWriter_3_3;
    typedef CellVolumesWriter<2, 2> CellVolumesWriter_2_2;
    typedef CellVolumesWriter<3, 3> CellVolumesWriter_3_3;
    typedef ImmersedBoundaryBoundaryCellWriter<2, 2> ImmersedBoundaryBoundaryCellWriter_2_2;
    typedef ImmersedBoundaryBoundaryCellWriter<3, 3> ImmersedBoundaryBoundaryCellWriter_3_3;
    typedef ImmersedBoundaryNeighbourNumberWriter<2, 2> ImmersedBoundaryNeighbourNumberWriter_2_2;
    typedef ImmersedBoundaryNeighbourNumberWriter<3, 3> ImmersedBoundaryNeighbourNumberWriter_3_3;
    typedef LegacyCellProliferativeTypesWriter<2, 2> LegacyCellProliferativeTypesWriter_2_2;
    typedef LegacyCellProliferativeTypesWriter<3, 3> LegacyCellProliferativeTypesWriter_3_3;
    typedef AbstractCellPopulationWriter<2, 2> AbstractCellPopulationWriter_2_2;
    typedef AbstractCellPopulationWriter<3, 3> AbstractCellPopulationWriter_3_3;
    typedef BoundaryNodeWriter<2, 2> BoundaryNodeWriter_2_2;
    typedef BoundaryNodeWriter<3, 3> BoundaryNodeWriter_3_3;
    typedef CellPopulationAdjacencyMatrixWriter<2, 2> CellPopulationAdjacencyMatrixWriter_2_2;
    typedef CellPopulationAdjacencyMatrixWriter<3, 3> CellPopulationAdjacencyMatrixWriter_3_3;
    typedef CellPopulationAreaWriter<2, 2> CellPopulationAreaWriter_2_2;
    typedef CellPopulationAreaWriter<3, 3> CellPopulationAreaWriter_3_3;
    typedef CellPopulationElementWriter<2, 2> CellPopulationElementWriter_2_2;
    typedef CellPopulationElementWriter<3, 3> CellPopulationElementWriter_3_3;
    typedef HeterotypicBoundaryLengthWriter<2, 2> HeterotypicBoundaryLengthWriter_2_2;
    typedef HeterotypicBoundaryLengthWriter<3, 3> HeterotypicBoundaryLengthWriter_3_3;
    typedef NodeLocationWriter<2, 2> NodeLocationWriter_2_2;
    typedef NodeLocationWriter<3, 3> NodeLocationWriter_3_3;
    typedef NodeVelocityWriter<2, 2> NodeVelocityWriter_2_2;
    typedef NodeVelocityWriter<3, 3> NodeVelocityWriter_3_3;
    typedef RadialCellDataDistributionWriter<2, 2> RadialCellDataDistributionWriter_2_2;
    typedef RadialCellDataDistributionWriter<3, 3> RadialCellDataDistributionWriter_3_3;
    typedef VertexIntersectionSwapLocationsWriter<2, 2> VertexIntersectionSwapLocationsWriter_2_2;
    typedef VertexIntersectionSwapLocationsWriter<3, 3> VertexIntersectionSwapLocationsWriter_3_3;
    typedef VertexT1SwapLocationsWriter<2, 2> VertexT1SwapLocationsWriter_2_2;
    typedef VertexT1SwapLocationsWriter<3, 3> VertexT1SwapLocationsWriter_3_3;
    typedef VertexT2SwapLocationsWriter<2, 2> VertexT2SwapLocationsWriter_2_2;
    typedef VertexT2SwapLocationsWriter<3, 3> VertexT2SwapLocationsWriter_3_3;
    typedef VertexT3SwapLocationsWriter<2, 2> VertexT3SwapLocationsWriter_2_2;
    typedef VertexT3SwapLocationsWriter<3, 3> VertexT3SwapLocationsWriter_3_3;
    typedef VoronoiDataWriter<2, 2> VoronoiDataWriter_2_2;
    typedef VoronoiDataWriter<3, 3> VoronoiDataWriter_3_3;
    typedef AbstractCellPopulationCountWriter<2, 2> AbstractCellPopulationCountWriter_2_2;
    typedef AbstractCellPopulationCountWriter<3, 3> AbstractCellPopulationCountWriter_3_3;
    typedef CellMutationStatesCountWriter<2, 2> CellMutationStatesCountWriter_2_2;
    typedef CellMutationStatesCountWriter<3, 3> CellMutationStatesCountWriter_3_3;
    typedef CellProliferativePhasesCountWriter<2, 2> CellProliferativePhasesCountWriter_2_2;
    typedef CellProliferativePhasesCountWriter<3, 3> CellProliferativePhasesCountWriter_3_3;
    typedef CellProliferativeTypesCountWriter<2, 2> CellProliferativeTypesCountWriter_2_2;
    typedef CellProliferativeTypesCountWriter<3, 3> CellProliferativeTypesCountWriter_3_3;
    typedef AbstractCellPopulationEventWriter<2, 2> AbstractCellPopulationEventWriter_2_2;
    typedef AbstractCellPopulationEventWriter<3, 3> AbstractCellPopulationEventWriter_3_3;
    typedef CellDivisionLocationsWriter<2, 2> CellDivisionLocationsWriter_2_2;
    typedef CellDivisionLocationsWriter<3, 3> CellDivisionLocationsWriter_3_3;
    typedef CellRemovalLocationsWriter<2, 2> CellRemovalLocationsWriter_2_2;
    typedef CellRemovalLocationsWriter<3, 3> CellRemovalLocationsWriter_3_3;
    typedef AttractingPlaneBoundaryCondition<2, 2> AttractingPlaneBoundaryCondition_2_2;
    typedef AttractingPlaneBoundaryCondition<3, 3> AttractingPlaneBoundaryCondition_3_3;
    typedef VtkSceneModifier<2> VtkSceneModifier_2;
    typedef VtkSceneModifier<3> VtkSceneModifier_3;
    typedef PythonSimulationModifier<2> PythonSimulationModifier_2;
    typedef PythonSimulationModifier<3> PythonSimulationModifier_3;
    typedef VtkScene<2> VtkScene_2;
    typedef VtkScene<3> VtkScene_3;
    typedef AbstractPyChasteActorGenerator<2> AbstractPyChasteActorGenerator_2;
    typedef AbstractPyChasteActorGenerator<3> AbstractPyChasteActorGenerator_3;
    typedef CellPopulationPyChasteActorGenerator<2> CellPopulationPyChasteActorGenerator_2;
    typedef CellPopulationPyChasteActorGenerator<3> CellPopulationPyChasteActorGenerator_3;
} // namespace cppwg

#endif // pychaste_HEADERS_HPP_
