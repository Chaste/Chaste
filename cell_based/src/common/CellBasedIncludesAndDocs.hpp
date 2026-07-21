/*

Copyright (c) 2005-2026, University of Oxford.
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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef CELLBASEDINCLUDESANDDOCS_HPP_
#define CELLBASEDINCLUDESANDDOCS_HPP_

/* We include the TestSuite and a wide variety of useful includes
 * for cell-based simulations, with documentation. This file may
 * be included into any new project looking to explore the cell-based
 * functionality offered by Chaste.
 */

#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"

#include "CellBasedSimulationArchiver.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"

#include "AdhesionPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"

#include "CellsGenerator.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"

// Standard family of CellCycleModels, with explanation:

/* Cell cycle models agnostic to G1/S/G2/M i.e. not phase-based
 * Only relevant quantity is cycle duration, when age exceeds this value
 * ReadyToDivide() returns true by default
 * Default Bernoulli division 'probability' is p=0.1
 */
#include "AlwaysDivideCellCycleModel.hpp" // always returns true for ReadyToDivide
#include "NoCellCycleModel.hpp" // always returns false for ReadyToDivide
#include "UniformCellCycleModel.hpp" // duration ~ U[12, 14] hours then is ReadyToDivide
#include "BernoulliTrialCellCycleModel.hpp" // waits one hour, then P(ReadyToDivide)=p*GetTimeStep()
#include "BiasedBernoulliTrialCellCycleModel.hpp" // ditto, with P=CellData()["bias"] * p*GetTimeStep()
#include "LabelDependentBernoulliTrialCellCycleModel.hpp" // ditto, with another p used if CellLabel exists

/* Models which expose individual phases G1/S/G2/M
 * DifferentiatedCellProliferativeType has a G0 phase as well, which is eternal
 * (without user interference) and in which mitosis never occurs
 *
 * Current (perhaps dubious) defaults (in hours):
 * (idea is to have G1 stochastically vary or be set, all others are hard-coded)
 *  Cell cycle phase = M_PHASE
 *  Minimum gap = 0.01
 *  Stem cell generic G1 duration = 14.0
 *  Transit cell generic G1 duration = 2.0
 *  G1 duration = Unset (derived classes compute it using the above and their own logic)
 *
 *  S duration = 5.0
 *  G2 duration = 4.0
 *  M duration = 1.0
 *
 *  For generational models, after the default n=3 generations (incremented in mitosis)
 *  the daughter cells differentiate and mitosis ceases
 */
#include "UniformG1GenerationalCellCycleModel.hpp" // G1 ~ U[2,4] or U[14,18] for Transit, Stem
#include "GammaG1CellCycleModel.hpp" // G1 ~ Gamma(user-specified shape, scale)
#include "ExponentialG1GenerationalCellCycleModel.hpp" // G1 ~ Exp(1/transit-cell g1 duration)
#include "FixedSequenceCellCycleModel.hpp" // G1 ~ fixed, singleton vector of random Exponential-G1 times

// More exotic phase-based models
#include "ContactInhibitionCellCycleModel.hpp" /* Has a 'quiescence' (dormancy) period
 * i.e. a G0 phase is included in the cell cycle. Quiescence induced by physical stress:
 * compressed/stressed cells remain dormant for longer.
 *
 * User-specified equilibrium volume in G1 phase, and non-dimensionalised quiescent volume fraction.
 *
 * These cells are labelled when quiescent, unlabelled when not
 * and when the CellData "volume" property is below the quiescence threshold,
 * the G1 duration parameter is incremented from the default baseline to prevent
 * passage into the next phase; cells stressed in their infancy then take longer to mature
 *
 * If the cell is differentiated, it is permanently in G0-phase.
 */
#include "SimpleOxygenBasedCellCycleModel.hpp" /* Has a quiescence feature as well
 * affected by oxygen levels, as well as a hypoxia feature which triggers cell death
 * defaults: 2 hours of hypoxia will trigger apoptosis
 * non-dimensionalised concentration thresholds:
 * 1.0 (for quiescence), 0.4 (for hypoxia)
 *
 * Oxygen concentration gotten through CellData()->GetItem("oxygen"), so the user's
 * simulation will need to assign this property.
 *
 * If the cell is not apoptotic, then when in G1 phase its G1 duration is incremented
 * when [O] < Q, Q being the threshold concentration for quiescence, with increment
 * (timestep) * (1 - [O]/Q), which is usually less than a timestep.
 * G1 does have a headstart of, say, 2 hours, but we see the cell will naturally pass into S-phase
 * even in low oxygen conditions, just more slowly. This is in contrast with the ContactInhibition
 * model, which increments its G1 duration by a full timestep each timestep, when stressed.
 *
 * When [O] < H, the hypoxia threshold (default 0.4), the cell is labelled with the apoptotic property
 * i.e. sentenced to death with probability p = 0.9 - (1/2) * [O]/H if the duration of hypoxia exceeds
 * the critical value (default, 2 hours)
 */

// ODE Cell cycle models

/* Update various chemical values and parameters
 * each timestep, and have their own logic for defining a "stopping event"
 * i.e. a significant chemical event which triggers the passage from G1 to S-phase
 *
 * The abstract ode cell cycle models use the "stopping event" as the mechanic when
 * overwriting the UpdateCellCyclePhase() method. The WNT signalling pathway ode models
 * are more complex and have additional considerations, such as a WNT concentration object,
 * and most of them are stored in 'crypt' rather than 'cell-based'.
 */
#include "Alarcon2004OxygenBasedCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"

#include "RandomCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "AbstractForce.hpp"
#include "NagaiHondaForce.hpp"
#include "RepulsionForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "VoronoiDataWriter.hpp"

#include "SmartPointers.hpp" // MAKE_PTR(type, variable name) must be used whenever creating a new class instance


#endif // CELLBASEDINCLUDESANDDOCS_HPP_
