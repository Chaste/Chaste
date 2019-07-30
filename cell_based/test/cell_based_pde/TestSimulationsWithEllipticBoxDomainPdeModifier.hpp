/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef TESTSIMULATIONSWITHELLIPTICBOXDOMAINPDEMODIFIER_HPP_
#define TESTSIMULATIONSWITHELLIPTICBOXDOMAINPDEMODIFIER_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellIdWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "VolumeTrackingModifier.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"
#include "CellLabel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "CaBasedCellPopulation.hpp"
#include "DiffusionCaUpdateRule.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

static const double M_TIME_FOR_SIMULATION = 1.0;

static const double M_NUM_CELLS_ACROSS = 5; // this is 5^2 initial cells

class TestSimulationWithEllipticBoxDomainPdeModifier : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, double EquilibriumVolume)
    {
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        for (unsigned i=0; i<num_cells; i++)
        {
            SimpleOxygenBasedCellCycleModel* p_cycle_model = new SimpleOxygenBasedCellCycleModel();
            p_cycle_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            rCells.push_back(p_cell);
        }
     }

public:

    void TestEllipticBoxDomainPdeModifierWithVertexBasedMonolayer()
    {
        // Create Mesh
        HoneycombVertexMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create Cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells, 1.0);

        // Create Population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EllipticBoxMonolayers/Vertex");
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Create Forces and pass to simulation NOTE: PARAMETERS CHOSEN TO GET CIRCULAR MONOLAYER
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(6.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(12.0);
        simulator.AddForce(p_force);

        // Create Modifiers and pass to simulation

        // Create a PDE modifier and pass it to the simulation Add this first so in place for SimpleTargetArea one (calls cell pop update)

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();

        // Test some simulation statistics
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 25u); //No birth yet

        // Test new locations and concentrations (i.e that the results haven't changed since last revision.)
        std::vector<double> node_5_location = simulator.GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.7247, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 0.7715, 1e-4);

        // Note this is cell associated with element 5 not node 5
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(5))->GetCellData()->GetItem("oxygen"), 0.6191, 1e-4);
    }

    void TestEllipticBoxDomainPdeModifierWithNodeBasedMonolayer()
    {
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS,0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells,1.0);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EllipticBoxMonolayers/Node");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(RepulsionForce<2>, p_force);
        simulator.AddForce(p_force);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();

        // Test some simulation statistics
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 25u); // No birth yet

        // Test new locations and concentrations (i.e that the results haven't changed since last revision.)
        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 4.0, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.0, 1e-4);
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(4))->GetCellData()->GetItem("oxygen"), 0.6786, 1e-4);

        delete p_mesh; // to stop memory leaks
    }

    void TestEllipticBoxDomainPdeModifierWithMeshBasedMonolayer()
    {
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS,M_NUM_CELLS_ACROSS,0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells,1.0);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EllipticBoxMonolayers/MeshPoint");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();

        // Test some simulation statistics
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 25u); // No birth yet

        // Test new locations and concentrations (i.e that the results haven't changed since last revision.)
        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 4.0, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.0, 1e-4);
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(4))->GetCellData()->GetItem("oxygen"), 0.6786, 1e-4);
    }

    void TestEllipticBoxDomainPdeModifierWithMeshBasedWithGhostNodesBasedMonolayer()
    {
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS,M_NUM_CELLS_ACROSS,2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cell population
        std::vector<CellPtr> cells;
        GenerateCells(location_indices.size(),cells,1.0);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EllipticBoxMonolayers/MeshWithGhosts");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();

        // Test some simulation statistics
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 25u); // No birth yet

        // Test new locations and concentrations (i.e that the results haven't changed since last revision.)
        std::vector<double> real_node_4_location = simulator.GetNodeLocation(24);
        TS_ASSERT_DELTA(real_node_4_location[0], 4.0, 1e-4);
        TS_ASSERT_DELTA(real_node_4_location[1], 0.0, 1e-4);
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(24))->GetCellData()->GetItem("oxygen"), 0.6786, 1e-4);
    }

    void TestEllipticBoxDomainPdeModifierWithPottsBasedMonolayer()
    {
        unsigned cell_width = 4;
        unsigned domain_width = 200;
        PottsMeshGenerator<2> generator(domain_width, M_NUM_CELLS_ACROSS, cell_width, domain_width, M_NUM_CELLS_ACROSS, cell_width);
        PottsMesh<2>* p_mesh = generator.GetMesh();
        p_mesh->Translate(-0.5*(domain_width-(M_NUM_CELLS_ACROSS-1)*cell_width),-0.5*(domain_width-(M_NUM_CELLS_ACROSS-1)*cell_width));
        p_mesh->Scale(0.25,0.25);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells,16);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(0.1);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EllipticBoxMonolayers/Potts");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule);
        simulator.AddUpdateRule(p_surface_area_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();

        // Test some simulation statistics
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 25u); // No birth yet

        // Test new locations and concentrations (i.e that the results haven't changed since last revision.)
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(1))->GetCellData()->GetItem("oxygen"), 0.6606, 1e-4);
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(4))->GetCellData()->GetItem("oxygen"), 0.7123, 1e-4);
    }

    void TestEllipticBoxDomainPdeModifierWithCaBasedMonolayer()
    {
        // Create a simple 2D PottsMesh
        unsigned domain_width = 5*M_NUM_CELLS_ACROSS;

        PottsMeshGenerator<2> generator(domain_width, 0, 0, domain_width, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();
        p_mesh->Translate(-0.5*(domain_width-M_NUM_CELLS_ACROSS),-0.5*(domain_width-M_NUM_CELLS_ACROSS));

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<M_NUM_CELLS_ACROSS; i++)
        {
            for (unsigned j=0; j<M_NUM_CELLS_ACROSS; j++)
            {
                unsigned offset = (domain_width+1) * (domain_width-M_NUM_CELLS_ACROSS)/2;
                location_indices.push_back(offset + j + i * domain_width );
            }
        }

        std::vector<CellPtr> cells;
        GenerateCells(location_indices.size(),cells,1);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EllipticBoxMonolayers/CA");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Adding update rule
        MAKE_PTR(DiffusionCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.1);
        simulator.AddUpdateRule(p_diffusion_update_rule);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Run simulation
        simulator.Solve();

        // Test some simulation statistics
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 25u);

        // Test new locations and concentrations (i.e that the results haven't changed since last revision.)
        AbstractCellPopulation<2, 2>::Iterator cell_iter=simulator.rGetCellPopulation().Begin();
        TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("oxygen"), 0.6754, 1e-4);
        ++cell_iter;
        ++cell_iter;
        ++cell_iter;
        ++cell_iter;
        TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("oxygen"), 0.7106, 1e-4);
    }
};

#endif /*TESTSIMULATIONSWITHELLIPTICBOXDOMAINPDEMODIFIER_HPP_*/
