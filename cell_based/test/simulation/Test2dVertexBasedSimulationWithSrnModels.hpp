/*

Copyright (c) 2005-2021, University of Oxford.
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

#ifndef TESTDELTANOTCHEDGEODESIMULATION_HPP_
#define TESTDELTANOTCHEDGEODESIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "WelikyOsterForce.hpp"
#include "FarhadifarForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellVolumesWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"

#include "CellSrnModel.hpp"

#include "DeltaNotchEdgeSrnModel.hpp"
#include "DeltaNotchEdgeTrackingModifier.hpp"

#include "RDEdgeSrnModel.hpp"
#include "RDEdgeTrackingModifier.hpp"


/**
 * These tests check and demonstrate simulation of vertex based models with edge  Srn models
 */
class TestDeltaNotchEdgeOnlyODESimulation : public AbstractCellBasedTestSuite
{
public:
    /*
     * Test vertex based simulations when both edge AND interior SRN models are specified
     */
    // void TestRunningMultiODECellWithEdgesAndInterior()
    // {
    //     EXIT_IF_PARALLEL;
    //     /* First we create a regular vertex mesh. */
    //     HoneycombVertexMeshGenerator generator(6, 6);
    //     MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

    //     std::vector<CellPtr> cells;
    //     MAKE_PTR(WildTypeCellMutationState, p_state);
    //     MAKE_PTR(TransitCellProliferativeType, p_diff_type);

    //     for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
    //     {
    //         /* Initalise cell cycle */
    //         UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
    //         p_cc_model->SetDimension(2);

    //         auto p_element = p_mesh->GetElement(elem_index);
    //         /* Initialise edge based SRN */
    //         auto p_cell_edge_srn_model = new CellSrnModel();
    //         /* We choose to initialise the total concentrations to random levels */
    //         auto delta_concentration = RandomNumberGenerator::Instance()->ranf();
    //         auto notch_concentration = RandomNumberGenerator::Instance()->ranf();

    //         double total_edge_length = 0.0;
    //         for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
    //         {
    //             total_edge_length += p_element->GetEdge(i)->rGetLength();
    //         }

    //         /* Gets the edges of the element and create an SRN for each edge */
    //         for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
    //         {
    //             auto p_elem_edge = p_element->GetEdge(i);
    //             auto p_edge_length = p_elem_edge->rGetLength();
    //             std::vector<double> initial_conditions;

    //             /* Initial concentration of delta and notch vary depending on the edge length */
    //             initial_conditions.push_back( p_edge_length/total_edge_length * delta_concentration);
    //             initial_conditions.push_back( p_edge_length/total_edge_length * notch_concentration);

    //             MAKE_PTR(DeltaNotchEdgeSrnModel, p_srn_model);
    //             p_srn_model->SetInitialConditions(initial_conditions);
    //             p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
    //         }
    //         //Add interior SRN models to cells
    //         MAKE_PTR(DeltaNotchInteriorSrnModel, p_cell_srn_model);
    //         std::vector<double> zero_conditions(2);
    //         p_cell_srn_model->SetInitialConditions(zero_conditions);
    //         p_cell_edge_srn_model->SetInteriorSrnModel(p_cell_srn_model);

    //         CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
    //         p_cell->SetCellProliferativeType(p_diff_type);

    //         double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
    //         p_cell->SetBirthTime(birth_time);
    //         cells.push_back(p_cell);
    //     }
    //     /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
    //      * output to file. */
    //     VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
    //     cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
    //     cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
    //     cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
    //     cell_population.AddCellWriter<CellProliferativePhasesWriter>();
    //     cell_population.AddCellWriter<CellAgesWriter>();
    //     cell_population.AddCellWriter<CellVolumesWriter>();


    //     /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
    //      * and run the simulation.*/
    //     OffLatticeSimulation<2> simulator(cell_population);
    //     simulator.SetOutputDirectory("TestDeltaNotchEdgeInteriorODESimulation");
    //     simulator.SetSamplingTimestepMultiple(10);
    //     simulator.SetEndTime(10.0);

    //     /* Update CellData and CellEdgeData so that SRN simulations can run properly */
    //     MAKE_PTR(DeltaNotchEdgeInteriorTrackingModifier<2>, p_cell_modifier);
    //     simulator.AddSimulationModifier(p_cell_modifier);
    //     MAKE_PTR(DeltaNotchEdgeTrackingModifier<2>, p_edge_modifier);
    //     simulator.AddSimulationModifier(p_edge_modifier);

    //     MAKE_PTR(NagaiHondaForce<2>, p_force);
    //     simulator.AddForce(p_force);

    //     /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
    //      */
    //     MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
    //     simulator.AddSimulationModifier(p_growth_modifier);
    //     TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    // }

// /*
//      * Test vertex based simulations when both edge AND interior SRN models are specified
//      */
//     void TestRunningMultiODECellWithEdgesAndInteriorCellSplit()
//     {
//         EXIT_IF_PARALLEL;
//         /* First we create a regular vertex mesh. */
//         HoneycombVertexMeshGenerator generator(6, 6);
//         MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

//         p_mesh->SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);

//         std::vector<CellPtr> cells;
//         MAKE_PTR(WildTypeCellMutationState, p_state);
//         MAKE_PTR(TransitCellProliferativeType, p_diff_type);

//         for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
//         {
//             /* Initalise cell cycle */
//             UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
//             p_cc_model->SetDimension(2);

//             auto p_element = p_mesh->GetElement(elem_index);
//             /* Initialise edge based SRN */
//             auto p_cell_edge_srn_model = new CellSrnModel();
//             /* We choose to initialise the total concentrations to random levels */
//             //auto delta_concentration = RandomNumberGenerator::Instance()->ranf();
//             //auto notch_concentration = RandomNumberGenerator::Instance()->ranf();

//             auto delta_concentration = 1.0;
//             auto notch_concentration = 1.0;

//             double total_edge_length = 0.0;
//             for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
//             {
//                 total_edge_length += p_element->GetEdge(i)->rGetLength();
//             }

//             /* Gets the edges of the element and create an SRN for each edge */
//             for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
//             {
//                 auto p_elem_edge = p_element->GetEdge(i);
//                 auto p_edge_length = p_elem_edge->rGetLength();
//                 std::vector<double> initial_conditions;

//                 /* Initial concentration of delta and notch vary depending on the edge length */
//                 initial_conditions.push_back( p_edge_length/total_edge_length * delta_concentration);
//                 initial_conditions.push_back( p_edge_length/total_edge_length * notch_concentration);

//                 MAKE_PTR(DeltaNotchEdgeSrnModel, p_srn_model);
//                 p_srn_model->SetInitialConditions(initial_conditions);
//                 p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
//             }
//             //Add interior SRN models to cells
//             MAKE_PTR(DeltaNotchInteriorSrnModel, p_cell_srn_model);
//             std::vector<double> zero_conditions(2);
//             p_cell_srn_model->SetInitialConditions(zero_conditions);
//             p_cell_edge_srn_model->SetInteriorSrnModel(p_cell_srn_model);

//             CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
//             p_cell->SetCellProliferativeType(p_diff_type);

//             //double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
//             double birth_time = -1.0;
//             p_cell->SetBirthTime(birth_time);
//             cells.push_back(p_cell);
//         }
//         /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
//          * output to file. */
//         VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
//         cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
//         cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
//         cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
//         cell_population.AddCellWriter<CellProliferativePhasesWriter>();
//         cell_population.AddCellWriter<CellAgesWriter>();
//         cell_population.AddCellWriter<CellVolumesWriter>();


//         /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
//          * and run the simulation.*/
//         OffLatticeSimulation<2> simulator(cell_population);
//         simulator.SetOutputDirectory("TestDeltaNotchEdgeInteriorODESimulationEdgeSplit");
//         simulator.SetSamplingTimestepMultiple(5);
//         simulator.SetEndTime(10.0);

//         /* Update CellData and CellEdgeData so that SRN simulations can run properly */
//         MAKE_PTR(DeltaNotchEdgeInteriorTrackingModifier<2>, p_cell_modifier);
//         simulator.AddSimulationModifier(p_cell_modifier);
//         MAKE_PTR(DeltaNotchEdgeTrackingModifier<2>, p_edge_modifier);
//         simulator.AddSimulationModifier(p_edge_modifier);

//         MAKE_PTR(WelikyOsterForce<2>, p_force);
//         simulator.AddForce(p_force);

//         /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
//          */
//         MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
//         simulator.AddSimulationModifier(p_growth_modifier);
//         TS_ASSERT_THROWS_NOTHING(simulator.Solve());
//     }

//     /*
//      * Test vertex based simulations when both edge AND interior SRN models are specified and checking for T1 swaps ( Found at 1725 - 1726 on LHS of sim)
//      */
//     void TestRunningMultiODECellWithEdgesAndInteriorCellT1()
//     {
//         EXIT_IF_PARALLEL;
//         /* First we create a regular vertex mesh. */
//         HoneycombVertexMeshGenerator generator(2, 2);
//         MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

//         p_mesh->SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);

//         std::vector<CellPtr> cells;
//         MAKE_PTR(WildTypeCellMutationState, p_state);
//         MAKE_PTR(TransitCellProliferativeType, p_diff_type);

//         for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
//         {
//             /* Initalise cell cycle */
//             UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
//             p_cc_model->SetDimension(2);

//             auto p_element = p_mesh->GetElement(elem_index);
//             /* Initialise edge based SRN */
//             auto p_cell_edge_srn_model = new CellSrnModel();
//             /* We choose to initialise the total concentrations to random levels */
//             //auto delta_concentration = RandomNumberGenerator::Instance()->ranf();
//             //auto notch_concentration = RandomNumberGenerator::Instance()->ranf();

//             auto delta_concentration = 1.0;
//             auto notch_concentration = 1.0;

//             double total_edge_length = 0.0;
//             for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
//             {
//                 total_edge_length += p_element->GetEdge(i)->rGetLength();
//             }

//             /* Gets the edges of the element and create an SRN for each edge */
//             for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
//             {
//                 auto p_elem_edge = p_element->GetEdge(i);
//                 auto p_edge_length = p_elem_edge->rGetLength();
//                 std::vector<double> initial_conditions;

//                 /* Initial concentration of delta and notch vary depending on the edge length */
//                 initial_conditions.push_back( p_edge_length/total_edge_length * delta_concentration);
//                 initial_conditions.push_back( p_edge_length/total_edge_length * notch_concentration);

//                 MAKE_PTR(DeltaNotchEdgeSrnModel, p_srn_model);
//                 p_srn_model->SetInitialConditions(initial_conditions);
//                 p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
//             }
//             //Add interior SRN models to cells
//             MAKE_PTR(DeltaNotchInteriorSrnModel, p_cell_srn_model);
//             std::vector<double> zero_conditions(2);
//             p_cell_srn_model->SetInitialConditions(zero_conditions);
//             p_cell_edge_srn_model->SetInteriorSrnModel(p_cell_srn_model);

//             CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
//             p_cell->SetCellProliferativeType(p_diff_type);

//             //double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
//             double birth_time = -5.0;
//             p_cell->SetBirthTime(birth_time);
//             cells.push_back(p_cell);
//         }
//         /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
//          * output to file. */
//         VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
//         cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
//         cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
//         cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
//         cell_population.AddCellWriter<CellProliferativePhasesWriter>();
//         cell_population.AddCellWriter<CellAgesWriter>();
//         cell_population.AddCellWriter<CellVolumesWriter>();


//         /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
//          * and run the simulation.*/
//         OffLatticeSimulation<2> simulator(cell_population);
//         simulator.SetOutputDirectory("TestDeltaNotchEdgeInteriorODESimulationEdgeT1");
//         simulator.SetSamplingTimestepMultiple(10);
//         simulator.SetEndTime(50.0);

//         /* Update CellData and CellEdgeData so that SRN simulations can run properly */
//         MAKE_PTR(DeltaNotchEdgeInteriorTrackingModifier<2>, p_cell_modifier);
//         simulator.AddSimulationModifier(p_cell_modifier);
//         MAKE_PTR(DeltaNotchEdgeTrackingModifier<2>, p_edge_modifier);
//         simulator.AddSimulationModifier(p_edge_modifier);

//         MAKE_PTR(FarhadifarForce<2>, p_force);
//         simulator.AddForce(p_force);

//         /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
//          */
//         MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
//         simulator.AddSimulationModifier(p_growth_modifier);
//         TS_ASSERT_THROWS_NOTHING(simulator.Solve());
//     }

    /*
     * Test whether running vertex based model with edge based SRN.
     */
    void TestRunningMultiODECellWithEdges()
    {
        /* First we create a regular vertex mesh. */
        HoneycombVertexMeshGenerator generator(4, 4);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
            /* Initalise cell cycle */
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);

            auto p_element = p_mesh->GetElement(elem_index);

            /* Initialise edge based SRN */
            auto p_cell_edge_srn_model = new CellSrnModel();

            /* We choose to initialise the total concentrations to random levels */
            auto delta_concentration = RandomNumberGenerator::Instance()->ranf();
            auto notch_concentration = RandomNumberGenerator::Instance()->ranf();

            double total_edge_length = 0.0;
            for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
            {
                total_edge_length += p_element->GetEdge(i)->rGetLength();
            }

            /* Gets the edges of the element and create an SRN for each edge */
            for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
            {
                auto p_elem_edge = p_element->GetEdge(i);
                auto p_edge_length = p_elem_edge->rGetLength();
                std::vector<double> initial_conditions;

                /* Initial concentration of delta and notch vary depending on the edge length */
                initial_conditions.push_back( p_edge_length/total_edge_length * delta_concentration);
                initial_conditions.push_back( p_edge_length/total_edge_length * notch_concentration);

                MAKE_PTR(DeltaNotchEdgeSrnModel, p_srn_model);
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
            }

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchEdgeOnlyODESimulation");
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(50.0);

        /* Update CellEdgeData so that SRN simulations can run properly */
        MAKE_PTR(DeltaNotchEdgeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        simulator.Solve();
        //TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }

};


#endif /*TESTDELTANOTCHEDGEODESIMULATION_HPP_*/
