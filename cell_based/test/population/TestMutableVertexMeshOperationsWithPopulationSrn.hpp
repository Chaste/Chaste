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

#ifndef TESTMUTABLEVERTEXMESHOPERATIONSWITHPOPULATIONSRN_HPP_
#define TESTMUTABLEVERTEXMESHOPERATIONSWITHPOPULATIONSRN_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "VertexMeshWriter.hpp"
#include "MutableVertexMesh.hpp"
#include "FileComparison.hpp"
#include "Warnings.hpp"

#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"
#include "AlwaysDivideCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"

#include "VertexBasedCellPopulation.hpp"
#include "DeltaNotchSrnEdgeModel.hpp"

#include "SrnCellModel.hpp"
#include "DeltaNotchSrnInteriorModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "DeltaNotchCellEdgeTrackingModifier.hpp"
#include "DeltaNotchEdgeInteriorTrackingModifier.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

// This test is designed to check if SRN update after mesh operations works correctly
class TestMutableVertexMeshOperationsWithPopulationSrn : public AbstractCellBasedTestSuite
{
public:

    void TestPerformNodeMergeWithSrn()
    {
        /*
         * Create a mesh comprising a single triangular element, as shown below.
         * We will test that the nodes marked with an x are merged correctly.
         *
         *      /|
         *     / |
         *    /  |
         *   /   |
         *  /    |
         *  --xx-
         */
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.4, 0.0));
        nodes.push_back(new Node<2>(4, true, 0.6, 0.0));

        unsigned node_indices_elem_0[5] = {0, 3, 4, 1, 2};
        std::vector<Node<2>*> nodes_elem_0;
        for (unsigned i=0; i<5; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Create cell edge SRN with four edges
        auto p_cell_srn_model = new SrnCellModel();
        for (int i = 0; i < 5; i++)
        {
            boost::shared_ptr<DeltaNotchSrnEdgeModel> p_delta_notch_edge_srn_model(new DeltaNotchSrnEdgeModel());

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(1.0);
            starter_conditions.push_back(1.0);
            p_delta_notch_edge_srn_model->SetInitialConditions(starter_conditions);
            p_cell_srn_model->AddEdgeSrnModel(p_delta_notch_edge_srn_model);
        }
        boost::shared_ptr<DeltaNotchSrnInteriorModel> p_delta_notch_interior_srn_model(new DeltaNotchSrnInteriorModel());

        // Create a vector of initial conditions
        std::vector<double> starter_conditions;
        starter_conditions.push_back(1.0);
        starter_conditions.push_back(1.0);
        p_delta_notch_interior_srn_model->SetInitialConditions(starter_conditions);
        p_cell_srn_model->SetInteriorSrnModel(p_delta_notch_interior_srn_model);

        UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_cell_srn_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->InitialiseSrnModel();
        std::vector<CellPtr> cells;
        cells.push_back(p_cell);
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Merge nodes 3 and 4
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(3), vertex_mesh.GetNode(4));

        //Update population srns
        VertexElementMap element_map(1);
        element_map.ResetToIdentity();
        VertexBasedPopulationSrn<2>* population_srn(&cell_population.rGetVertexBasedPopulationSrn());
        population_srn->UpdateSrnAfterBirthOrDeath(element_map);

        const unsigned int n_edge_srns = p_cell_srn_model->GetNumEdgeSrn();
        TS_ASSERT_EQUALS(n_edge_srns, 4u);
        TS_ASSERT_DELTA(p_delta_notch_interior_srn_model->GetDelta(), 1.5,1e-6);
        TS_ASSERT_DELTA(p_delta_notch_interior_srn_model->GetNotch(), 1.5,1e-6);
        auto p_delta_notch_edge_0
        = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(static_cast<SrnCellModel*>(p_cell_srn_model)->GetEdgeSrn(0));
        auto p_delta_notch_edge_1
        = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(static_cast<SrnCellModel*>(p_cell_srn_model)->GetEdgeSrn(1));
        auto p_delta_notch_edge_2
        = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(static_cast<SrnCellModel*>(p_cell_srn_model)->GetEdgeSrn(2));
        auto p_delta_notch_edge_3
        = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(static_cast<SrnCellModel*>(p_cell_srn_model)->GetEdgeSrn(3));
        TS_ASSERT_DELTA(p_delta_notch_edge_0->GetDelta(), 1.25,1e-6);
        TS_ASSERT_DELTA(p_delta_notch_edge_0->GetNotch(), 1.25,1e-6);
        TS_ASSERT_DELTA(p_delta_notch_edge_1->GetDelta(), 1.25,1e-6);
        TS_ASSERT_DELTA(p_delta_notch_edge_1->GetNotch(), 1.25,1e-6);

        TS_ASSERT_DELTA(p_delta_notch_edge_2->GetDelta(), 1.0,1e-6);
        TS_ASSERT_DELTA(p_delta_notch_edge_2->GetNotch(), 1.0,1e-6);
        TS_ASSERT_DELTA(p_delta_notch_edge_3->GetDelta(), 1.0,1e-6);
        TS_ASSERT_DELTA(p_delta_notch_edge_3->GetNotch(), 1.0,1e-6);

    }

    void TestPerformT1SwapAndIdentifySwapTypeWithSrn()
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap of the two central nodes is correctly implemented.
         *  _____
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         */
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true,  0.0, 0.0));
        nodes.push_back(new Node<2>(1, true,  1.0, 0.0));
        nodes.push_back(new Node<2>(2, true,  1.0, 1.0));
        nodes.push_back(new Node<2>(3, true,  0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.5, 0.4));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.6));

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {2, 3, 5};
        unsigned node_indices_elem_1[4] = {2, 5, 4, 1};
        unsigned node_indices_elem_2[3] = {1, 4, 0};
        unsigned node_indices_elem_3[4] = {0, 4, 5, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
                nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            }
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        const std::vector<double> init_edge_vars(2, 1.0);
        const std::vector<double> init_interior_vars(2, 1.0);
        for (unsigned elem_index=0; elem_index < vertex_mesh.GetNumElements(); elem_index++)
        {
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

            /* Initialise edge based SRN */
            auto p_element = vertex_mesh.GetElement(elem_index);

            auto p_cell_srn_model = new SrnCellModel();

            /* Gets the edges of the element and create an SRN for each edge */
            for (unsigned i = 0; i < p_element->GetNumEdges(); i++)
            {
                boost::shared_ptr<DeltaNotchSrnEdgeModel> p_edge_model(new DeltaNotchSrnEdgeModel());
                p_edge_model->SetInitialConditions(init_edge_vars);

                p_cell_srn_model->AddEdgeSrnModel(p_edge_model);
            }
            MAKE_PTR(DeltaNotchSrnInteriorModel, p_interior_srn_model);
            p_interior_srn_model->SetInitialConditions(init_interior_vars);
            p_cell_srn_model->SetInteriorSrnModel(p_interior_srn_model);

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->SetBirthTime(0);
            p_cell->InitialiseSrnModel();
            cells.push_back(p_cell);
        }
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);
        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);

        // Perform a T1 swap on nodes 4 and 5
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));

        //Update population srns
        VertexElementMap element_map(4);
        element_map.ResetToIdentity();
        VertexBasedPopulationSrn<2>* population_srn(&cell_population.rGetVertexBasedPopulationSrn());
        population_srn->UpdateSrnAfterBirthOrDeath(element_map);
        //Testing if the SRN quantities have been updated correctly
        //Cells 0 and 2 have new edges
        // New edge in cell 0
        {
            CellPtr cell = cell_population.GetCellUsingLocationIndex(0);
            auto p_cell_model = static_cast<SrnCellModel*>(cell->GetSrnModel());
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_interior_model
            = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_model->GetInteriorSrn());
            const unsigned int n_edge_srns = p_cell_model->GetNumEdgeSrn();
            std::vector<boost::shared_ptr<DeltaNotchSrnEdgeModel> > edge_srn_models(n_edge_srns);
            TS_ASSERT_EQUALS(n_edge_srns, 4u);
            for (unsigned int i=0; i<n_edge_srns; ++i)
            {
                edge_srn_models[i]
                                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
            }
            TS_ASSERT_DELTA(edge_srn_models[0]->GetDelta(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[0]->GetNotch(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetDelta(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetNotch(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetDelta(), 0.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetNotch(), 0.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[3]->GetDelta(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[3]->GetNotch(), 1.0,1e-6);
        }
        // New edge in cell 2
        {
            CellPtr cell = cell_population.GetCellUsingLocationIndex(2);
            auto p_cell_model = static_cast<SrnCellModel*>(cell->GetSrnModel());
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_interior_model
            = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_model->GetInteriorSrn());
            const unsigned int n_edge_srns = p_cell_model->GetNumEdgeSrn();
            std::vector<boost::shared_ptr<DeltaNotchSrnEdgeModel> > edge_srn_models(n_edge_srns);
            TS_ASSERT_EQUALS(n_edge_srns, 4u);
            for (unsigned int i=0; i<n_edge_srns; ++i)
            {
                edge_srn_models[i]
                                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
            }
            TS_ASSERT_DELTA(edge_srn_models[0]->GetDelta(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[0]->GetNotch(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetDelta(), 0.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetNotch(), 0.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetDelta(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetNotch(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[3]->GetDelta(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[3]->GetNotch(), 1.0,1e-6);
        }

        //Edge shrinkage in cells 1 and 3
        {
            CellPtr cell = cell_population.GetCellUsingLocationIndex(1);
            auto p_cell_model = static_cast<SrnCellModel*>(cell->GetSrnModel());
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_interior_model
            = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_model->GetInteriorSrn());
            const unsigned int n_edge_srns = p_cell_model->GetNumEdgeSrn();
            std::vector<boost::shared_ptr<DeltaNotchSrnEdgeModel> > edge_srn_models(n_edge_srns);
            TS_ASSERT_EQUALS(n_edge_srns, 3u);
            for (unsigned int i=0; i<n_edge_srns; ++i)
            {
                edge_srn_models[i]
                                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
            }
            TS_ASSERT_DELTA(edge_srn_models[0]->GetDelta(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[0]->GetNotch(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetDelta(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetNotch(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetDelta(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetNotch(), 1.0,1e-6);

            TS_ASSERT_DELTA(p_interior_model->GetDelta(), 1.5,1e-6);
            TS_ASSERT_DELTA(p_interior_model->GetNotch(), 1.5,1e-6);
        }

        {
            CellPtr cell = cell_population.GetCellUsingLocationIndex(3);
            auto p_cell_model = static_cast<SrnCellModel*>(cell->GetSrnModel());
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_interior_model
            = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_model->GetInteriorSrn());
            const unsigned int n_edge_srns = p_cell_model->GetNumEdgeSrn();
            std::vector<boost::shared_ptr<DeltaNotchSrnEdgeModel> > edge_srn_models(n_edge_srns);
            TS_ASSERT_EQUALS(n_edge_srns, 3u);
            for (unsigned int i=0; i<n_edge_srns; ++i)
            {
                edge_srn_models[i]
                                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
            }
            TS_ASSERT_DELTA(edge_srn_models[0]->GetDelta(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[0]->GetNotch(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetDelta(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetNotch(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetDelta(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetNotch(), 1.0,1e-6);

            TS_ASSERT_DELTA(p_interior_model->GetDelta(), 1.5,1e-6);
            TS_ASSERT_DELTA(p_interior_model->GetNotch(), 1.5,1e-6);
        }
    }

    void TestPerformT2SwapWithSrn()
    {
        /*
         * Create a mesh comprising six nodes contained in three trapezium element and
         * a central triangle element, as shown below. We will test that a T2 swap
         * correctly removes the triangle element from the mesh.
         *
         *      /|\
         *     / | \
         *    /  |  \    (the triangular element has index zero)
         *   /2 /_\ 1\
         *  /  /   \  \
         * /__/__3__\__\
         */
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.5, 0.5));
        nodes.push_back(new Node<2>(3, false, 0.4, 0.25));
        nodes.push_back(new Node<2>(4, false, 0.6, 0.25));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.3));

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {3, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 5, 4};
        unsigned node_indices_elem_2[4] = {2, 0, 3, 5};
        unsigned node_indices_elem_3[4] = {0, 1, 4, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            }
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        const std::vector<double> init_edge_vars(2, 1.0);
        const std::vector<double> init_interior_vars(2, 1.0);
        for (unsigned elem_index=0; elem_index < vertex_mesh.GetNumElements(); elem_index++)
        {
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

            /* Initialise edge based SRN */
            auto p_element = vertex_mesh.GetElement(elem_index);

            auto p_cell_srn_model = new SrnCellModel();

            /* Gets the edges of the element and create an SRN for each edge */
            for (unsigned i = 0; i < p_element->GetNumEdges(); i++)
            {
                boost::shared_ptr<DeltaNotchSrnEdgeModel> p_edge_model(new DeltaNotchSrnEdgeModel());
                p_edge_model->SetInitialConditions(init_edge_vars);

                p_cell_srn_model->AddEdgeSrnModel(p_edge_model);
            }
            MAKE_PTR(DeltaNotchSrnInteriorModel, p_interior_srn_model);
            p_interior_srn_model->SetInitialConditions(init_interior_vars);
            p_cell_srn_model->SetInteriorSrnModel(p_interior_srn_model);

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->SetBirthTime(0);
            p_cell->InitialiseSrnModel();
            cells.push_back(p_cell);
        }
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Perform a T2 swap on the central triangle element
        VertexElement<2,2>* p_element_0 = vertex_mesh.GetElement(0);
        c_vector<double, 2> centroid_of_element_0_before_swap = vertex_mesh.GetCentroidOfElement(0);
        vertex_mesh.PerformT2Swap(*p_element_0);

        //Update population srns
        VertexElementMap element_map(4);
        element_map.SetDeleted(0);
        cell_population.GetCellUsingLocationIndex(0)->Kill();
        cell_population.RemoveDeadCells();
        vertex_mesh.RemoveDeletedNodesAndElements(element_map);

        std::map<Cell*, unsigned> old_map = cell_population.mCellLocationMap;

        cell_population.mCellLocationMap.clear();
        cell_population.mLocationCellMap.clear();

        for (std::list<CellPtr>::iterator cell_iter = cell_population.mCells.begin();
                cell_iter != cell_population.mCells.end();
                ++cell_iter)
        {
            // The cell vector should only ever contain living cells
            unsigned old_elem_index = old_map[(*cell_iter).get()];
            assert(!element_map.IsDeleted(old_elem_index));

            unsigned new_elem_index = element_map.GetNewIndex(old_elem_index);
            cell_population.SetCellUsingLocationIndex(new_elem_index, *cell_iter);
        }

        VertexBasedPopulationSrn<2>* population_srn(&cell_population.rGetVertexBasedPopulationSrn());
        population_srn->UpdateSrnAfterBirthOrDeath(element_map);
        //Testing if the SRN quantities have been updated correctly
        for (unsigned int i=0; i<3; ++i)
        {
            CellPtr cell = cell_population.GetCellUsingLocationIndex(i);
            auto p_cell_model = static_cast<SrnCellModel*>(cell->GetSrnModel());
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_interior_model
            = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_model->GetInteriorSrn());
            const unsigned int n_edge_srns = p_cell_model->GetNumEdgeSrn();
            std::vector<boost::shared_ptr<DeltaNotchSrnEdgeModel> > edge_srn_models(n_edge_srns);
            TS_ASSERT_EQUALS(n_edge_srns, 3u);
            for (unsigned int i=0; i<n_edge_srns; ++i)
            {
                edge_srn_models[i]
                                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
            }
            TS_ASSERT_DELTA(edge_srn_models[0]->GetDelta(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[0]->GetNotch(), 1.0,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetDelta(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetNotch(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetDelta(), 1.25,1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetNotch(), 1.25,1e-6);

            TS_ASSERT_DELTA(p_interior_model->GetDelta(), 1.5,1e-6);
            TS_ASSERT_DELTA(p_interior_model->GetNotch(), 1.5,1e-6);
        }
    }

    void TestPerformT3SwapWithSrn()
    {
        /*
         * Create a mesh comprising 13 nodes containing in five elements, as shown below.
         * We will test that a T3 swap is correctly performed.
         *         _____
         *    |\  |     |  /|
         *    | \ |     | /_|
         *    | / |     | \ |
         *    |/__|_____|__\|
         *    |             |
         *    |_____________|
         */
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true,   0.0,  0.0));
        nodes.push_back(new Node<2>(1, true,   1.0,  0.0));
        nodes.push_back(new Node<2>(2, true,   1.0,  1.0));
        nodes.push_back(new Node<2>(3, true,   0.0,  1.0));
        nodes.push_back(new Node<2>(4, true,   2.0,  0.0));
        nodes.push_back(new Node<2>(5, true,   2.0,  1.0));
        nodes.push_back(new Node<2>(6, true,   1.1,  0.5));
        nodes.push_back(new Node<2>(7, true,  -1.0,  0.0));
        nodes.push_back(new Node<2>(8, true,  -0.1,  0.5));
        nodes.push_back(new Node<2>(9, true,  -1.0,  1.0));
        nodes.push_back(new Node<2>(10, true, -1.0, -1.0));
        nodes.push_back(new Node<2>(11, true,  2.0, -1.0));
        nodes.push_back(new Node<2>(12, true,  2.0,  0.5));

        std::vector<Node<2>*> nodes_in_element0, nodes_in_element1, nodes_in_element2, nodes_in_element3, nodes_in_element4;
        unsigned node_indices_element_0[4] = {0, 1, 2, 3};
        unsigned node_indices_element_1[3] = {4, 12, 6};
        unsigned node_indices_element_2[3] = {12, 5, 6};
        unsigned node_indices_element_3[3] = {7, 8, 9};
        unsigned node_indices_element_4[6] = {10, 11, 4, 1, 0, 7};
        for (unsigned i=0; i<6; i++)
        {
            if (i < 4)
            {
                nodes_in_element0.push_back(nodes[node_indices_element_0[i]]);
            }
            if (i < 3)
            {
                nodes_in_element1.push_back(nodes[node_indices_element_1[i]]);
                nodes_in_element2.push_back(nodes[node_indices_element_2[i]]);
                nodes_in_element3.push_back(nodes[node_indices_element_3[i]]);
            }
            nodes_in_element4.push_back(nodes[node_indices_element_4[i]]);
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_in_element0));
        elements.push_back(new VertexElement<2,2>(1, nodes_in_element1));
        elements.push_back(new VertexElement<2,2>(2, nodes_in_element2));
        elements.push_back(new VertexElement<2,2>(3, nodes_in_element3));
        elements.push_back(new VertexElement<2,2>(4, nodes_in_element4));

        MutableVertexMesh<2,2> mesh(nodes, elements);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        const std::vector<double> init_edge_vars(2, 1.0);
        const std::vector<double> init_interior_vars(2, 1.0);
        for (unsigned elem_index=0; elem_index < mesh.GetNumElements(); elem_index++)
        {
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

            /* Initialise edge based SRN */
            auto p_element = mesh.GetElement(elem_index);

            auto p_cell_srn_model = new SrnCellModel();

            /* Gets the edges of the element and create an SRN for each edge */
            for (unsigned i = 0; i < p_element->GetNumEdges(); i++)
            {
                boost::shared_ptr<DeltaNotchSrnEdgeModel> p_edge_model(new DeltaNotchSrnEdgeModel());
                p_edge_model->SetInitialConditions(init_edge_vars);

                p_cell_srn_model->AddEdgeSrnModel(p_edge_model);
            }
            MAKE_PTR(DeltaNotchSrnInteriorModel, p_interior_srn_model);
            p_interior_srn_model->SetInitialConditions(init_interior_vars);
            p_cell_srn_model->SetInteriorSrnModel(p_interior_srn_model);

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->SetBirthTime(0);
            p_cell->InitialiseSrnModel();
            cells.push_back(p_cell);
        }
        VertexBasedCellPopulation<2> cell_population(mesh, cells);

        // Set the threshold distance between vertices for a T3 swap as follows, to ease calculations
        mesh.SetCellRearrangementThreshold(0.1*1.0/1.5);

        // Move node 6 to the left so that it overlaps element 1
        ChastePoint<2> point = mesh.GetNode(6)->GetPoint();
        point.SetCoordinate(0u, 0.9);
        mesh.SetNode(6, point);

        // Move node 8 to the left so that it overlaps element 1
        point.SetCoordinate(0u, 0.1);
        mesh.SetNode(8, point);

        // Call method to update mesh in this situation
        mesh.ReMesh();

        // Test if the swap has been recorded properly
        auto operation_recorder = mesh.GetOperationRecorder();
        std::vector<EdgeOperation*> edge_operations = operation_recorder->GetEdgeOperations();
        const unsigned int n_operations = edge_operations.size();
        //Two node merging operations in two elements and two new edge operations in the other two elements
        TS_ASSERT_EQUALS(n_operations, 8u);
        unsigned n_edge_splits= 0, n_new_edges= 0;
        std::vector<std::vector<unsigned int> > element_to_operations(5);
        for (unsigned int i=0; i<n_operations; ++i)
        {
            if (edge_operations[i]->GetOperation() == EDGE_OPERATION_SPLIT)
                n_edge_splits++;
            if (edge_operations[i]->GetOperation() == EDGE_OPERATION_ADD)
                n_new_edges++;
            //Determine operations that an element underwent
            const unsigned int elem_index = edge_operations[i]->GetElementIndex();
            element_to_operations[elem_index].push_back(edge_operations[i]->GetOperation());
        }
        TS_ASSERT_EQUALS(n_edge_splits, 5u);
        TS_ASSERT_EQUALS(n_new_edges, 3);

        //Update population srns
        VertexElementMap element_map(5);
        element_map.ResetToIdentity();
        VertexBasedPopulationSrn<2>* population_srn(&cell_population.rGetVertexBasedPopulationSrn());
        population_srn->UpdateSrnAfterBirthOrDeath(element_map);

        for (unsigned int elem_index = 0; elem_index<mesh.GetNumElements(); ++elem_index)
        {
            CellPtr cell = cell_population.GetCellUsingLocationIndex(elem_index);
            auto p_cell_model = static_cast<SrnCellModel*>(cell->GetSrnModel());
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_interior_model
            = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_model->GetInteriorSrn());
            const unsigned int n_edge_srns = p_cell_model->GetNumEdgeSrn();
            std::vector<boost::shared_ptr<DeltaNotchSrnEdgeModel> > edge_srn_models(n_edge_srns);
            unsigned int n_splits= 0;
            switch (elem_index)
            {
            case 0:
                TS_ASSERT_EQUALS(element_to_operations[elem_index].size(), 5u);
                //all edge splits occured in element 0
                for (unsigned int i=0; i<5; ++i)
                {
                    if (element_to_operations[elem_index][i]==EDGE_OPERATION_SPLIT)
                        n_splits++;
                }
                TS_ASSERT_EQUALS(n_splits, 5u);
                TS_ASSERT_DELTA(p_interior_model->GetDelta(), 1.0,1e-6);
                TS_ASSERT_DELTA(p_interior_model->GetNotch(), 1.0,1e-6);
                for (unsigned int i=0; i<n_edge_srns; ++i)
                {
                    edge_srn_models[i]
                                    = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
                }
                TS_ASSERT_DELTA(edge_srn_models[0]->GetDelta(), 1.0,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[0]->GetNotch(), 1.0,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[1]->GetDelta(), 0.4,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[1]->GetNotch(), 0.4,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[2]->GetDelta(), 0.1,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[2]->GetNotch(), 0.1,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[3]->GetDelta(), 0.1,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[3]->GetNotch(), 0.1,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[4]->GetDelta(), 0.4,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[4]->GetNotch(), 0.4,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[5]->GetDelta(), 1.0,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[5]->GetNotch(), 1.0,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[6]->GetDelta(), 0.45,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[6]->GetNotch(), 0.45,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[7]->GetDelta(), 0.1,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[7]->GetNotch(), 0.1,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[8]->GetDelta(), 0.45,1e-4);
                TS_ASSERT_DELTA(edge_srn_models[8]->GetNotch(), 0.45,1e-4);
                break;
            case 1:
                //New edge has been added to element_1
                TS_ASSERT_EQUALS(element_to_operations[elem_index].size(), 1u);
                TS_ASSERT_EQUALS(element_to_operations[elem_index][0], EDGE_OPERATION_ADD);
                TS_ASSERT_DELTA(p_interior_model->GetDelta(), 1.0,1e-6);
                TS_ASSERT_DELTA(p_interior_model->GetNotch(), 1.0,1e-6);
                for (unsigned int i=0; i<n_edge_srns; ++i)
                {
                    edge_srn_models[i]
                                    = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
                    if (i!=2)
                    {
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetDelta(), 1.0,1e-6);
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetNotch(), 1.0,1e-6);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetDelta(), 0.0,1e-6);
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetNotch(), 0.0,1e-6);
                    }
                }
                break;
            case 2:
                //New edge has been added to element_2
                TS_ASSERT_EQUALS(element_to_operations[elem_index].size(), 1u);
                TS_ASSERT_EQUALS(element_to_operations[elem_index][0], EDGE_OPERATION_ADD);
                TS_ASSERT_DELTA(p_interior_model->GetDelta(), 1.0,1e-6);
                TS_ASSERT_DELTA(p_interior_model->GetNotch(), 1.0,1e-6);
                for (unsigned int i=0; i<n_edge_srns; ++i)
                {
                    edge_srn_models[i]
                                    = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
                    if (i!=2)
                    {
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetDelta(), 1.0,1e-6);
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetNotch(), 1.0,1e-6);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetDelta(), 0.0,1e-6);
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetNotch(), 0.0,1e-6);
                    }
                }
                break;
            case 3:
                //New edge has been added to element_3
                TS_ASSERT_EQUALS(element_to_operations[elem_index].size(), 1u);
                TS_ASSERT_EQUALS(element_to_operations[elem_index][0], EDGE_OPERATION_ADD);
                TS_ASSERT_DELTA(p_interior_model->GetDelta(), 1.0,1e-6);
                TS_ASSERT_DELTA(p_interior_model->GetNotch(), 1.0,1e-6);
                for (unsigned int i=0; i<n_edge_srns; ++i)
                {
                    edge_srn_models[i]
                                    = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
                    if (i!=1)
                    {
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetDelta(), 1.0,1e-6);
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetNotch(), 1.0,1e-6);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetDelta(), 0.0,1e-6);
                        TS_ASSERT_DELTA(edge_srn_models[i]->GetNotch(), 0.0,1e-6);
                    }
                }
                break;
            case 4:
                //element_4 has not been modified
                TS_ASSERT_EQUALS(element_to_operations[elem_index].size(), 0u);
                break;
            }
        }

    }

    void TestDivideElementWithSrn()
    {
        //This test is based on TestAddCellWithSimpleMesh() in TestVertexBasedCellPopulation
        // Make some nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 2.0, -1.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, true, -2.0, 1.0));
        nodes.push_back(new Node<2>(3, true, -2.0, -1.0));
        nodes.push_back(new Node<2>(4, true, 0.0, 2.0));

        // Make a rectangular element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);

        // Make a triangular element out of nodes 1,4,2
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[2]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_2));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Create cells
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        std::vector<CellPtr> cells;
        for (unsigned elem_index=0; elem_index < vertex_mesh.GetNumElements(); elem_index++)
        {
            /* Initalise cell cycle */
            auto p_cc_model = new AlwaysDivideCellCycleModel();
            p_cc_model->SetDimension(2);

            /* Initialise edge based SRN */
            auto p_element = vertex_mesh.GetElement(elem_index);
            auto p_cell_edge_srn_model = new SrnCellModel();

            /* Gets the edges of the element and create an SRN for each edge */
            for (unsigned i = 0; i < p_element->GetNumEdges() ; i ++)
            {
                std::vector<double> initial_conditions;

                /* Initial concentration of delta and notch is the same */
                initial_conditions.push_back(2.0);
                initial_conditions.push_back(2.0);

                MAKE_PTR(DeltaNotchSrnEdgeModel, p_srn_model);
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
            }
            MAKE_PTR(DeltaNotchSrnInteriorModel, p_interior_srn);
            p_interior_srn->SetInitialConditions(std::vector<double>(2, 4.0));
            p_cell_edge_srn_model->SetInteriorSrnModel(p_interior_srn);

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->SetBirthTime(0.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();
            cells.push_back(p_cell);
        }
        // Create cell population
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        unsigned old_num_nodes = vertex_mesh.GetNumNodes();

        /* Create an edge tracking modifier. DeltaNotch models require information obtained from
           CellData and CellEdgeData */
        MAKE_PTR(DeltaNotchEdgeInteriorTrackingModifier<2>, p_modifier);
        p_modifier->SetupSolve(cell_population,"TestVertexMeshCellDivisionWithSrn");
        /*Divide the 0th cell*/
        {
            auto p_cell = cell_population.GetCellUsingLocationIndex(0);
            p_cell->ReadyToDivide();
            auto p_new_cell = p_cell->Divide();
            cell_population.AddCell(p_new_cell, p_cell);
        }

        // Check the location of the new nodes
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes)->rGetLocation()[1], 1.0, 1e-12);

        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes+1)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes+1)->rGetLocation()[1], -1.0, 1e-12);

        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(3), 6u);

        TS_ASSERT_EQUALS(cell_population.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(cell_population.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(cell_population.GetElement(1)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(cell_population.GetElement(1)->GetNodeGlobalIndex(3), 5u);

        TS_ASSERT_EQUALS(cell_population.GetElement(2)->GetNodeGlobalIndex(0), 5u);
        TS_ASSERT_EQUALS(cell_population.GetElement(2)->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(cell_population.GetElement(2)->GetNodeGlobalIndex(2), 3u);
        TS_ASSERT_EQUALS(cell_population.GetElement(2)->GetNodeGlobalIndex(3), 6u);

        cell_population.Update(true);

        // Test if interior/edge quantities are update properly
        // Check the 0th cell and its edges//
        {
            auto p_cell = cell_population.GetCellUsingLocationIndex(0);
            auto p_cell_model = static_cast<SrnCellModel*>(p_cell->GetSrnModel());
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_interior_model
            = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_model->GetInteriorSrn());
            const unsigned int n_edge_srns = p_cell_model->GetNumEdgeSrn();
            std::vector<boost::shared_ptr<DeltaNotchSrnEdgeModel> > edge_srn_models(n_edge_srns);
            TS_ASSERT_EQUALS(n_edge_srns, 4);

            for (unsigned int i=0; i<n_edge_srns; ++i)
            {
                edge_srn_models[i]
                                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
            }

            TS_ASSERT_DELTA(p_interior_model->GetDelta(), 2.0, 1e-6);
            TS_ASSERT_DELTA(p_interior_model->GetNotch(), 2.0, 1e-6);

            TS_ASSERT_DELTA(edge_srn_models[0]->GetDelta(), 2.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[0]->GetNotch(), 2.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetDelta(), 1.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetNotch(), 1.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetDelta(), 0.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetNotch(), 0.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[3]->GetDelta(), 1.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[3]->GetNotch(), 1.0, 1e-6);
        }

        // Check the 2nd cell //
        {
            auto p_cell = cell_population.GetCellUsingLocationIndex(2);
            auto p_cell_model = static_cast<SrnCellModel*>(p_cell->GetSrnModel());
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_interior_model
            = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_model->GetInteriorSrn());
            const unsigned int n_edge_srns = p_cell_model->GetNumEdgeSrn();
            std::vector<boost::shared_ptr<DeltaNotchSrnEdgeModel> > edge_srn_models(n_edge_srns);
            TS_ASSERT_EQUALS(n_edge_srns, 4);

            for (unsigned int i=0; i<n_edge_srns; ++i)
            {
                edge_srn_models[i]
                                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
            }

            TS_ASSERT_DELTA(p_interior_model->GetDelta(), 2.0, 1e-6);
            TS_ASSERT_DELTA(p_interior_model->GetNotch(), 2.0, 1e-6);

            TS_ASSERT_DELTA(edge_srn_models[0]->GetDelta(), 1.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[0]->GetNotch(), 1.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetDelta(), 2.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetNotch(), 2.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetDelta(), 1.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetNotch(), 1.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[3]->GetDelta(), 0.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[3]->GetNotch(), 0.0, 1e-6);
        }

        // Check the 1st cell //
        {
            auto p_cell = cell_population.GetCellUsingLocationIndex(1);
            auto p_cell_model = static_cast<SrnCellModel*>(p_cell->GetSrnModel());
            boost::shared_ptr<DeltaNotchSrnInteriorModel> p_interior_model
            = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_model->GetInteriorSrn());
            const unsigned int n_edge_srns = p_cell_model->GetNumEdgeSrn();
            std::vector<boost::shared_ptr<DeltaNotchSrnEdgeModel> > edge_srn_models(n_edge_srns);
            TS_ASSERT_EQUALS(n_edge_srns, 4);

            for (unsigned int i=0; i<n_edge_srns; ++i)
            {
                edge_srn_models[i]
                                = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_model->GetEdgeSrn(i));
            }

            TS_ASSERT_DELTA(p_interior_model->GetDelta(), 4.0, 1e-6);
            TS_ASSERT_DELTA(p_interior_model->GetNotch(), 4.0, 1e-6);

            TS_ASSERT_DELTA(edge_srn_models[0]->GetDelta(), 2.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[0]->GetNotch(), 2.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetDelta(), 2.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[1]->GetNotch(), 2.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetDelta(), 1.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[2]->GetNotch(), 1.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[3]->GetDelta(), 1.0, 1e-6);
            TS_ASSERT_DELTA(edge_srn_models[3]->GetNotch(), 1.0, 1e-6);
        }
    }


};

#endif /*TESTMUTABLEVERTEXMESHOPERATIONSWITHPOPULATIONSRN_HPP_*/
