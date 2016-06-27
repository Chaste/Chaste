/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTNUMERICALMETHODS_HPP_
#define TESTNUMERICALMETHODS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "CellAncestor.hpp"
#include "CellId.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "PopulationTestingForce.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "Warnings.hpp"


#include "PetscSetupAndFinalize.hpp"

class TestNumericalMethods : public AbstractCellBasedTestSuite
{
public:
    
    void TestUpdateAllNodePositionsWithMeshBased() throw(Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population, with no ghost nodes at the moment
        MeshBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetDampingConstantNormal(1.1);

        // Create a force collection
        std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
        MAKE_PTR(PopulationTestingForce<2>, p_test_force);
        force_collection.push_back(p_test_force);

        // Create numerical method for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);
        
        double dt = 0.01;
        
        p_fe_method->SetCellPopulation(&cell_population);
        p_fe_method->SetForceCollection(&force_collection);

        // Save starting positions
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        for(unsigned j=0; j<cell_population.GetNumNodes(); j++){
            old_posns[j][0] = cell_population.GetNode(j)->rGetLocation()[0];
            old_posns[j][1] = cell_population.GetNode(j)->rGetLocation()[1];
        }
     
        // Update positions and check the answer   
        p_fe_method->UpdateAllNodePositions(dt);

        for(unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            c_vector<double, 2> actualLocation = cell_population.GetNode(j)->rGetLocation();
            
            double damping =  cell_population.GetDampingConstant(j);
            c_vector<double, 2> expectedLocation;
            expectedLocation = p_test_force->GetExpectedOneStepLocationFE(j, damping, old_posns[j], dt);
            
            TS_ASSERT_DELTA(norm_2(actualLocation - expectedLocation), 0, 1e-6);
        }
    }

    void TestUpdateAllNodePositionsWithMeshBasedWithGhosts() throw(Exception)
    {
        HoneycombMeshGenerator generator(3, 3, 1);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        CellPropertyRegistry::Instance()->Clear();
        RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

        // Set up cells
        std::vector<CellPtr> cells;
        cells.clear();
        unsigned num_cells = location_indices.empty() ? p_mesh->GetNumNodes() : location_indices.size();
        cells.reserve(num_cells);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            unsigned generation;
            double y = 0.0;

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            }

            FixedDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel;
            p_cell_cycle_model->SetDimension(2);

            double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
            double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            boost::shared_ptr<AbstractCellProperty> p_stem_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
            boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
            boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            double birth_time = -p_random_num_gen->ranf();

            if (y <= 0.3)
            {
                p_cell->SetCellProliferativeType(p_stem_type);
                generation = 0;
                birth_time *= typical_stem_cycle_time; // hours
            }
            else if (y < 2.0)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 1;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 3.0)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 2;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 4.0)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 3;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else
            {
                if (p_cell_cycle_model->CanCellTerminallyDifferentiate())
                {
                    p_cell->SetCellProliferativeType(p_diff_type);
                }
                else
                {
                    p_cell->SetCellProliferativeType(p_transit_type);
                }
                generation = 4;
                birth_time *= typical_transit_cycle_time; // hours
            }

            p_cell_cycle_model->SetGeneration(generation);
            p_cell->SetBirthTime(birth_time);

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                cells.push_back(p_cell);
            }
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.SetDampingConstantNormal(1.1);

        // Create a force collection
        std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
        MAKE_PTR(PopulationTestingForce<2>, p_test_force);
        force_collection.push_back(p_test_force);

        // Create numerical methods for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);
        
        double dt = 0.01;
       
        p_fe_method->SetCellPopulation(&cell_population);
        p_fe_method->SetForceCollection(&force_collection);

        // Save starting positions
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        for(unsigned j=0; j<cell_population.GetNumNodes(); j++){
            old_posns[j][0] = cell_population.GetNode(j)->rGetLocation()[0];
            old_posns[j][1] = cell_population.GetNode(j)->rGetLocation()[1];
        }
     
        // Update positions   
        p_fe_method->UpdateAllNodePositions(dt);

        //Check the answer (for cell associated nodes only)
        for(AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            int j = cell_population.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, 2> actualLocation = cell_population.GetNode(j)->rGetLocation();
            
            double damping =  cell_population.GetDampingConstant(j);
            c_vector<double, 2> expectedLocation;
            expectedLocation = p_test_force->GetExpectedOneStepLocationFE(j, damping, old_posns[j], dt);
            TS_ASSERT_DELTA(norm_2(actualLocation - expectedLocation), 0, 1e-9);
        }
    }
    void TestUpdateAllNodePositionsWithNodeBased() throw(Exception)
    {
    }

    void TestUpdateAllNodePositionsWithNodeBasedWithParticles() throw(Exception)
    {
        EXIT_IF_PARALLEL;    // This test doesn't work in parallel.

        HoneycombMeshGenerator generator(3, 3, 1);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        MAKE_PTR(NodesOnlyMesh<2>, p_mesh);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        CellPropertyRegistry::Instance()->Clear();
        RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

        // Set up cells
        std::vector<CellPtr> cells;
        cells.clear();
        unsigned num_cells = location_indices.empty() ? p_mesh->GetNumNodes() : location_indices.size();
        cells.reserve(num_cells);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double y = 0.0;
            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            }

            FixedDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel;
            p_cell_cycle_model->SetDimension(2);

            double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
            double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

            unsigned generation;
            if (y <= 0.3)
            {
                generation = 0;
            }
            else if (y < 2.0)
            {
                generation = 1;
            }
            else if (y < 3.0)
            {
                generation = 2;
            }
            else if (y < 4.0)
            {
                generation = 3;
            }
            else
            {
                generation = 4;
            }
            p_cell_cycle_model->SetGeneration(generation);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            boost::shared_ptr<AbstractCellProperty> p_stem_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
            boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
            boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            if (y <= 0.3)
            {
                p_cell->SetCellProliferativeType(p_stem_type);
            }
            else
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                if (y >= 4.0 && p_cell_cycle_model->CanCellTerminallyDifferentiate())
                {
                    p_cell->SetCellProliferativeType(p_diff_type);
                }
            }

            double birth_time = -p_random_num_gen->ranf();
            if (y <= 0.3)
            {
                birth_time *= typical_stem_cycle_time; // hours
            }
            else
            {
                birth_time *= typical_transit_cycle_time; // hours
            }
            p_cell->SetBirthTime(birth_time);

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                cells.push_back(p_cell);
            }
        }

        NodeBasedCellPopulationWithParticles<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.SetDampingConstantNormal(1.1);

        // Create a force collection
        std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
        MAKE_PTR(PopulationTestingForce<2>, p_test_force);
        force_collection.push_back(p_test_force);

        // Create numerical method for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);
        
        double dt = 0.01;
        
        p_fe_method->SetCellPopulation(&cell_population);
        p_fe_method->SetForceCollection(&force_collection);

        // Save starting positions
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        for(unsigned j=0; j<cell_population.GetNumNodes(); j++){
            old_posns[j][0] = cell_population.GetNode(j)->rGetLocation()[0];
            old_posns[j][1] = cell_population.GetNode(j)->rGetLocation()[1];
        }
     
        // Update positions and check the answer   
        p_fe_method->UpdateAllNodePositions(dt);

        for(unsigned j=0; j<cell_population.GetNumNodes(); j++)
        {
            c_vector<double, 2> actualLocation = cell_population.GetNode(j)->rGetLocation();
            
            double damping =  cell_population.GetDampingConstant(j);
            c_vector<double, 2> expectedLocation;
            expectedLocation = p_test_force->GetExpectedOneStepLocationFE(j, damping, old_posns[j], dt);
            
            TS_ASSERT_DELTA(norm_2(actualLocation - expectedLocation), 0, 1e-12);
        }
    }

    void TestUpdateAllNodePositionsWithVertexBased() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(5, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Impose a larger cell rearrangement threshold so that motion is uninhibited (see #1376)
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create a cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.1);

        // Create a force collection
        std::vector<boost::shared_ptr<AbstractForce<2,2> > > force_collection;
        MAKE_PTR(PopulationTestingForce<2>, p_test_force);
        force_collection.push_back(p_test_force);

        // Create numerical methods for testing
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_fe_method);
        
        double dt = 0.01;
        
        p_fe_method->SetCellPopulation(&cell_population);
        p_fe_method->SetForceCollection(&force_collection);

        // Save starting positions
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        for(unsigned j=0; j<cell_population.GetNumNodes(); j++){
            old_posns[j][0] = cell_population.GetNode(j)->rGetLocation()[0];
            old_posns[j][1] = cell_population.GetNode(j)->rGetLocation()[1];
        }
     
        // Update positions and check the answer   
        p_fe_method->UpdateAllNodePositions(dt);

        for(unsigned j=0; j<cell_population.GetNumNodes(); j++){
            
            c_vector<double, 2> actualLocation = cell_population.GetNode(j)->rGetLocation();
            
            double damping =  cell_population.GetDampingConstant(j);
            c_vector<double, 2> expectedLocation;
            expectedLocation = p_test_force->GetExpectedOneStepLocationFE(j, damping, old_posns[j], dt);
            
            TS_ASSERT_DELTA(norm_2(actualLocation - expectedLocation), 0, 1e-12);
        }
    }
};

#endif /*TESTMESHBASEDCELLPOPULATION_HPP_*/
