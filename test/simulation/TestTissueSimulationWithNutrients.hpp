/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "TissueSimulationArchiver.hpp"

#include "TissueSimulationWithNutrients.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "FixedDurationGenerationBasedCellCycleModelCellsGenerator.hpp"
#include "SimpleNutrientPde.hpp"
#include "CellwiseNutrientSinkPde.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ReplicatableVector.hpp"
#include "NumericFileComparison.hpp"

class SimplePdeForTesting : public AbstractLinearEllipticPde<2,2>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>& x)
    {
        return -1.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>& x, Element<2,2>*)
    {
        return 0.0;
    }

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }

};


class TestTissueSimulationWithNutrients : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    /*
     * A two-part test for the PostSolve() method.
     *
     * Firstly, test the PDE solver using the problem del squared C = 0.1
     * on the unit disc, with boundary condition C=1 on r=1, which has
     * analytic solution C = 1-0.25*(1-r^2).
     *
     * Secondly, test that cells' hypoxic durations are correctly updated when a
     * nutrient distribution is prescribed.
     */
    void TestPostSolve() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Change the hypoxic concentration, just for this test
        TissueConfig::Instance()->SetHepaOneCellHypoxicConcentration(0.9);
        TissueConfig::Instance()->SetHepaOneCellQuiescentConcentration(0.9);
        TissueConfig::Instance()->SetHepaOneParameters();

        // Set up mesh
        MutableMesh<2,2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<TissueCell> cells;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel(2));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                    (TissueConfig::Instance()->GetHepaOneCellG1Duration()
                                    +TissueConfig::Instance()->GetSG2MDuration());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Set up tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Set up CellwiseData and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);
        p_data->SetTissue(tissue);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i));
        }

        // Set up PDE
        SimplePdeForTesting pde;

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        // Use an extremely small cutoff so that no cells interact
        // - this is to ensure that in the Solve method, the cells don't move
        // (we need to call Solve to set up the .viznutrient file)
        linear_force.UseCutoffPoint(0.0001);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, force_collection, &pde);
        simulator.SetOutputDirectory("TestPostSolveMethod");
        simulator.SetEndTime(2.0/120.0);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&tissue);
        simulator.AddCellKiller(&killer);

        // Run tissue simulation
        simulator.Solve();

        // Check the correct solution was obtained
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            double radius = norm_2(tissue.GetLocationOfCellCentre(*cell_iter));
            double analytic_solution = 1 - 0.25*(1 - pow(radius,2.0));

            // Get cell model
            AbstractCellCycleModel* p_abstract_model = cell_iter->GetCellCycleModel();
            SimpleOxygenBasedCellCycleModel* p_oxygen_model = static_cast<SimpleOxygenBasedCellCycleModel*> (p_abstract_model);

            // First part of test - check that PDE solver is working correctly
            TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), analytic_solution, 1e-2);

            // Second part of test - check that each cell's hypoxic duration is correctly updated
            if ( p_data->GetValue(*cell_iter) >= TissueConfig::Instance()->GetHepaOneCellHypoxicConcentration() )
            {
                TS_ASSERT_DELTA(p_oxygen_model->GetCurrentHypoxicDuration(), 0.0, 1e-5);
            }
            else
            {
                TS_ASSERT_DELTA(p_oxygen_model->GetCurrentHypoxicDuration(), 2/120.0, 1e-5);
            }
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    void TestWithOxygen() throw(Exception)
    {
        EXIT_IF_PARALLEL; //defined in PetscTools

        TissueConfig::Instance()->SetHepaOneParameters();

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel(2));
            double birth_time = -1.0 - ( (double) i/p_mesh->GetNumNodes() )*
                                    (TissueConfig::Instance()->GetHepaOneCellG1Duration()
                                    +TissueConfig::Instance()->GetSG2MDuration());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Set up tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        TissueConfig::Instance()->SetOutputTissueAreas(true); // record the spheroid radius and apoptotic radius

        // Set up CellwiseData and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(),1);
        p_data->SetTissue(tissue);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i));
        }

        // Set up PDE
        SimpleNutrientPde<2> pde(0.1);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, force_collection, &pde);
        simulator.SetOutputDirectory("TissueSimulationWithOxygen");
        simulator.SetEndTime(0.5);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&tissue);
        simulator.AddCellKiller(&killer);

        // Run tissue simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    /*
     * This test compares the visualizer output from the previous test
     * with a known file.
     *
     * Note: if the previous test is changed we need to update the file
     * this test refers to.
     */
    void TestVisualizerOutput() throw (Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Work out where one of the previous tests wrote its files
        OutputFileHandler handler("TissueSimulationWithOxygen", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";
       
        NumericFileComparison comp_nut(results_dir + "/results.viznutrient", "notforrelease_cell_based/test/data/TissueSimulationWithOxygen/results.viznutrient");
        TS_ASSERT(comp_nut.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.viznutrient notforrelease_cell_based/test/data/TissueSimulationWithOxygen/results.viznutrient").c_str()), 0);

        NumericFileComparison comp_ele(results_dir + "/results.vizelements", "notforrelease_cell_based/test/data/TissueSimulationWithOxygen/results.vizelements");
        TS_ASSERT(comp_ele.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizelements notforrelease_cell_based/test/data/TissueSimulationWithOxygen/results.vizelements").c_str()), 0);

        NumericFileComparison comp_nodes(results_dir + "/results.viznodes", "notforrelease_cell_based/test/data/TissueSimulationWithOxygen/results.viznodes");
        TS_ASSERT(comp_nodes.CompareFiles(1e-15));

        NumericFileComparison comp_celltypes(results_dir + "/results.vizcelltypes", "notforrelease_cell_based/test/data/TissueSimulationWithOxygen/results.vizcelltypes");
        TS_ASSERT(comp_celltypes.CompareFiles(1e-15));

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizsetup notforrelease_cell_based/test/data/TissueSimulationWithOxygen/results.vizsetup").c_str()), 0);
    }


    void TestWithPointwiseNutrientSink() throw(Exception)
    {
        EXIT_IF_PARALLEL; //defined in PetscTools

        TissueConfig::Instance()->SetHepaOneParameters();

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel(2));
            double birth_time = -1.0 - ( (double) i/p_mesh->GetNumNodes() )*
                                    (TissueConfig::Instance()->GetHepaOneCellG1Duration()
                                    +TissueConfig::Instance()->GetSG2MDuration());
            cell.SetBirthTime(birth_time);

            // Make the cell apoptotic if near the centre
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            double y = p_mesh->GetNode(i)->rGetLocation()[1];
            double dist_from_centre = sqrt( (x-2.5)*(x-2.5) + (y-2.5)*(y-2.5) );
            if (dist_from_centre < 1.5)
            {
                cell.SetCellProliferativeType(APOPTOTIC);
            }

            cells.push_back(cell);
        }

        // Set up tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        TissueConfig::Instance()->SetOutputTissueAreas(true); // record the spheroid radius and apoptotic radius

        // Set up CellwiseData and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(),1);
        p_data->SetTissue(tissue);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i));
        }

        // Set up PDE
        CellwiseNutrientSinkPde<2> pde(tissue, 0.1);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, force_collection, &pde);
        simulator.SetOutputDirectory("TissueSimulationWithPointwiseNutrientSink");
        simulator.SetEndTime(0.5);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&tissue);
        simulator.AddCellKiller(&killer);

        // Run tissue simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // A few hardcoded tests to check nothing has changed
        std::vector<double> node_5_location = simulator.GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.6576, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 1.1358, 1e-4);
        TS_ASSERT_DELTA(p_data->GetValue(simulator.rGetTissue().rGetCellUsingLocationIndex(5)), 0.9702, 1e-4);

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    void TestSpheroidStatistics() throw (Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        TissueConfig::Instance()->SetHepaOneParameters();

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel(2));
            cell.SetBirthTime(-0.1);

            // Label three neighbouring cells as apoptotic
            if (i==12 || i==13 || i==17)
            {
                cell.SetCellProliferativeType(APOPTOTIC);
            }
            cells.push_back(cell);
        }

        // Set up tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        TissueConfig::Instance()->SetOutputTissueAreas(true); // record the spheroid radius and apoptotic radius

        // Set up CellwiseData and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(),1);
        p_data->SetTissue(tissue);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i));
        }

        // Set up PDE
        SimpleNutrientPde<2> pde(0.1);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, force_collection, &pde);
        simulator.SetOutputDirectory("TestSpheroidStatistics");
        simulator.SetEndTime(1.0/120.0);
        simulator.SetWriteAverageRadialNutrientResults(5);

        // Add an oxygen-dependent cell killer to the tissue simulation
        OxygenBasedCellKiller<2> killer(&tissue);
        simulator.AddCellKiller(&killer);

        // Run the tissue simulation for one timestep
        simulator.Solve();

        // Just check that we do indeed have three apoptotic cells
        unsigned num_apoptotic_cells = 0;
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            if (cell_iter->GetCellProliferativeType()==APOPTOTIC)
            {
                num_apoptotic_cells++;
            }
        }
        TS_ASSERT_EQUALS(num_apoptotic_cells, 3u);

        /**
         * We have 25 cells. Adding up the boundary cell areas, we
         * should have the equivalent area of 16 full regular hexagonal
         * cells.
         *
         * The area of a single hexagonal cell is sqrt(3)/2, so the
         * correct spheroid radius is given by sqrt((16*sqrt(3)/2)/pi).
         *
         * Since there are 3 apoptotic cells, the correct apoptotic radius
         * is given by sqrt((3*sqrt(3)/2)/pi).
         */

        // Work out where the previous test wrote its files
        OutputFileHandler handler("TestSpheroidStatistics", false);
        std::string areas_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/tissueareas.dat";
        TS_ASSERT_EQUALS(system(("diff " + areas_results_file + " notforrelease_cell_based/test/data/TestSpheroidStatistics/Areas.dat").c_str()), 0);

        std::string dist_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/radial_dist.dat";
        TS_ASSERT_EQUALS(system(("diff " + dist_results_file + " notforrelease_cell_based/test/data/TestSpheroidStatistics/radial_dist.dat").c_str()), 0);

        // Coverage
        TS_ASSERT_THROWS_NOTHING(simulator.WriteAverageRadialNutrientDistribution(SimulationTime::Instance()->GetTime(),5));

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    void TestCoarseNutrientMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        TissueConfig::Instance()->SetHepaOneParameters();

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel(2));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                    (TissueConfig::Instance()->GetHepaOneCellG1Duration()
                                    +TissueConfig::Instance()->GetSG2MDuration());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Set up tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Set up CellwiseData and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(),1);
        p_data->SetTissue(tissue);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i));
        }

        // Set up PDE
        AveragedSinksPde<2> pde(tissue, -0.1);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, force_collection, NULL, &pde);
        simulator.SetOutputDirectory("TestCoarseNutrientMesh");
        simulator.SetEndTime(0.05);

        // Coverage
        simulator.SetAveragedSinksPde(&pde);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&tissue);
        simulator.AddCellKiller(&killer);

        // Test creation of mpCoarseNutrientMesh
        simulator.UseCoarseNutrientMesh(10.0);

        // Find centre of tissue
        c_vector<double,2> centre_of_tissue = zero_vector<double>(2);

        for (unsigned i=0; i<simulator.rGetTissue().GetNumNodes(); i++)
        {
            centre_of_tissue += simulator.rGetTissue().GetNode(i)->rGetLocation();
        }
        centre_of_tissue /= simulator.rGetTissue().GetNumNodes();

        // Find centre of nutrient mesh
        c_vector<double,2> centre_of_nutrient_mesh = zero_vector<double>(2);

        for (unsigned i=0; i<simulator.mpCoarseNutrientMesh->GetNumNodes(); i++)
        {
            centre_of_nutrient_mesh += simulator.mpCoarseNutrientMesh->GetNode(i)->rGetLocation();
        }
        centre_of_nutrient_mesh /= simulator.mpCoarseNutrientMesh->GetNumNodes();

        // Test that the two centres match
        TS_ASSERT_DELTA(centre_of_tissue[0], centre_of_nutrient_mesh[0], 1e-4);
        TS_ASSERT_DELTA(centre_of_tissue[1], centre_of_nutrient_mesh[1], 1e-4);

        // Test FindElementContainingCell and initialisation of mCellNutrientElementMap

        simulator.InitialiseCoarseNutrientMesh(); // coverage

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
            cell_iter != tissue.End();
            ++cell_iter)
        {
            unsigned containing_element_index = simulator.mCellNutrientElementMap[&(*cell_iter)];
            TS_ASSERT(containing_element_index < simulator.mpCoarseNutrientMesh->GetNumElements());
            TS_ASSERT_EQUALS(containing_element_index, simulator.FindElementContainingCell(*cell_iter));
        }

        // Run tissue simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        TS_ASSERT(simulator.mpCoarseNutrientMesh != NULL);

        ReplicatableVector nutrient_conc(simulator.GetNutrientSolution());

        // Test the nutrient concentration at the coarse mesh nodes is
        // equal to 1.0 if the nodes are away from the cells
        for (unsigned i=0; i<nutrient_conc.GetSize(); i++)
        {
            c_vector<double,2> centre;
            centre(0) = 2.5; // assuming 5 by 5 honeycomb mesh
            centre(1) = 2.5;
            c_vector<double,2> posn = simulator.mpCoarseNutrientMesh->GetNode(i)->rGetLocation();
            double dist = norm_2(centre-posn);
            double u = nutrient_conc[i];

            if (dist > 4.0)
            {
                TS_ASSERT_DELTA(u, 1.0, 1e-5);
            }
        }

        // Loop over cells, find the coarse mesh element containing it, then
        // check the interpolated nutrient concentration is between the min
        // and max of the nutrient concentrations on the nodes of that element
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
            cell_iter != tissue.End();
            ++cell_iter)
        {
            unsigned elem_index = simulator.mpCoarseNutrientMesh->GetContainingElementIndex(tissue.GetLocationOfCellCentre(*cell_iter));
            Element<2,2>* p_element = simulator.mpCoarseNutrientMesh->GetElement(elem_index);


            double max = std::max(nutrient_conc[p_element->GetNodeGlobalIndex(0)],
                                  nutrient_conc[p_element->GetNodeGlobalIndex(1)]);

            max = std::max(max, nutrient_conc[p_element->GetNodeGlobalIndex(2)]);

            double min = std::min(nutrient_conc[p_element->GetNodeGlobalIndex(0)],
                                  nutrient_conc[p_element->GetNodeGlobalIndex(1)]);

            min = std::min(min, nutrient_conc[p_element->GetNodeGlobalIndex(2)]);


            double value_at_cell = CellwiseData<2>::Instance()->GetValue(*cell_iter, 0);

            TS_ASSERT_LESS_THAN_EQUALS(min, value_at_cell);
            TS_ASSERT_LESS_THAN_EQUALS(value_at_cell, max);
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    void TestCoarseNutrientMeshBoundaryConditionImplementation() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        TissueConfig::Instance()->SetHepaOneParameters();

        // Create a cigar-shaped mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        MutableMesh<2,2>* p_mesh = new MutableMesh<2,2>;
        p_mesh->ConstructFromMeshReader(mesh_reader);
        p_mesh->Scale(5.0,1.0);

        // Set up cells
        std::vector<TissueCell> cells;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel(2));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                    (TissueConfig::Instance()->GetHepaOneCellG1Duration()
                                    +TissueConfig::Instance()->GetSG2MDuration());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Set up tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Set up CellwiseData and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(), 1);
        p_data->SetTissue(tissue);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i));
        }

        // Set up PDE
        AveragedSinksPde<2> pde(tissue, -0.01);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation to use a coarse nutrient mesh
        TissueSimulationWithNutrients<2> simulator(tissue, force_collection, NULL, &pde);
        simulator.SetOutputDirectory("TestCoarseNutrientMeshBoundaryConditionImplementation");
        simulator.SetEndTime(0.01);
        simulator.UseCoarseNutrientMesh(2.0);

        // Run tissue simulation
        simulator.Solve();

        // Test that boundary cells experience the right boundary condition
        for (AbstractTissue<2>::Iterator cell_iter = simulator.rGetTissue().Begin();
             cell_iter != simulator.rGetTissue().End();
             ++cell_iter)
        {
            if ( (static_cast<AbstractCellCentreBasedTissue<2>*>(&(simulator.rGetTissue())))->GetNodeCorrespondingToCell(*cell_iter)->IsBoundaryNode() )
            {
                TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), 1.0, 1e-1);
            }
        }

        // Tidy up
        CellwiseData<2>::Destroy();
        delete p_mesh;
    }


    void TestArchivingWithSimplePde() throw (Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        TissueConfig::Instance()->SetHepaOneParameters();

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel(2));
            double birth_time = -1.0 - ( (double) i/p_mesh->GetNumNodes() )*
                                            (TissueConfig::Instance()->GetHepaOneCellG1Duration()
                                             +TissueConfig::Instance()->GetSG2MDuration());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Set up tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Set up CellwiseData and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(),1);
        p_data->SetTissue(tissue);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i));
        }

        // Set up PDE
        SimpleNutrientPde<2> pde(0.1);

        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, force_collection, &pde);
        simulator.SetOutputDirectory("TissueSimulationWithNutrientsSaveAndLoad");
        simulator.SetEndTime(0.2);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&tissue);
        simulator.AddCellKiller(&killer);

        // Run tissue simulation
        simulator.Solve();

        // Save tissue simulation
        TissueSimulationArchiver<2, TissueSimulationWithNutrients<2> >::Save(&simulator);

        TissueSimulationWithNutrients<2>* p_simulator
            = TissueSimulationArchiver<2, TissueSimulationWithNutrients<2> >::Load("TissueSimulationWithNutrientsSaveAndLoad", 0.2);
        p_simulator->SetPde(&pde);

        p_simulator->SetEndTime(0.5);
        p_simulator->Solve();

        // These results are from time 0.5 in TestWithOxygen.
        std::vector<double> node_5_location = p_simulator->GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.4968, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 0.8635, 1e-4);

        std::vector<double> node_15_location = p_simulator->GetNodeLocation(15);
        TS_ASSERT_DELTA(node_15_location[0], 0.4976, 1e-4);
        TS_ASSERT_DELTA(node_15_location[1], 2.5977, 1e-4);

        // Test CellwiseData was set up correctly
        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);

        // Test the CellwiseData result
        TS_ASSERT_DELTA(p_data->GetValue(p_simulator->rGetTissue().rGetCellUsingLocationIndex(5)), 0.9604, 1e-4);
        TS_ASSERT_DELTA(p_data->GetValue(p_simulator->rGetTissue().rGetCellUsingLocationIndex(15)), 0.9584, 1e-4);

        // Run tissue simulation
        delete p_simulator;
        CellwiseData<2>::Destroy();
    }


    /**
     * This test demonstrates how to archive a TissueSimulationWithNutrients
     * in the case where the nutrient PDE has the tissue as a member variable
     */
    void TestArchivingWithCellwisePde() throw (Exception)
    {
        if (!PetscTools::IsSequential())
        {
            TS_TRACE("This test does not pass in parallel yet.");
            return;
        }

        TissueConfig::Instance()->SetHepaOneParameters();

        std::string output_directory = "TestArchivingWithCellwisePde";
        double end_time = 0.1;

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel(2));
            double birth_time = -1.0 - ( (double) i/p_mesh->GetNumNodes() )*
                                    (TissueConfig::Instance()->GetHepaOneCellG1Duration()
                                    +TissueConfig::Instance()->GetSG2MDuration());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Set up tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Set up CellwiseData and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(), 1);
        p_data->SetTissue(tissue);

        // Set initial conditions for CellwiseData (needed to avoid memory errors)
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i));
        }

        // Set up PDE
        CellwiseNutrientSinkPde<2> pde(tissue, 0.03);

        // Set up mechanics system
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(3.0);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, force_collection, &pde);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(end_time);

        // Run tissue simulation
        simulator.Solve();

        // Save tissue simulation
        TissueSimulationArchiver<2, TissueSimulationWithNutrients<2> >::Save(&simulator);

        // Load simulation
        TissueSimulationWithNutrients<2>* p_simulator
            = TissueSimulationArchiver<2, TissueSimulationWithNutrients<2> >::Load(output_directory, end_time);

        /**
         * In this case, the PDE had a reference to the tissue. To avoid a
         * segmentation fault, we need to first get the archived tissue, pass
         * it in to a new instance of our PDE, taking care to use the same
         * consumption rate as before. We then pass this PDE into the tissue
         * simulation.
         */
        MeshBasedTissue<2>* p_tissue = static_cast<MeshBasedTissue<2>*>(&(p_simulator->rGetTissue()));
        CellwiseNutrientSinkPde<2> pde2(*p_tissue, 0.03);
        p_simulator->SetPde(&pde2);
        p_simulator->SetEndTime(2.0*end_time);

        // Run tissue simulation
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

        // Tidy up
        delete p_simulator;
        CellwiseData<2>::Destroy();
    }


    void Test3DTissueSimulationWithNutrients() throw(Exception)
    {
        EXIT_IF_PARALLEL; //defined in PetscTools

        TissueConfig::Instance()->SetHepaOneParameters();

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("TestSolveMethodSpheroidSimulation3DMesh","StartMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<3> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Set up tissue
        MeshBasedTissue<3> tissue(mesh, cells);

        // Set up CellwiseData and associate it with the tissue
        CellwiseData<3>* p_data = CellwiseData<3>::Instance();
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(),1);
        p_data->SetTissue(tissue);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i));
        }

        // Set up PDE
        SimpleNutrientPde<3> pde(0.1);

        // Set up force law
        GeneralisedLinearSpringForce<3> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<3>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulationWithNutrients<3> simulator(tissue, force_collection, &pde);
        simulator.SetOutputDirectory("TissueSimulationWithOxygen3d");
        simulator.SetEndTime(0.5);

        // Set up cell killer and pass into simulation
        AbstractCellKiller<3>* p_killer = new OxygenBasedCellKiller<3>(&tissue);
        simulator.AddCellKiller(p_killer);

        // Coverage
        TS_ASSERT_THROWS_THIS(simulator.CreateCoarseNutrientMesh(10.0), "This method is only implemented in 2D");

        // Run tissue simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Tidy up
        CellwiseData<3>::Destroy();

        delete p_killer;
    }

};
#endif /*TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_*/
