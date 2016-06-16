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

// Needed for the test environment
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "OffLatticeSimulation.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"
#include "OutputFileHandler.hpp"

#include <boost/lexical_cast.hpp>

#include "Debug.hpp"
#include "Timer.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestNumericsPaperSimulations : public AbstractCellBasedTestSuite
{
public:

    void xTestEllipseRelaxing() throw(Exception)
    {
        /*
         * 1: num nodes
         * 2: superellipse exponent
         * 3: cell width
         * 4: cell height
         * 5: bottom left x
         * 6: bottom left y
         */
        SuperellipseGenerator* p_gen = new SuperellipseGenerator(128, 1.0, 0.4, 0.6, 0.3, 0.2);
        std::vector<c_vector<double, 2> > locations = p_gen->GetPointsAsVectors();

        std::vector<Node<2>* > nodes;
        std::vector<ImmersedBoundaryElement<2,2>* > elements;

        for (unsigned location = 0; location < locations.size(); location++)
        {
            nodes.push_back(new Node<2>(location, locations[location], true));
        }

        elements.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));
        elements[0]->rGetCornerNodes().push_back(nodes[0]);
        elements[0]->rGetCornerNodes().push_back(nodes[1]);
        elements[0]->rGetCornerNodes().push_back(nodes[2]);
        elements[0]->rGetCornerNodes().push_back(nodes[3]);

        ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2,2>(nodes, elements);
        p_mesh->SetNumGridPtsXAndY(32);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIfPopulationHasActiveSources(false);

        OffLatticeSimulation<2> simulator(cell_population);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force laws
        MAKE_PTR(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force);
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetSpringConstant(1e8);

        std::string output_directory = "numerics_paper/ellipse_relaxing";
        simulator.SetOutputDirectory(output_directory);

        // Write the state of the immersed boundary mesh to file
        ImmersedBoundaryMeshWriter<2,2> mesh_at_start(output_directory, "example_simulation_mesh_at_start");
        ImmersedBoundaryMeshWriter<2,2> mesh_at_end(output_directory, "example_simulation_mesh_at_end");

        OutputFileHandler results_handler(output_directory, false);
        out_stream results_file = results_handler.OpenOutputFile("example_simulation_esf.dat");

        // Output summary statistics to results file
        (*results_file) << "time,esf\n";
        (*results_file) << 0.0 << "," << p_mesh->GetElongationShapeFactorOfElement(0) << "\n";

        // Set simulation properties
        double dt = 0.05;
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(1);

        for (unsigned i=0; i< 100; i++)
        {
            double new_end_time = dt * (1.0 + i);

            simulator.SetEndTime(new_end_time);
            simulator.Solve();

            (*results_file) << new_end_time << "," << p_mesh->GetElongationShapeFactorOfElement(0) << "\n";

            if (i==0)
            {
                mesh_at_start.WriteFilesUsingMesh(*p_mesh);
            }
        }

        mesh_at_end.WriteFilesUsingMesh(*p_mesh);

        // Tidy up
        results_file->close();
        delete(p_mesh);
    }

    void xTestSingleCellVolumeChangeWithNodeSpacing() throw(Exception)
    {
        /**
         * This test simulates a single circular cell for a fixed simulation time.
         *
         * We test node spacing ratios 0.1, 0.2, ..., 3.9, 4.0
         *
         * All parameters are fixed except the number of nodes, and the following are exported to a csv file:
         *  * Number of nodes in the cell
         *  * Node spacing to mesh spacing ratio
         *  * Change in volume as ratio: absolute change / initial volume
         */

        unsigned num_sims = 40;
        unsigned num_grid_pts = 256;
        double radius = 0.4;

        std::string output_directory = "numerics_paper/node_spacing_ratio";

        OutputFileHandler results_handler(output_directory, false);
        out_stream results_file = results_handler.OpenOutputFile("numerical_results_nsr.dat");

        // Output summary statistics to results file
        (*results_file) << "id,node_spacing_ratio,absolute_volume_change\n";

        for (unsigned sim_idx = 0; sim_idx < num_sims; sim_idx++)
        {
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            double target_spacing_ratio = 0.1 * (1 + sim_idx);
            unsigned num_nodes = (unsigned)(M_PI / asin(0.5 * target_spacing_ratio / (radius * num_grid_pts)));

            /*
             * Create an immersed boundary mesh using a SuperellipseGenerator
             *
             * 1: Num nodes
             * 2: Superellipse exponent
             * 3: Width
             * 4: Height
             * 5: Bottom-left x
             * 6: Botton-left y
             */
            SuperellipseGenerator *p_gen = new SuperellipseGenerator(num_nodes, 1.0, 2.0*radius, 2.0*radius, 0.5-radius, 0.5-radius);

            // Generate a mesh using this superellipse
            std::vector<c_vector<double, 2> > locations = p_gen->GetPointsAsVectors();
            delete p_gen;

            std::vector<Node<2> *> nodes;
            for (unsigned node_idx = 0; node_idx < locations.size(); node_idx++)
            {
                nodes.push_back(new Node<2>(node_idx, locations[node_idx], true));
            }

            std::vector<ImmersedBoundaryElement < 2, 2>* > elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements);
            mesh.SetNumGridPtsXAndY(256);

            // The real node spacing ratio (it won't be exactly the target due to needing an integer number of nodes)
            double node_spacing_ratio = mesh.GetSpacingRatio();

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            OffLatticeSimulation<2> simulation(cell_population);
            cell_population.SetIfPopulationHasActiveSources(false);

            // Add main immersed boundary simulation modifier
            MAKE_PTR(ImmersedBoundarySimulationModifier < 2 >, p_main_modifier);
            simulation.AddSimulationModifier(p_main_modifier);

            // Add force law
            MAKE_PTR(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force);
            p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
            p_boundary_force->SetSpringConstant(1e9);
            p_boundary_force->SetRestLengthMultiplier(0.5);

            std::string sim_output_dir = output_directory + "/sim";
            sim_output_dir += boost::lexical_cast<std::string>(sim_idx);

            // Set simulation properties
            double dt = 0.01;
            simulation.SetOutputDirectory(sim_output_dir);
            simulation.SetDt(dt);
            simulation.SetSamplingTimestepMultiple(5);
            simulation.SetEndTime(1000.0 * dt);

            // Run the simulation
            simulation.Solve();

            double new_radius = mesh.GetElement(0)->GetNode(0)->rGetLocation()[0] - 0.5;
            double absolute_volume_change = fabs(radius * radius - new_radius * new_radius) / (radius * radius);

            (*results_file) << boost::lexical_cast<std::string>(sim_idx) << ","
                            << boost::lexical_cast<std::string>(node_spacing_ratio) << ","
                            << boost::lexical_cast<std::string>(absolute_volume_change) << "\n";
        }

        results_file->close();
    }


    void TestIntraCellularParameterScaling() throw(Exception)
    {
        /**
         * This test runs two simulations of the same elliptical membrane twice with different
         * numbers of nodes, in order to demonstrate the inbuilt scaling correctly calculates
         * the necessary intra-cellular spring constant to compensate.
         *
         * All parameters are fixed except the number of nodes, and the following are exported to a csv file:
         *  * The time points at which the ESF is sampled
         *  * The ESF for the first scenario, at each time point
         *  * The ESF for the second scenario, at each time point
         */

        std::string output_directory = "numerics_paper/intra_cellular_scaling";

        OutputFileHandler results_handler(output_directory, false);
        out_stream results_file = results_handler.OpenOutputFile("numerical_results_intra_scaling.dat");

        // Output summary statistics to results file
        (*results_file) << "time,esf_256,esf_512\n";

        // Vectors to store summary statistics
        std::vector<double> time_points;
        std::vector<double> esf_256;
        std::vector<double> esf_512;

        // Initial time value at start of simulation
        time_points.push_back(0.0);

        unsigned num_time_pts = 20;
        unsigned num_steps_per_sample = 50;

        // Sim with 256 nodes
        {
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            /*
             * 1: Num cells
             * 2: Num nodes per cell
             * 3: Superellipse exponent
             * 4: Superellipse aspect ratio
             * 5: Random y-variation
             * 6: Include membrane
             */
            ImmersedBoundaryPalisadeMeshGenerator gen(1, 256, 1.0, 2.0, 0.0, false);
            ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

            p_mesh->SetNumGridPtsXAndY(256);

            // Initial value at start of simulation
            esf_256.push_back(p_mesh->GetElongationShapeFactorOfElement(0));

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

            OffLatticeSimulation<2> simulation(cell_population);
            cell_population.SetIfPopulationHasActiveSources(false);

            // Add main immersed boundary simulation modifier
            MAKE_PTR(ImmersedBoundarySimulationModifier < 2 >, p_main_modifier);
            simulation.AddSimulationModifier(p_main_modifier);

            // Add force law
            MAKE_PTR(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force);
            p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
            p_boundary_force->SetSpringConstant(1e7);
            p_boundary_force->SetRestLengthMultiplier(0.5);

            // Simulation output directory
            std::string output_dir_256 = output_directory + '/' + boost::lexical_cast<std::string>(256);

            // Set simulation properties
            double dt = 0.05;
            simulation.SetDt(dt);
            simulation.SetSamplingTimestepMultiple(1);
            simulation.SetOutputDirectory(output_dir_256);

            for (unsigned i=0; i<num_time_pts; i++)
            {
                double new_end_time = num_steps_per_sample * dt * (1.0 + i);

                simulation.SetEndTime(new_end_time);
                simulation.Solve();

                time_points.push_back(new_end_time);
                esf_256.push_back(p_mesh->GetElongationShapeFactorOfElement(0));
            }
        }

        // Sim with 512 nodes
        {
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            /*
             * 1: Num cells
             * 2: Num nodes per cell
             * 3: Superellipse exponent
             * 4: Superellipse aspect ratio
             * 5: Random y-variation
             * 6: Include membrane
             */
            ImmersedBoundaryPalisadeMeshGenerator gen(1, 512, 1.0, 2.0, 0.0, false);
            ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

            p_mesh->SetNumGridPtsXAndY(256);

            // Initial value at start of simulation
            esf_512.push_back(p_mesh->GetElongationShapeFactorOfElement(0));

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

            OffLatticeSimulation<2> simulation(cell_population);
            cell_population.SetIfPopulationHasActiveSources(false);

            // Add main immersed boundary simulation modifier
            MAKE_PTR(ImmersedBoundarySimulationModifier < 2 >, p_main_modifier);
            simulation.AddSimulationModifier(p_main_modifier);

            // Add force law
            MAKE_PTR(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force);
            p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
            p_boundary_force->SetSpringConstant(1e7);
            p_boundary_force->SetRestLengthMultiplier(0.5);

            // Simulation output directory
            std::string output_dir_512 = output_directory + '/' + boost::lexical_cast<std::string>(512);

            // Set simulation properties
            double dt = 0.05;
            simulation.SetDt(dt);
            simulation.SetSamplingTimestepMultiple(1);
            simulation.SetOutputDirectory(output_dir_512);

            for (unsigned i=0; i<num_time_pts; i++)
            {
                double new_end_time = num_steps_per_sample * dt * (1.0 + i);

                simulation.SetEndTime(new_end_time);
                simulation.Solve();

                esf_512.push_back(p_mesh->GetElongationShapeFactorOfElement(0));
            }
        }

        // Check we output arrays all contain the same number of elements
        TS_ASSERT_EQUALS(time_points.size(), esf_256.size());
        TS_ASSERT_EQUALS(time_points.size(), esf_512.size());

        for(unsigned time_pt = 0 ; time_pt < time_points.size() ; time_pt++)
        {
            (*results_file) << boost::lexical_cast<std::string>(time_points[time_pt]) << ","
                            << boost::lexical_cast<std::string>(esf_256[time_pt]) << ","
                            << boost::lexical_cast<std::string>(esf_512[time_pt]) << "\n";
        }

        results_file->close();
    }
};
