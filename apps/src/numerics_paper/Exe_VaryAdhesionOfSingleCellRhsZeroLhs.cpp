/*

Copyright (c) 2005-2015, University of Oxford.
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

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

#include "SimulationTime.hpp"
#include "RandomNumberGenerator.hpp"
#include "ExecutableSupport.hpp"

#include "Debug.hpp"

/*
 * These headers handle passing parameters to the executable.
 */
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include "OffLatticeSimulation.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"
#include "ImmersedBoundaryCellCellInteractionForce.hpp"

/*
 * Prototype functions
 */
void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(unsigned kick, unsigned localSpringConst, unsigned globalSpringConst);
void OutputOnCompletion(unsigned kick, unsigned localSpringConst, unsigned globalSpringConst);

int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StandardStartup(&argc, &argv);

    // Define command line options
    boost::program_options::options_description general_options("This is a Chaste Immersed Boundary executable.\n");
    general_options.add_options()
                    ("help", "produce help message")
                    ("K", boost::program_options::value<unsigned>()->default_value(0),"Amount of kick for the simulation")
                    ("L", boost::program_options::value<unsigned>()->default_value(0),"Local multiplier for the cell-cell spring const")
                    ("G", boost::program_options::value<unsigned>()->default_value(0),"Global multiplier for the cell-cell spring const");

    // define parse command line into variables_map
    boost::program_options::variables_map variables_map;
    boost::program_options::store(parse_command_line(argc, argv, general_options), variables_map);

    // print help message if wanted
    if (variables_map.count("help"))
    {
        std::cout << setprecision(3) << general_options << "\n";
        std::cout << general_options << "\n";
        return 1;
    }

    // get id and name from command line
    unsigned kick = variables_map["K"].as<unsigned>();
    unsigned local = variables_map["L"].as<unsigned>();
    unsigned global = variables_map["G"].as<unsigned>();

    SetupSingletons();
    SetupAndRunSimulation(kick, local, global);
    DestroySingletons();
    OutputOnCompletion(kick, local, global);
}

void SetupSingletons()
{
    // Set up what the test suite would do
    SimulationTime::Instance()->SetStartTime(0.0);

    // Reseed with 0 for same random numbers each time, or time(NULL) or simulation_id to change each realisation
    RandomNumberGenerator::Instance()->Reseed(0);
    CellPropertyRegistry::Instance()->Clear();
    CellId::ResetMaxCellId();
}

void DestroySingletons()
{
    // this is from the tearDown method of the test suite
    SimulationTime::Destroy();
    RandomNumberGenerator::Destroy();
    CellPropertyRegistry::Instance()->Clear();
}

void OutputOnCompletion(unsigned kick, unsigned localSpringConst, unsigned globalSpringConst)
{
    // Compose the message
    std::stringstream message;
    message << "Completed simulation with global SC " << globalSpringConst << ", local SC " << localSpringConst << ", and kick " << kick << std::endl;

    // Send it to the console
    std::cout << message.str() << std::flush;
}
void SetupAndRunSimulation(unsigned kick, unsigned localSpringConst, unsigned globalSpringConst)
{
    double reference_spring_const = 5.0 * 1e6;                      // the reference global spring const

    double fp_global_sc = 0.5 + 0.25 * (double)(globalSpringConst); // global cell-cell spring constant multiplier
    double fp_local_sc = 0.2 * (double)localSpringConst;            // local cell-cell spring constant multiplier
    double fp_sim_kick = 1.0 + 0.075 * (double)kick;                 // the amount we kick the appropriate cell

    /*
     * 1: Num cells
     * 2: Num nodes per cell
     * 3: Superellipse exponent
     * 4: Superellipse aspect ratio
     * 5: Random y-variation
     * 6: Include membrane
     */
    ImmersedBoundaryPalisadeMeshGenerator gen(7, 128, 0.1, 2.5, 0.0, true);
    ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

    p_mesh->SetNumGridPtsXAndY(256);

    std::vector<CellPtr> cells;
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
    CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
    cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

    ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
    cell_population.SetIfPopulationHasActiveSources(false);

    OffLatticeSimulation<2> simulator(cell_population);

    // Add main immersed boundary simulation modifier
    MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
    simulator.AddSimulationModifier(p_main_modifier);

    // Add force laws
    MAKE_PTR_ARGS(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force, (cell_population));
    p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
    p_boundary_force->SetSpringConstant(reference_spring_const);

    // Create and set an output directory that is different for each simulation
    std::stringstream output_directory;
    output_directory << "numerics_paper/Exe_VaryAdhesionOfSingleCellRhsZeroLhs/sim/"
                     << globalSpringConst << "_" << localSpringConst << "_" << kick;
    simulator.SetOutputDirectory(output_directory.str());

    // Set simulation properties
    double dt = 0.075;
    simulator.SetDt(dt);
    simulator.SetSamplingTimestepMultiple(100);
    simulator.SetEndTime(100.0 * dt);
    simulator.Solve();

    // Now we have relaxed the mesh, we reset the start time to zero.  This overwrites the output, so we end up with
    // less clutter
    SimulationTime::Instance()->Destroy();
    SimulationTime::Instance()->SetStartTime(0.0);

    // Add a cell-cell interaction force with the same intrinsic strength as the membrane force
    MAKE_PTR_ARGS(ImmersedBoundaryCellCellInteractionForce<2>, p_cell_cell_force, (cell_population));
    p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
    p_cell_cell_force->SetSpringConstant(reference_spring_const * fp_global_sc);

    // Get height of basement lamina
    double lamina_height = 0.0;
    for (unsigned node_idx = 0 ; node_idx < p_mesh->GetElement(0)->GetNumNodes() ; node_idx++)
    {
        lamina_height += p_mesh->GetElement(0)->GetNode(node_idx)->rGetModifiableLocation()[1];
    }
    lamina_height /= p_mesh->GetElement(0)->GetNumNodes();

    // Kick the second cell in from the left and set its E-cad level
    unsigned e_cad_location = p_cell_cell_force->rGetProteinNodeAttributeLocations()[0];

    double x_centroid_before = p_mesh->GetCentroidOfElement(3)[0];

    for (unsigned node_idx = 0 ; node_idx < p_mesh->GetElement(3)->GetNumNodes() ; node_idx++)
    {
        // Scale the node to its new height based on the kick parameter
        double new_height = lamina_height + fp_sim_kick * (p_mesh->GetElement(3)->GetNode(node_idx)->rGetLocation()[1] - lamina_height);
        p_mesh->GetElement(3)->GetNode(node_idx)->rGetModifiableLocation()[1] = new_height;

        // Scale the e_cad level so as to alter the adhesion just of the RHS specific cell
        if (p_mesh->GetElement(3)->GetNode(node_idx)->rGetLocation()[0] > x_centroid_before)
        {
            p_mesh->GetElement(3)->GetNode(node_idx)->rGetNodeAttributes()[e_cad_location] = fp_local_sc;
        }
        else
        {
            p_mesh->GetElement(3)->GetNode(node_idx)->rGetNodeAttributes()[e_cad_location] = 0.0;
        }
    }

    ChasteCuboid<2> bounding_box_before = p_mesh->CalculateBoundingBoxOfElement(3);

    simulator.SetSamplingTimestepMultiple(160);
    simulator.SetEndTime(16000.0 * dt);
    simulator.Solve();

    double x_centroid_after = p_mesh->GetCentroidOfElement(3)[0];
    ChasteCuboid<2> bounding_box_after = p_mesh->CalculateBoundingBoxOfElement(3);

    OutputFileHandler results_handler(output_directory.str(), false);
    out_stream results_file = results_handler.OpenOutputFile("results.dat");

    // Calculate summary statistics

    double width = bounding_box_after.GetWidth(0);

    double ss_centroid = x_centroid_after - x_centroid_before;
    double ss_symmetry = (bounding_box_after.rGetUpperCorner()[0] - x_centroid_before) / width;

    // Output summary statistics to results file
    (*results_file) << fp_global_sc << "," << fp_local_sc << "," << fp_sim_kick << ","
                    << ss_centroid << "," << ss_symmetry;

    // Tidy up
    results_file->close();
}
