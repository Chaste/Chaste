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

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

#include "SimulationTime.hpp"
#include "RandomNumberGenerator.hpp"
#include "ExecutableSupport.hpp"

/*
 * These headers handle passing parameters to the executable.
 */
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <boost/lexical_cast.hpp>

#include "OffLatticeSimulation.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"

/*
 * Prototype functions
 */
void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(unsigned simulationId, unsigned numMeshPoints);
void OutputOnCompletion(unsigned simulationId, unsigned numMeshPoints);

int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StandardStartup(&argc, &argv);

    // Define command line options
    boost::program_options::options_description general_options("This is a Chaste Immersed Boundary executable.\n");
    general_options.add_options()
                    ("help", "produce help message")
                    ("ID", boost::program_options::value<unsigned>()->default_value(0),"Index of the simulation")
                    ("MP", boost::program_options::value<unsigned>()->default_value(32),"Number of fluid mesh points");

    // Parse command line into variables_map
    boost::program_options::variables_map variables_map;
    boost::program_options::store(parse_command_line(argc, argv, general_options), variables_map);

    // Print help message if wanted
    if (variables_map.count("help"))
    {
        std::cout << setprecision(3) << general_options << "\n";
        std::cout << general_options << "\n";
        return 1;
    }

    // Get ID and name from command line
    unsigned simulation_id = variables_map["ID"].as<unsigned>();
    unsigned num_mesh_points = variables_map["MP"].as<unsigned>();

    SetupSingletons();
    SetupAndRunSimulation(simulation_id, num_mesh_points);
    DestroySingletons();
    OutputOnCompletion(simulation_id, num_mesh_points);
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
    // This is from the tearDown method of the test suite
    SimulationTime::Destroy();
    RandomNumberGenerator::Destroy();
    CellPropertyRegistry::Instance()->Clear();
}

void OutputOnCompletion(unsigned simulationId, unsigned numMeshPoints)
{
    // Compose the message
    std::stringstream message;
    message << "Completed simulation with ID " << simulationId << " and mesh spacing " << 1.0 / numMeshPoints << std::endl;

    // Send it to the console
    std::cout << message.str() << std::flush;
}
void SetupAndRunSimulation(unsigned simulationId, unsigned numMeshPoints)
{
    /*
     * 1: Num cells
     * 2: Num nodes per cell
     * 3: Superellipse exponent
     * 4: Superellipse aspect ratio
     * 5: Random y-variation
     * 6: Include membrane
     */
    ImmersedBoundaryPalisadeMeshGenerator gen(1, 8192, 1.0, 2.0, 0.0, false);
    ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

    p_mesh->SetNumGridPtsXAndY(numMeshPoints);

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
    p_boundary_force->SetSpringConstant(1e7);
    p_boundary_force->SetRestLengthMultiplier(0.5);

    // Create and set an output directory that is different for each simulation
    std::stringstream output_directory;
    output_directory << "convergence/mesh/sim/" << simulationId;
    simulator.SetOutputDirectory(output_directory.str());

    double end_time = 10.0;
    double dt = 0.01;

    // Set simulation properties and solve
    simulator.SetDt(dt);
    simulator.SetSamplingTimestepMultiple(10);
    simulator.SetEndTime(end_time);
    simulator.Solve();

    // Prepare results file to record simulation statistics
    OutputFileHandler results_handler(output_directory.str(), false);
    out_stream results_file = results_handler.OpenOutputFile("results.dat");

    double esf_at_end = p_mesh->GetElongationShapeFactorOfElement(0);

    // Output summary statistics to results file.  lexical_cast is a convenient way to output doubles at max precision.
    (*results_file) << simulationId << ","
                    << boost::lexical_cast<std::string>(1.0 / numMeshPoints) << ","
                    << boost::lexical_cast<std::string>(esf_at_end);

    // Tidy up
    results_file->close();
}
