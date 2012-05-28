/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTCELLBASEDPDEHANDLER_HPP_
#define TESTCELLBASEDPDEHANDLER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <ctime>

#include "ArchiveOpener.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "CellBasedPdeHandler.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellwiseSourcePde.hpp"
#include "SimpleUniformSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AveragedSourcePde.hpp"
#include "NumericFileComparison.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"



class SimplePdeForTesting : public AbstractLinearEllipticPde<2,2>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>&, Element<2,2>* pElement)
    {
        return -1.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>&, Element<2,2>*)
    {
        return 0.0;
    }

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }
};

class TestCellBasedPdeHandler : public AbstractCellBasedTestSuite
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

    void TestConstructor() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a cell population
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Test that member variables are initialised correctly
        TS_ASSERT_EQUALS(pde_handler.GetCellPopulation(), &cell_population);
        TS_ASSERT_EQUALS(pde_handler.GetWriteAverageRadialPdeSolution(), false);
        TS_ASSERT_EQUALS(pde_handler.GetWriteDailyAverageRadialPdeSolution(), false);
        TS_ASSERT_EQUALS(pde_handler.GetImposeBcsOnCoarseBoundary(), true);
        TS_ASSERT_EQUALS(pde_handler.GetNumRadialIntervals(), UNSIGNED_UNSET);
        TS_ASSERT(pde_handler.GetCoarsePdeMesh() == NULL);
    }

    void TestSetMethods() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a cell population
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Test set and get methods
        pde_handler.SetWriteAverageRadialPdeSolution("averaged quantity");

        TS_ASSERT_EQUALS(pde_handler.GetWriteAverageRadialPdeSolution(), true);
        TS_ASSERT_EQUALS(pde_handler.GetWriteDailyAverageRadialPdeSolution(), false);
        TS_ASSERT_EQUALS(pde_handler.GetNumRadialIntervals(), 10u);

        pde_handler.SetWriteAverageRadialPdeSolution("averaged quantity", 5, true);

        TS_ASSERT_EQUALS(pde_handler.GetWriteAverageRadialPdeSolution(), true);
        TS_ASSERT_EQUALS(pde_handler.GetWriteDailyAverageRadialPdeSolution(), true);
        TS_ASSERT_EQUALS(pde_handler.GetNumRadialIntervals(), 5u);

        pde_handler.SetImposeBcsOnCoarseBoundary(false);
        TS_ASSERT_EQUALS(pde_handler.GetImposeBcsOnCoarseBoundary(), false);

        // Test AddPdeAndBc()
        SimpleUniformSourcePde<2> pde(-0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("averaged quantity");

        unsigned num_nodes = mesh.GetNumNodes();
        std::vector<double> data(num_nodes);
        for (unsigned i=0; i<num_nodes; i++)
        {
            data[i] = i + 0.45;
        }

        Vec vector = PetscTools::CreateVec(data);
        pde_and_bc.SetSolution(vector);

        pde_handler.AddPdeAndBc(&pde_and_bc);

        TS_ASSERT_EQUALS(pde_handler.mPdeAndBcCollection[0]->IsNeumannBoundaryCondition(), false);
        TS_ASSERT_EQUALS(pde_handler.mPdeAndBcCollection[0]->HasAveragedSourcePde(), false);

        ReplicatableVector solution(pde_handler.GetPdeSolution());

        TS_ASSERT_EQUALS(solution.GetSize(), num_nodes);
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(solution[i], i + 0.45, 1e-4);
        }
    }

    void TestCellBasedPdeHandlerOutputParameters() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        std::string output_directory = "TestCellBasedPdeHandlerOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Create a cell population
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Test output methods
        TS_ASSERT_EQUALS(pde_handler.GetIdentifier(), "CellBasedPdeHandler-2");

        out_stream pde_handler_parameter_file = output_file_handler.OpenOutputFile("CellBasedPdeHandler.parameters");
        pde_handler.OutputParameters(pde_handler_parameter_file);
        pde_handler_parameter_file->close();

        std::string pde_handler_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + pde_handler_results_dir + "CellBasedPdeHandler.parameters cell_based/test/data/TestCellBasedPdeHandler/CellBasedPdeHandler.parameters").c_str()), 0);
    }

    void TestArchiving() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "CellBasedPdeHandler.arch";
        ArchiveLocationInfo::SetMeshFilename("pde_handler_mesh");

        {
            // Create a cell population
            HoneycombMeshGenerator generator(2, 2, 0);
            MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

            NodeBasedCellPopulation<2> cell_population(mesh, cells);

            // Create a PDE handler object using this cell population
            CellBasedPdeHandler<2>* const p_pde_handler = new CellBasedPdeHandler<2>(&cell_population);

            // Set member variables for testing
            p_pde_handler->SetWriteAverageRadialPdeSolution("averaged quantity", 5, true);
            p_pde_handler->SetImposeBcsOnCoarseBoundary(false);

            // Set up PDE and pass to handler
            AveragedSourcePde<2> pde(cell_population, -0.1);
            ConstBoundaryCondition<2> bc(1.0);
            PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
            pde_and_bc.SetDependentVariableName("averaged quantity");
            p_pde_handler->AddPdeAndBc(&pde_and_bc);

            // Test UseCoarsePdeMesh() again
            ChastePoint<2> lower(0.0, 0.0);
            ChastePoint<2> upper(9.0, 9.0);
            ChasteCuboid<2> cuboid(lower, upper);
            p_pde_handler->UseCoarsePdeMesh(3.0, cuboid, true);

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Archive object
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 11);
            (*p_arch) << static_cast<const SimulationTime&>(*p_simulation_time);
            (*p_arch) << p_pde_handler;

            // Tidy up
            SimulationTime::Destroy();
            delete p_pde_handler;
        }

        {
            CellBasedPdeHandler<2>* p_pde_handler;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore object from the archive
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            (*p_arch) >> *p_simulation_time;
            (*p_arch) >> p_pde_handler;

            // Test that the member variables were archived correctly
            TS_ASSERT_EQUALS(p_pde_handler->mpCellPopulation->GetNumRealCells(), 4u);
            TS_ASSERT_EQUALS(p_pde_handler->GetWriteAverageRadialPdeSolution(), true);
            TS_ASSERT_EQUALS(p_pde_handler->GetWriteDailyAverageRadialPdeSolution(), true);
            TS_ASSERT_EQUALS(p_pde_handler->GetImposeBcsOnCoarseBoundary(), false);
            TS_ASSERT_EQUALS(p_pde_handler->GetNumRadialIntervals(), 5u);
            TS_ASSERT_EQUALS(p_pde_handler->mAverageRadialSolutionVariableName, "averaged quantity");

            ///\todo we currently do not archive mpCoarsePdeMesh - consider doing this (#1891)
            TS_ASSERT(p_pde_handler->GetCoarsePdeMesh() == NULL);

            TS_ASSERT_EQUALS(p_pde_handler->mPdeAndBcCollection.size(), 1u);
            TS_ASSERT_EQUALS(p_pde_handler->mPdeAndBcCollection[0]->IsNeumannBoundaryCondition(), false);

            // Tidy up
            delete p_pde_handler->mpCellPopulation;
            delete p_pde_handler;
        }
    }

    void TestUseCoarsePdeMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a cell population
        HoneycombMeshGenerator generator(4, 4, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Test that UseCoarsePdeMesh() throws an exception if no PDEs are specified
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(9.0, 9.0);
        ChasteCuboid<2> cuboid(lower, upper);
        TS_ASSERT_THROWS_THIS(pde_handler.UseCoarsePdeMesh(3.0, cuboid, true),
            "mPdeAndBcCollection should be populated prior to calling UseCoarsePdeMesh().");

        // Set up PDE and pass to handler
        AveragedSourcePde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Test UseCoarsePdeMesh() again
        pde_handler.UseCoarsePdeMesh(3.0, cuboid, true);

        // Test that the coarse mesh has the correct number of nodes and elements
        TetrahedralMesh<2,2>* p_coarse_mesh = pde_handler.GetCoarsePdeMesh();
        TS_ASSERT_EQUALS(p_coarse_mesh->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(p_coarse_mesh->GetNumElements(), 18u);

        // Find centre of cell population
        c_vector<double,2> centre_of_cell_population = cell_population.GetCentroidOfCellPopulation();

        // Find centre of coarse PDE mesh
        c_vector<double,2> centre_of_coarse_pde_mesh = zero_vector<double>(2);
        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();

        // Test that the two centres match
        c_vector<double,2> centre_diff = centre_of_cell_population - centre_of_coarse_pde_mesh;
        TS_ASSERT_DELTA(norm_2(centre_diff), 0.0, 1e-4);

        // Test that UseCoarsePdeMesh()  throws an exception if the wrong type of PDE is specified
        SimpleUniformSourcePde<2> pde2(-0.1);
        ConstBoundaryCondition<2> bc2(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc2(&pde2, &bc2, false);
        pde_and_bc2.SetDependentVariableName("second variable");
        pde_handler.AddPdeAndBc(&pde_and_bc2);

        TS_ASSERT_THROWS_THIS(pde_handler.UseCoarsePdeMesh(3.0, cuboid, true),
            "UseCoarsePdeMesh() should only be called if averaged-source PDEs are specified.");

        // Now test the 1D case
        std::vector<Node<1>*> nodes_1d;
        nodes_1d.push_back(new Node<1>(0, true,  0.0));

        NodesOnlyMesh<1> mesh_1d;
        mesh_1d.ConstructNodesWithoutMesh(nodes_1d);

        std::vector<CellPtr> cells_1d;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 1> cells_generator_1d;
        cells_generator_1d.GenerateBasic(cells_1d, mesh_1d.GetNumNodes());

        NodeBasedCellPopulation<1> cell_population_1d(mesh_1d, cells_1d);
        cell_population_1d.SetMechanicsCutOffLength(1.5);

        CellBasedPdeHandler<1> pde_handler_1d(&cell_population_1d);

        AveragedSourcePde<1> pde_1d(cell_population_1d, -0.1);
        ConstBoundaryCondition<1> bc_1d(1.0);
        PdeAndBoundaryConditions<1> pde_and_bc_1d(&pde_1d, &bc_1d, false);
        pde_handler_1d.AddPdeAndBc(&pde_and_bc_1d);

        ChastePoint<1> lower1(0.0);
        ChastePoint<1> upper1(9.0);
        ChasteCuboid<1> cuboid1(lower1, upper1);
        pde_handler_1d.UseCoarsePdeMesh(3.0, cuboid1, true);

        // Test that the coarse mesh has the correct number of nodes and elements
        TetrahedralMesh<1,1>* p_coarse_mesh_1d = pde_handler_1d.GetCoarsePdeMesh();
        TS_ASSERT_EQUALS(p_coarse_mesh_1d->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(p_coarse_mesh_1d->GetNumElements(), 3u);

        // Now test the 3D case
        std::vector<Node<3>*> nodes_3d;
        nodes_3d.push_back(new Node<3>(0, true,  0.0));

        NodesOnlyMesh<3> mesh_3d;
        mesh_3d.ConstructNodesWithoutMesh(nodes_3d);

        std::vector<CellPtr> cells_3d;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator_3d;
        cells_generator_3d.GenerateBasic(cells_3d, mesh_3d.GetNumNodes());

        NodeBasedCellPopulation<3> cell_population_3d(mesh_3d, cells_3d);
        cell_population_3d.SetMechanicsCutOffLength(1.5);

        CellBasedPdeHandler<3> pde_handler_3d(&cell_population_3d);

        AveragedSourcePde<3> pde_3d(cell_population_3d, -0.1);
        ConstBoundaryCondition<3> bc_3d(1.0);
        PdeAndBoundaryConditions<3> pde_and_bc_3d(&pde_3d, &bc_3d, false);
        pde_handler_3d.AddPdeAndBc(&pde_and_bc_3d);

        ChastePoint<3> lower3(0.0, 0.0, 0.0);
        ChastePoint<3> upper3(9.0, 9.0, 9.0);
        ChasteCuboid<3> cuboid3(lower3, upper3);

        pde_handler_3d.UseCoarsePdeMesh(3.0, cuboid3, true);

        // Test that the coarse mesh has the correct number of nodes and elements
        TetrahedralMesh<3,3>* p_coarse_mesh_3d = pde_handler_3d.GetCoarsePdeMesh();
        TS_ASSERT_EQUALS(p_coarse_mesh_3d->GetNumNodes(), 64u);
        TS_ASSERT_EQUALS(p_coarse_mesh_3d->GetNumElements(), 162u);

        // Tidy up
        for (unsigned i=0; i<nodes_1d.size(); i++)
        {
            delete nodes_1d[i];
            delete nodes_3d[i];
        }
    }

    void TestUseCoarsePdeMeshNotCentredOnPopulation() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a cell population
        HoneycombMeshGenerator generator(4, 4, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Test that UseCoarsePdeMesh() throws an exception if no PDEs are specified
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(9.0, 9.0);
        ChasteCuboid<2> cuboid(lower, upper);

        // Set up PDE and pass to handler
        AveragedSourcePde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Test UseCoarsePdeMesh() again
        pde_handler.UseCoarsePdeMesh(3.0, cuboid);

        // Check that the centre of the mesh co-incides with centre of the cuboid.
        c_vector<double,2> centre_of_coarse_mesh = zero_vector<double>(2);
        for (unsigned i=0; i<pde_handler.GetCoarsePdeMesh()->GetNumNodes(); i++)
        {
            centre_of_coarse_mesh += pde_handler.GetCoarsePdeMesh()->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_mesh /= pde_handler.GetCoarsePdeMesh()->GetNumNodes();

        c_vector<double,2> centre_of_cuboid = 0.5*(lower.rGetLocation() + upper.rGetLocation());

        TS_ASSERT_DELTA(norm_2(centre_of_cuboid - centre_of_coarse_mesh), 0.0,  1e-4);

    }
    void TestInitialiseCellPdeElementMapAndFindCoarseElementContainingCell() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a cell population
        HoneycombMeshGenerator generator(4, 4, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        TS_ASSERT_EQUALS(pde_handler.mCellPdeElementMap.size(), 0u);

        // Test that InitialiseCellPdeElementMap() throws an exception if mpCoarsePdeMesh is not set up
        TS_ASSERT_THROWS_THIS(pde_handler.InitialiseCellPdeElementMap(),
            "InitialiseCellPdeElementMap() should only be called if mpCoarsePdeMesh is set up.");

        // Set up PDE and pass to handler
        AveragedSourcePde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Use a coarse PDE mesh
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(9.0, 9.0);
        ChasteCuboid<2> cuboid(lower, upper);

        pde_handler.UseCoarsePdeMesh(3.0, cuboid, true);

        // Test that mCellPdeElementMap is initialised correctly
        pde_handler.InitialiseCellPdeElementMap();

        TS_ASSERT_EQUALS(pde_handler.mCellPdeElementMap.size(), cell_population.GetNumRealCells());
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned containing_element_index = pde_handler.mCellPdeElementMap[*cell_iter];
            TS_ASSERT_LESS_THAN(containing_element_index, pde_handler.GetCoarsePdeMesh()->GetNumElements());
            TS_ASSERT_EQUALS(containing_element_index, pde_handler.FindCoarseElementContainingCell(*cell_iter));
        }
    }

    void TestWritingToFile() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        std::string output_directory = "TestCellBasedPdeHandlerWritingToFile";

        // Create a cell population
        HoneycombMeshGenerator generator(4, 4, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        TS_ASSERT_EQUALS(cells.size(), mesh.GetNumNodes());
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        cell_population.SetDataOnAllCells("variable", 1.0);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        TS_ASSERT_THROWS_THIS(pde_handler.OpenResultsFiles(output_directory),
                "Trying to solve a PDE on a cell population that doesn't have a mesh. Try calling UseCoarsePdeMesh().");

        // Use a coarse PDE mesh since we are using a node-based cell population
        AveragedSourcePde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");
        pde_handler.AddPdeAndBc(&pde_and_bc);

        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(9.0, 9.0);
        ChasteCuboid<2> cuboid(lower, upper);

        pde_handler.UseCoarsePdeMesh(3.0, cuboid, true);

        // For coverage, call SetWriteAverageRadialPdeSolution() prior to output
        pde_handler.SetWriteAverageRadialPdeSolution("variable", 5, true);

        // Test that output files are opened correctly
        pde_handler.OpenResultsFiles(output_directory);

        FileFinder file_finder(output_directory + "/results.vizcoarsepdesolution", RelativeTo::ChasteTestOutput);
        TS_ASSERT(file_finder.Exists());
        TS_ASSERT(file_finder.IsFile());

        FileFinder file_finder2(output_directory + "/radial_dist.dat", RelativeTo::ChasteTestOutput);
        TS_ASSERT(file_finder2.Exists());
        TS_ASSERT(file_finder2.IsFile());

        TS_ASSERT_THROWS_NOTHING(pde_handler.CloseResultsFiles());

        // For coverage, also test that output files are opened correctly when not using a coarse PDE mesh
        HoneycombMeshGenerator generator2(5, 5, 0);
        MutableMesh<2,2>* p_mesh2 = generator2.GetMesh();

        std::vector<CellPtr> cells2;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator2;
        cells_generator2.GenerateBasic(cells2, p_mesh2->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population2(*p_mesh2, cells2);

        cell_population2.SetDataOnAllCells("another variable", 1.0);

        CellBasedPdeHandler<2> pde_handler2(&cell_population2);
        pde_handler2.OpenResultsFiles(output_directory);

        FileFinder file_finder3(output_directory + "/results.vizpdesolution", RelativeTo::ChasteTestOutput);
        TS_ASSERT(file_finder3.Exists());
        TS_ASSERT(file_finder3.IsFile());

        pde_handler2.CloseResultsFiles();
    }

    void TestWriteAverageRadialPdeSolution() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a cell population using a circular mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //Put random data on the cells
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("averaged quantity", RandomNumberGenerator::Instance()->ranf());
        }

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        pde_handler.SetWriteAverageRadialPdeSolution("averaged quantity", 2);

        // Open result file ourselves
        OutputFileHandler output_file_handler("TestWriteAverageRadialPdeSolution", false);
        pde_handler.mpAverageRadialPdeSolutionResultsFile = output_file_handler.OpenOutputFile("radial_dist.dat");

        // Write average radial PDE solution to file
        pde_handler.WriteAverageRadialPdeSolution(SimulationTime::Instance()->GetTime());

        // Test that this is correct by comparing with an existing results file
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        NumericFileComparison comparison(results_dir + "/radial_dist.dat", "cell_based/test/data/TestCellBasedPdeHandler/radial_dist.dat");
        TS_ASSERT(comparison.CompareFiles());

        // Close result file ourselves
        pde_handler.mpAverageRadialPdeSolutionResultsFile->close();
    }

    void TestWritePdeSolution() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //Put random data on the cells
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("variable", RandomNumberGenerator::Instance()->ranf());
        }

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Create a single PDE and pass to the handler
        SimpleUniformSourcePde<2> pde(-0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Open result file ourselves
        OutputFileHandler output_file_handler("TestWritePdeSolution", false);
        pde_handler.mpVizPdeSolutionResultsFile = output_file_handler.OpenOutputFile("results.vizpdesolution");

        // Write average radial PDE solution to file
        pde_handler.WritePdeSolution(SimulationTime::Instance()->GetTime());

        // Test that this is correct by comparing with an existing results file
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        NumericFileComparison comparison(results_dir + "/results.vizpdesolution", "cell_based/test/data/TestCellBasedPdeHandler/results.vizpdesolution");
        TS_ASSERT(comparison.CompareFiles());

        // Close result file ourselves
        pde_handler.mpVizPdeSolutionResultsFile->close();
    }

    void TestSolvePdeAndWriteResultsToFileWithoutCoarsePdeMeshDirichlet() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.5, 6);

        // Set up mesh
        MutableMesh<2,2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        cell_population.SetDataOnAllCells("variable", 1.0);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Create a single PDE and pass to the handler
        SimplePdeForTesting pde;
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");

        // For coverage, provide an initial guess for the solution
        std::vector<double> data(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            data[i] = 1.0;
        }

        Vec vector = PetscTools::CreateVec(data);
        pde_and_bc.SetSolution(vector);

        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Solve PDE (set sampling timestep multiple to be large to avoid writing results to file)
        pde_handler.SolvePdeAndWriteResultsToFile(10);

        // Check the correct solution was obtained
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double radius = norm_2(cell_population.GetLocationOfCellCentre(*cell_iter));
            double analytic_solution = 1.0 - 0.25*(1 - pow(radius,2.0));

            // Test that PDE solver is working correctly
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable"), analytic_solution, 0.02);
        }
    }

    void TestSolvePdeAndWriteResultsToFileWithoutCoarsePdeMeshNeumann() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.5, 6);

        // Set up mesh
        MutableMesh<2,2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        cell_population.SetDataOnAllCells("", 1.0);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Create a single PDE and pass to the handler
        CellwiseSourcePde<2> pde(cell_population, 0.0);
        ConstBoundaryCondition<2> bc(0.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, true);


        // For coverage, provide an initial guess for the solution
        std::vector<double> data(mesh.GetNumNodes()+1);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            data[i] = 1.0;
        }

        Vec vector = PetscTools::CreateVec(data);
        pde_and_bc.SetSolution(vector);

        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Solve PDE (set sampling timestep multiple to be large to avoid writing results to file)
        pde_handler.SolvePdeAndWriteResultsToFile(10);

        // Check the correct solution was obtained
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Test that PDE solver is working correctly
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem(""), 0.0, 0.02);
        }
    }

    void TestSolvePdeAndWriteResultsToFileCoarsePdeMeshBoundaryConditions() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.05, 6);

        // Create a cigar-shaped mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(5.0, 1.0);

        // Create a cell population
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        cell_population.SetDataOnAllCells("variable", 1.0);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Set up PDE and pass to handler
        AveragedSourcePde<2> pde(cell_population, -0.01);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");

        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Solve PDEs on a coarse mesh
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(50.0, 50.0);
        ChasteCuboid<2> cuboid(lower, upper);
        pde_handler.UseCoarsePdeMesh(10.0, cuboid, true);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        // For coverage, provide an initial guess for the solution
        std::vector<double> data(pde_handler.mpCoarsePdeMesh->GetNumNodes());
        for (unsigned i=0; i<pde_handler.mpCoarsePdeMesh->GetNumNodes(); i++)
        {
            data[i] = 1.0;
        }

        Vec vector = PetscTools::CreateVec(data);
        pde_and_bc.SetSolution(vector);

        // Solve PDEs (set sampling timestep multiple to be large to avoid writing results to file)
        pde_handler.SolvePdeAndWriteResultsToFile(10);

        // Test that boundary cells experience the right boundary condition
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (cell_population.GetNodeCorrespondingToCell(*cell_iter)->IsBoundaryNode())
            {
                TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable"), 1.0, 1e-1);
            }
        }
    }

    void TestSolvePdeAndWriteResultsToFileWithCoarsePdeMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.05, 6);

        // Create a cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.SetDataOnAllCells("quantity 1", 1.0);
        cell_population.SetDataOnAllCells("quantity 2", 1.0);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Set up PDE and pass to handler
        AveragedSourcePde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("quantity 1");
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Set up second PDE and pass to handler
        AveragedSourcePde<2> pde2(cell_population, -0.5);
        PdeAndBoundaryConditions<2> pde_and_bc2(&pde2, &bc, false);
        TS_ASSERT_THROWS_THIS(pde_handler.AddPdeAndBc(&pde_and_bc2), "When adding more than one PDE to CellBasedPdeHandler set the dependent variable name using SetDependentVariableName(name).");
        pde_and_bc2.SetDependentVariableName("quantity 1");
        TS_ASSERT_THROWS_THIS(pde_handler.AddPdeAndBc(&pde_and_bc2), "The name quantity 1 has already been used in the PDE collection");
        pde_and_bc2.SetDependentVariableName("quantity 2");
        pde_handler.AddPdeAndBc(&pde_and_bc2);

        // Solve PDEs on a coarse mesh
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(50.0, 50.0);
        ChasteCuboid<2> cuboid(lower, upper);
        pde_handler.UseCoarsePdeMesh(10.0, cuboid, true);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        // Open result file ourselves
        OutputFileHandler output_file_handler("TestSolvePdeAndWriteResultsToFileWithCoarsePdeMesh", false);
        pde_handler.mpVizPdeSolutionResultsFile = output_file_handler.OpenOutputFile("results.vizpdesolution");

        // Solve PDEs and (for coverage) write results to file
        pde_handler.SolvePdeAndWriteResultsToFile(1);

        // Close result file ourselves
        pde_handler.mpVizPdeSolutionResultsFile->close();

        TetrahedralMesh<2,2>* p_coarse_mesh = pde_handler.GetCoarsePdeMesh();
        TS_ASSERT(p_coarse_mesh != NULL);

        TS_ASSERT_THROWS_THIS(pde_handler.GetPdeSolution("quantity 3"), "The PDE collection does not contain a PDE named quantity 3");
        ReplicatableVector pde_solution0(pde_handler.GetPdeSolution("quantity 1"));
        ReplicatableVector pde_solution1(pde_handler.GetPdeSolution("quantity 2"));

        TS_ASSERT_EQUALS(pde_solution0.GetSize(), pde_solution1.GetSize());

        // Test that the solution is 1.0 at each coarse mesh node far from the cells
        for (unsigned i=0; i<pde_solution0.GetSize(); i++)
        {
            c_vector<double,2> centre;
            centre(0) = 2.5; // assuming 5x5 honeycomb mesh
            centre(1) = 2.5;

            c_vector<double,2> position = p_coarse_mesh->GetNode(i)->rGetLocation();
            double dist = norm_2(centre - position);
            double u0 = pde_solution0[i];
            double u1 = pde_solution1[i];

            if (dist > 4.0)
            {
                TS_ASSERT_DELTA(u0, 1.0, 1e-5);
                TS_ASSERT_DELTA(u1, 1.0, 1e-5);
            }
        }

        /*
         * Loop over cells, find the coarse mesh element containing it, then
         * check the interpolated PDE solution is between the min and max of
         * the PDE solution on the nodes of that element.
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            c_vector<double,2> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);

            unsigned elem_index = p_coarse_mesh->GetContainingElementIndex(cell_location);
            Element<2,2>* p_element = p_coarse_mesh->GetElement(elem_index);

            unsigned node_0_index = p_element->GetNodeGlobalIndex(0);
            unsigned node_1_index = p_element->GetNodeGlobalIndex(1);
            unsigned node_2_index = p_element->GetNodeGlobalIndex(2);

            double max0 = std::max(pde_solution0[node_0_index], pde_solution0[node_1_index]);
            max0 = std::max(max0, pde_solution0[node_2_index]);

            double max1 = std::max(pde_solution1[node_0_index], pde_solution1[node_1_index]);
            max1 = std::max(max1, pde_solution1[node_2_index]);

            double min0 = std::min(pde_solution0[node_0_index], pde_solution0[node_1_index]);
            min0 = std::min(min0, pde_solution0[node_2_index]);

            double min1 = std::min(pde_solution1[node_0_index], pde_solution1[node_1_index]);
            min1 = std::min(min1, pde_solution1[node_2_index]);

            double value0_at_cell = cell_iter->GetCellData()->GetItem("quantity 1");
            double value1_at_cell = cell_iter->GetCellData()->GetItem("quantity 2");

            TS_ASSERT_LESS_THAN_EQUALS(value1_at_cell, value0_at_cell);
            TS_ASSERT_LESS_THAN_EQUALS(min0, value0_at_cell + DBL_EPSILON);
            TS_ASSERT_LESS_THAN_EQUALS(value0_at_cell, max0 + DBL_EPSILON);
            TS_ASSERT_LESS_THAN_EQUALS(min1, value1_at_cell + DBL_EPSILON);
            TS_ASSERT_LESS_THAN_EQUALS(value1_at_cell, max1 + DBL_EPSILON);
        }
    }

    void TestCoarseSourceMeshWithNeumannIsNotImplemented() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.05, 6);

        // Create a cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.SetDataOnAllCells("first variable", 1.0);
        cell_population.SetDataOnAllCells("second variable", 1.0);

        // Create a PDE handler object using this cell population
        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Set up PDE and pass to handler
        AveragedSourcePde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(0.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, true);
        pde_and_bc.SetDependentVariableName("first variable");
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Set up second PDE and pass to handler
        AveragedSourcePde<2> pde2(cell_population, -0.5);
        ConstBoundaryCondition<2> bc2(0.0);
        PdeAndBoundaryConditions<2> pde_and_bc2(&pde2, &bc2, true);
        pde_and_bc2.SetDependentVariableName("second variable");
        pde_handler.AddPdeAndBc(&pde_and_bc2);

        // Solve PDEs on a coarse mesh
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(50.0, 50.0);
        ChasteCuboid<2> cuboid(lower, upper);
        pde_handler.UseCoarsePdeMesh(10.0, cuboid, true);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        // Test that the correct exception is thrown
        TS_ASSERT_THROWS_THIS(pde_handler.SolvePdeAndWriteResultsToFile(1),
            "Neumann BCs not yet implemented when using a coarse PDE mesh");
    }
};

#endif /*TESTCELLBASEDPDEHANDLER_HPP_*/
