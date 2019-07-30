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

#ifndef TESTELLIPTICGROWINGDOMAINPDEMODIFIER_HPP_
#define TESTELLIPTICGROWINGDOMAINPDEMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

// This macro prevents errors with GCC 4.8 of form "unable to find numeric literal operator 'operator"" Q'"
// when compiling with -std=gnu++11 (see #2929). \todo: remove when GCC 4.8 is no longer supported.
#define BOOST_MATH_DISABLE_FLOAT128
#include <boost/math/special_functions/bessel.hpp>

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "ApoptoticCellProperty.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "ReplicatableVector.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "VertexBasedCellPopulation.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * In this test suite we check the solution of the CellwisePdes against exact solutions.
 * In each case we are solving Laplacian U = f where f is constant in different regions.
 * We solve unit disc where the solutions are Bessel functions and logs.
 */
class TestEllipticGrowingDomainPdeModifier : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestEllipticConstructor()
    {
        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (-0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("averaged quantity");

        // Test that member variables are initialised correctly
        TS_ASSERT_EQUALS(p_pde_modifier->rGetDependentVariableName(), "averaged quantity");
    }

    void TestMeshGeneration()
    {
        // Create PDE and boundary condition objects to be used by all cell populations
        MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (-0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a CellsGenerator to be used by all cell populations
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("averaged quantity");
        {
            // Create a MeshBasedCellPopulation
            HoneycombMeshGenerator generator(10, 10, 0);
            MutableMesh<2,2>* p_mesh = generator.GetMesh();

            std::vector<CellPtr> mesh_cells;
            cells_generator.GenerateBasic(mesh_cells, p_mesh->GetNumNodes());

            MeshBasedCellPopulation<2> mesh_cell_population(*p_mesh, mesh_cells);

            // Now generate the finite element mesh
            p_pde_modifier->GenerateFeMesh(mesh_cell_population);

            // Check that the meshes have the same nodes
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[0], p_mesh->GetNode(i)->rGetLocation()[0], 1e-5);
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[1], p_mesh->GetNode(i)->rGetLocation()[1], 1e-5);
                TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNode(i)->IsBoundaryNode(), p_mesh->GetNode(i)->IsBoundaryNode());
            }
        }

        {
            // Make a NodeBasedCellPopulation
            HoneycombMeshGenerator generator(10, 10, 0);
            MutableMesh<2,2>* p_mesh = generator.GetMesh();
            NodesOnlyMesh<2> node_mesh;
            node_mesh.ConstructNodesWithoutMesh(*p_mesh, 1.5);

            std::vector<CellPtr> node_cells;
            cells_generator.GenerateBasic(node_cells, node_mesh.GetNumNodes());

            NodeBasedCellPopulation<2> node_cell_population(node_mesh, node_cells);

            // Now generate the finite element mesh
            p_pde_modifier->GenerateFeMesh(node_cell_population);

            // Check that the meshes have the same nodes
            for (unsigned i=0; i<node_mesh.GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[0], node_mesh.GetNode(i)->rGetLocation()[0],1e-5);
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[1], node_mesh.GetNode(i)->rGetLocation()[1],1e-5);
            }
        }

        {
            // Make a VertexBasedCellPopulation
            HoneycombVertexMeshGenerator generator(10, 10);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            std::vector<CellPtr> vertex_cells;
            cells_generator.GenerateBasic(vertex_cells, p_mesh->GetNumElements());

            VertexBasedCellPopulation<2> vertex_cell_population(*p_mesh, vertex_cells);

            // Now generate the finite element mesh
            p_pde_modifier->GenerateFeMesh(vertex_cell_population);

            // Check that the meshes have the same nodes
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[0], p_mesh->GetNode(i)->rGetLocation()[0],1e-5);
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[1], p_mesh->GetNode(i)->rGetLocation()[1],1e-5);

                TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNode(i)->IsBoundaryNode(), p_mesh->GetNode(i)->IsBoundaryNode());
            }
            // New node at every element centre
            for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
            {
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i+p_mesh->GetNumNodes())->rGetLocation()[0], p_mesh->GetCentroidOfElement(i)[0],1e-5);
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i+p_mesh->GetNumNodes())->rGetLocation()[1], p_mesh->GetCentroidOfElement(i)[1],1e-5);

                TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNode(i+p_mesh->GetNumNodes())->IsBoundaryNode(), false);
            }
        }

        {
            // Make a PottsBasedCellPopulation
            PottsMeshGenerator<2> generator(50,5,5,50,5,5);
            PottsMesh<2>* p_mesh = generator.GetMesh();

            std::vector<CellPtr> potts_cells;
            cells_generator.GenerateBasic(potts_cells, p_mesh->GetNumElements());

            PottsBasedCellPopulation<2> potts_cell_population(*p_mesh, potts_cells);

            // Now generate the finite element mesh
            p_pde_modifier->GenerateFeMesh(potts_cell_population);

            // Check that the meshes have the same nodes
            for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
            {
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[0], p_mesh->GetCentroidOfElement(i)[0],1e-5);
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[1], p_mesh->GetCentroidOfElement(i)[1],1e-5);
            }
        }

        {
            // Make a CaBasedCellPopulation
            PottsMeshGenerator<2> generator(50,0,0,50,0,0);
            PottsMesh<2>* p_mesh = generator.GetMesh();

            // Specify the location of each cell
            std::vector<unsigned> location_indices;
            for (unsigned i=0; i<10; i++)
            {
                for (unsigned j=0; j<10; j++)
                {
                    unsigned offset = (50+1) * (50-10)/2;
                    location_indices.push_back(offset + j + i * 50);
                }
            }

            std::vector<CellPtr> ca_cells;
            cells_generator.GenerateBasic(ca_cells, location_indices.size());

            // Create cell population
            CaBasedCellPopulation<2> ca_cell_population(*p_mesh, ca_cells, location_indices);

            // Now generate the finite element mesh
            p_pde_modifier->GenerateFeMesh(ca_cell_population);

            // Check that the mesh has nodes at the centre of the cells
            unsigned i=0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = ca_cell_population.Begin();
                 cell_iter != ca_cell_population.End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[0], ca_cell_population.GetLocationOfCellCentre(*cell_iter)[0],1e-5);
                TS_ASSERT_DELTA(p_pde_modifier->mpFeMesh->GetNode(i)->rGetLocation()[1], ca_cell_population.GetLocationOfCellCentre(*cell_iter)[1],1e-5);
                ++i;
            }
        }
    }

    void TestGrowingDomainPdeModifierExceptions()
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructFromMeshReader(mesh_reader);
        temp_mesh.Scale(5.0,1.0);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -1.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("nutrient");

        TS_ASSERT_THROWS_THIS(p_pde_modifier->SetupSolve(cell_population, "output_directory"),
            "EllipticGrowingDomainPdeModifier cannot be used with an AveragedSourceEllipticPde. Use an EllipticBoxDomainPdeModifier instead.");
    }

    void TestArchiveEllipticGrowingDomainPdeModifier()
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "EllipticGrowingDomainPdeModifier.arch";

        // Separate scope to write the archive
        {
            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (-0.1));
            MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

            // Create a PDE modifier and set the name of the dependent variable in the PDE
            std::vector<double> data(10);
            for (unsigned i=0; i<10; i++)
            {
                data[i] = i + 0.45;
            }
            Vec vector = PetscTools::CreateVec(data);
            EllipticGrowingDomainPdeModifier<2> modifier(p_pde, p_bc, false, vector);
            modifier.SetDependentVariableName("averaged quantity");

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            AbstractCellBasedSimulationModifier<2,2>* const p_modifier = &modifier;
            output_arch << p_modifier;
        }

        // Separate scope to read the archive
        {
            AbstractCellBasedSimulationModifier<2,2>* p_modifier2;

            // Restore the modifier
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_modifier2;

            // See whether we read out the correct variable name area
            std::string variable_name = (static_cast<EllipticGrowingDomainPdeModifier<2>*>(p_modifier2))->rGetDependentVariableName();
            TS_ASSERT_EQUALS(variable_name, "averaged quantity");

            Vec solution = (static_cast<EllipticGrowingDomainPdeModifier<2>*>(p_modifier2))->GetSolution();
            ReplicatableVector solution_repl(solution);

            TS_ASSERT_EQUALS(solution_repl.GetSize(), 10u);
            for (unsigned i=0; i<10; i++)
            {
                TS_ASSERT_DELTA(solution_repl[i], i + 0.45, 1e-6);
            }

            delete p_modifier2;
        }
    }

    /*
     * Here the exact solution is
     *
     * u = J0(r)/J0(1)
     *
     * where J0 is the zeroth order Bessel function.
     */
    void TestMeshBasedMonolayerWithEllipticPde()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, 1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        // For coverage output the solution gradient
        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseEllipticPdeWithMeshOnDisk");

        // Test the solution against the exact solution
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            c_vector<double,2> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double r = sqrt(cell_location(0)*cell_location(0) + cell_location(1)*cell_location(1));
            double u_exact = boost::math::cyl_bessel_j(0,r) / boost::math::cyl_bessel_j(0,1);

            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable"), u_exact, 1e-3);

            double ux_exact = -cell_location(0) / r * boost::math::cyl_bessel_j(1,r) / boost::math::cyl_bessel_j(0,1);
            double uy_exact = -cell_location(1) / r * boost::math::cyl_bessel_j(1,r) / boost::math::cyl_bessel_j(0,1);

            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable_grad_x"), ux_exact, 1e-1);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable_grad_y"), uy_exact, 1e-1);
        }
    }

    /*
     * Here the outer cells (r>1/2) are apoptotic so the solution is
     *
     * u = C*J0(r)  for r in [0,0.5]
     *     A*ln(r) + 1 for r in [0.5,1]
     *
     *  where J0 is the zeroth order bessel fn and C and A are constants.
     */
    void TestMeshBasedHeterogeneousMonolayerWithEllipticPde()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Make cells with r<1/2 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = mesh.GetNode(i)->rGetLocation();
            double r = sqrt(cell_location(0)*cell_location(0) + cell_location(1)*cell_location(1));
            if (r > 0.5)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, 1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseEllipticPdeWithMeshOnHeterogeneousDisk");

        // Test the solution against the exact solution
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {

            c_vector<double,2> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double r = sqrt(cell_location(0)*cell_location(0) + cell_location(1)*cell_location(1));

            double J005 = boost::math::cyl_bessel_j(0,0.5);
            double J105 = boost::math::cyl_bessel_j(1,0.5);

            double A = -1.0/(2.0*J005/J105 +log(0.5));
            double C = -2*A/J105;

            double u_exact = C*boost::math::cyl_bessel_j(0,r);
            if (r > 0.5)
            {
                u_exact = A*log(r)+1.0;
            }

            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable"), u_exact, 1e-2);
        }
    }

    // Now test on a square with half apoptotic cells to compare all the population types
    void TestMeshBasedSquareMonolayer()
    {
        HoneycombMeshGenerator generator(20,20,0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 10.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 200u);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseEllipticPdeWithMeshOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.4542, 1e-4);
    }

    void TestNodeBasedSquareMonolayer()
    {
        HoneycombMeshGenerator generator(20,20,0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 10.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),200u);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseEllipticPdeWithNodeOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA( p_cell_210->GetCellData()->GetItem("variable"), 0.4542, 1e-2); // Lower tolerance as slightly different meshes

        // Checking it doesn't change for this cell population
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.4476, 1e-4);

        // Clear memory
        delete p_mesh;
    }

    void TestVertexBasedSquareMonolayer()
    {
        HoneycombVertexMeshGenerator generator(20,20);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        p_mesh->Translate(-0.5,-sqrt(3.0)/3); // Shift so cells are on top of those in the above centre-based tests

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0) < 10.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),200u);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseEllipticPdeWithVertexOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.4542, 1e-1); // Low tolerance as mesh is slightly larger than for centre based models

        // Checking it doesn't change for this cell population
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.4654, 1e-4);
    }

    void TestPottsBasedSquareMonolayer()
    {
        PottsMeshGenerator<2> generator(100, 20, 4, 100, 20, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Translate and scale so cells are on top of those in the above centre based tests
        p_mesh->Translate(-11.5,-11.5);
        p_mesh->Scale(0.25,0.25 *sqrt(3.0)*0.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0) < 10.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),200u);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseEllipticPdeWithPottsOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA( p_cell_210->GetCellData()->GetItem("variable"), 0.4542, 2e-1); // Low tolerance as mesh is slightly larger than for centre based models

        // Check it doesn't change for this cell population
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.4338, 1e-4);
    }

    void TestCaBasedSquareMonolayer()
    {
        PottsMeshGenerator<2> generator(20, 0, 0, 20, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Scale so cells are on top of those in the above centre based tests
        p_mesh->Scale(1.0,sqrt(3.0)*0.5);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<400; i++)
        {
            location_indices.push_back(i);
        }

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 10.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),200u);

        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseEllipticPdeWithCaOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.4542, 2e-1); // Low tolerance as mesh is slightlty larger than for centre-based models

        // Check it doesn't change for this cell population
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.4338, 1e-3); // Note lower as slightly different answer with intel compiler
    }

    void TestEllipticGrowingDomainPdeModifierIn1d()
    {
        // Create mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, true,  0.0));
        nodes.push_back(new Node<1>(1, false, 1.0));
        nodes.push_back(new Node<1>(2, false, 2.0));
        nodes.push_back(new Node<1>(3, false, 3.0));
        nodes.push_back(new Node<1>(4, true,  4.0));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Set up cell population
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<1>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<1>, p_bc, (1.0));;

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<1>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestEllipticBoxDomainPdeModifierIn1d");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_2 = cell_population.GetCellUsingLocationIndex(2);
        TS_ASSERT_DELTA(p_cell_2->GetCellData()->GetItem("variable"), 0.8274, 1e-2);
        TS_ASSERT_DELTA(p_cell_2->GetCellData()->GetItem("variable_grad_x"), 0.0, 1e-2);

        // Clear memory
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestEllipticGrowingDomainPdeModifierIn3d()
    {
        // Create a simple mesh
        TetrahedralMesh<3,3> temp_mesh;
        temp_mesh.ConstructCuboid(3,3,3);
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Set up cell population
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<3>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestEllipticGrowingDomainPdeModifierIn3d");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_62 = cell_population.GetCellUsingLocationIndex(13);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable"), 1.0, 1e-2);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable_grad_x"), 0.0, 1e-2);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable_grad_y"), 0.0, 1e-2);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable_grad_z"), 0.0, 1e-2);
    }
};

#endif /*TESTELLIPTICGROWINGDOMAINPDEMODIFIER_HPP_*/
