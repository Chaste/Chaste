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

#ifndef TESTCELLBASEDELLIPTICPDES_HPP_
#define TESTCELLBASEDELLIPTICPDES_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "VolumeDependentAveragedSourceEllipticPde.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCellBasedEllipticPdes : public AbstractCellBasedTestSuite
{
public:

    void TestUniformSourceEllipticPde()
    {
        // Create PDE object
        UniformSourceEllipticPde<2> pde(0.05);

        // Test that the member variable has been initialised correctly
        TS_ASSERT_DELTA(pde.GetCoefficient(), 0.05, 1e-6);

        // Test ComputeDiffusionTerm() method
        ChastePoint<2> point;
        c_matrix<double,2,2> diffusion_matrix = pde.ComputeDiffusionTerm(point);
        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<2; j++)
            {
                double value = 0.0;
                if (i == j)
                {
                    value = 1.0;
                }
                TS_ASSERT_DELTA(diffusion_matrix(i,j), value, 1e-6);
            }
        }

        // Test ComputeConstantInUSourceTerm() method
        TS_ASSERT_DELTA(pde.ComputeConstantInUSourceTerm(point, NULL), 0.0, 1e-6);

        // Test ComputeLinearInUCoeffInSourceTerm() method
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        TS_ASSERT_DELTA(pde.ComputeLinearInUCoeffInSourceTerm(point, NULL), 0.05, 1e-6);
        TS_ASSERT_DELTA(pde.ComputeLinearInUCoeffInSourceTerm(point, p_mesh->GetElement(0)), 0.05, 1e-6);
    }

    void TestUniformSourceEllipticPdeArchiving()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "UniformSourceEllipticPde.arch";
        ArchiveLocationInfo::SetMeshFilename("UniformSourceEllipticPde");

        {
            // Create a PDE object
            AbstractLinearEllipticPde<2,2>* const p_pde = new UniformSourceEllipticPde<2>(0.05);

            // Create output archive and archive PDE object
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) << p_pde;

            delete p_pde;
        }

        {
            AbstractLinearEllipticPde<2,2>* p_pde;

            // Create an input archive and restore PDE object from archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) >> p_pde;

            // Test that the PDE and its member variable were archived correctly
            TS_ASSERT(dynamic_cast<UniformSourceEllipticPde<2>*>(p_pde) != NULL);

            UniformSourceEllipticPde<2>* p_static_cast_pde = static_cast<UniformSourceEllipticPde<2>*>(p_pde);
            TS_ASSERT_DELTA(p_static_cast_pde->GetCoefficient(), 0.05, 1e-6);

            delete p_pde;
        }
    }

    void TestCellwiseSourceEllipticPde()
    {
        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Make one cell apoptotic
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);
        cells[0]->AddCellProperty(p_apoptotic_state);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a PDE object
        CellwiseSourceEllipticPde<2> pde(cell_population, 0.05);

        // Test that the member variables have been initialised correctly
        TS_ASSERT_EQUALS(&(pde.rGetCellPopulation()), &cell_population);
        TS_ASSERT_DELTA(pde.GetCoefficient(), 0.05, 1e-6);

        // Test ComputeDiffusionTerm() method
        ChastePoint<2> point;
        c_matrix<double,2,2> diffusion_matrix = pde.ComputeDiffusionTerm(point);
        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<2; j++)
            {
                double value = 0.0;
                if (i == j)
                {
                    value = 1.0;
                }
                TS_ASSERT_DELTA(diffusion_matrix(i,j), value, 1e-6);
            }
        }

        // Test ComputeConstantInUSourceTerm() method
        TS_ASSERT_DELTA(pde.ComputeConstantInUSourceTerm(point, NULL), 0.0, 1e-6);

        Node<2>* p_node_0 = cell_population.GetNodeCorrespondingToCell(cell_population.GetCellUsingLocationIndex(0));
        TS_ASSERT_DELTA(pde.ComputeLinearInUCoeffInSourceTermAtNode(*p_node_0), 0.0, 1e-6);

        Node<2>* p_node_1 = cell_population.GetNodeCorrespondingToCell(cell_population.GetCellUsingLocationIndex(1));
        TS_ASSERT_DELTA(pde.ComputeLinearInUCoeffInSourceTermAtNode(*p_node_1), 0.05, 1e-6);
    }

    void TestCellwiseSourceEllipticPdeArchiving()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "CellwiseSourceEllipticPde.arch";
        ArchiveLocationInfo::SetMeshFilename("CellwiseSourceEllipticPde");

        {
            // Create a PDE object
            AbstractLinearEllipticPde<2,2>* const p_pde = new CellwiseSourceEllipticPde<2>(cell_population, 0.05);

            // Create output archive and archive PDE object
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) << p_pde;

            delete p_pde;
        }

        {
            AbstractLinearEllipticPde<2,2>* p_pde;

            // Create an input archive and restore PDE object from archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) >> p_pde;

            // Test that the PDE and its member variables were archived correctly
            TS_ASSERT(dynamic_cast<CellwiseSourceEllipticPde<2>*>(p_pde) != NULL);

            CellwiseSourceEllipticPde<2>* p_static_cast_pde = static_cast<CellwiseSourceEllipticPde<2>*>(p_pde);
            TS_ASSERT_DELTA(p_static_cast_pde->GetCoefficient(), 0.05, 1e-6);
            TS_ASSERT_EQUALS(p_static_cast_pde->mrCellPopulation.GetNumRealCells(), 25u);

            // Avoid memory leaks
            delete &(p_static_cast_pde->mrCellPopulation);
            delete p_pde;
        }
    }

    void TestAveragedSourceEllipticPde()
    {
        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // For simplicity we create a very large coarse mesh, so we know that all cells are contained in one element
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructFromMeshReader(mesh_reader);
        coarse_mesh.Scale(10.0, 10.0);

        // Create a PDE object
        AveragedSourceEllipticPde<2> pde(cell_population, 0.05);

        // Test that the member variables have been initialised correctly
        TS_ASSERT_EQUALS(&(pde.rGetCellPopulation()), &cell_population);
        TS_ASSERT_DELTA(pde.GetCoefficient(), 0.05, 1e-6);

        // Test ComputeDiffusionTerm() method
        ChastePoint<2> point;
        c_matrix<double,2,2> diffusion_matrix = pde.ComputeDiffusionTerm(point);
        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<2; j++)
            {
                double value = 0.0;
                if (i == j)
                {
                    value = 1.0;
                }
                TS_ASSERT_DELTA(diffusion_matrix(i,j), value, 1e-6);
            }
        }

        // Test SetupSourceTerms() when no map between cells and coarse mesh elements is supplied
        pde.SetupSourceTerms(coarse_mesh);

        TS_ASSERT_EQUALS(pde.mCellDensityOnCoarseElements.size(), 2u);

        // The first element has area 0.5*10*10 = 50 and there are 5*5 = 25 cells, so the cell density is 25/50 = 0.5
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[0], 0.5, 1e-6);

        // The first element doesn't contain any cells, so the cell density is zero
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[1], 0.0, 1e-6);

        // Now test SetupSourceTerms() when a map between cells and coarse mesh elements is supplied
        std::map<CellPtr, unsigned> cell_pde_element_map;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_pde_element_map[*cell_iter] = 0;
        }

        pde.SetupSourceTerms(coarse_mesh, &cell_pde_element_map);
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[1], 0.0, 1e-6);

        // Test ComputeLinearInUCoeffInSourceTerm() method
        TS_ASSERT_DELTA(pde.ComputeLinearInUCoeffInSourceTerm(point,coarse_mesh.GetElement(0)), 0.05*0.5, 1e-6);

        // Test ComputeConstantInUSourceTerm() method
        TS_ASSERT_DELTA(pde.ComputeConstantInUSourceTerm(point, NULL), 0.0, 1e-6);

        // Test GetUptakeRateForElement()
        TS_ASSERT_DELTA(pde.GetUptakeRateForElement(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(pde.GetUptakeRateForElement(1), 0.0, 1e-6);
    }

    void TestAveragedSourceEllipticPdeArchiving()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "AveragedSourceEllipticPde.arch";
        ArchiveLocationInfo::SetMeshFilename("AveragedSourceEllipticPde");

        {
            // Create a PDE object
            AbstractLinearEllipticPde<2,2>* const p_pde = new AveragedSourceEllipticPde<2>(cell_population, 0.05);

            // Create output archive and archive PDE object
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) << p_pde;

            delete p_pde;
        }

        {
            AbstractLinearEllipticPde<2,2>* p_pde;

            // Create an input archive and restore PDE object from archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) >> p_pde;

            // Test that the PDE and its member variables were archived correctly
            TS_ASSERT(dynamic_cast<AveragedSourceEllipticPde<2>*>(p_pde) != NULL);

            AveragedSourceEllipticPde<2>* p_static_cast_pde = static_cast<AveragedSourceEllipticPde<2>*>(p_pde);
            TS_ASSERT_DELTA(p_static_cast_pde->GetCoefficient(), 0.05, 1e-6);
            TS_ASSERT_EQUALS(p_static_cast_pde->mrCellPopulation.GetNumRealCells(), 25u);

            // Avoid memory leaks
            delete &(p_static_cast_pde->rGetCellPopulation());
            delete p_pde;
        }
    }

    void TestVolumeDependentAveragedSourceEllipticPde()
    {
        // Create a cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create a PDE object
        VolumeDependentAveragedSourceEllipticPde<2> pde(cell_population, 0.05);

        // Test that the member variables have been initialised correctly
        TS_ASSERT(pde.mpStaticCastCellPopulation != NULL);

        // For simplicity we create a very large coarse mesh, so we know that all cells are contained in one element
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructFromMeshReader(mesh_reader);
        coarse_mesh.Scale(10.0, 10.0);

        // Test SetupSourceTerms() when no map between cells and coarse mesh elements is supplied
        pde.SetupSourceTerms(coarse_mesh);

        TS_ASSERT_EQUALS(pde.mCellDensityOnCoarseElements.size(), 2u);

        // The first element has area 0.5*10*10 = 50 and there are 5*5 = 25 cells, so the cell density is 25/50 = 0.5
        // The radius of each node is 0.5 and the density is normalized with cell area (1/4.0) in `VolumeDependentAveragedSourceEllipticPde`.
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[0], 0.5/4.0, 1e-6);

        // The first element doesn't contain any cells, so the cell density is zero
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[1], 0.0, 1e-6);

        // Now test SetupSourceTerms() when a map between cells and coarse mesh elements is supplied
        std::map<CellPtr, unsigned> cell_pde_element_map;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_pde_element_map[*cell_iter] = 0;
        }

        pde.SetupSourceTerms(coarse_mesh, &cell_pde_element_map);
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[0], 0.5/4.0, 1e-6);
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[1], 0.0, 1e-6);
    }

    void TestVolumeDependentAveragedSourceEllipticPdeArchiving()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "VolumeDependentAveragedSourceEllipticPde.arch";
        ArchiveLocationInfo::SetMeshFilename("VolumeDependentAveragedSourceEllipticPde");

        {
            // Create a PDE object
            AbstractLinearEllipticPde<2,2>* const p_pde = new VolumeDependentAveragedSourceEllipticPde<2>(cell_population, 0.05);

            // Create output archive and archive PDE object
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) << p_pde;

            // Avoid memory leak
            delete p_pde;
        }

        {
            AbstractLinearEllipticPde<2,2>* p_pde;

            // Create an input archive and restore PDE object from archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) >> p_pde;

            // Test that the PDE and its member variables were archived correctly
            TS_ASSERT(dynamic_cast<VolumeDependentAveragedSourceEllipticPde<2>*>(p_pde) != NULL);

            VolumeDependentAveragedSourceEllipticPde<2>* p_static_cast_pde = static_cast<VolumeDependentAveragedSourceEllipticPde<2>*>(p_pde);
            TS_ASSERT_DELTA(p_static_cast_pde->GetCoefficient(), 0.05, 1e-6);
            TS_ASSERT_EQUALS(p_static_cast_pde->mrCellPopulation.GetNumRealCells(), 25u);
            TS_ASSERT(p_static_cast_pde->mpStaticCastCellPopulation != NULL);
            TS_ASSERT_EQUALS(p_static_cast_pde->mpStaticCastCellPopulation->GetNumRealCells(), 25u);

            // Avoid memory leaks
            delete &(p_static_cast_pde->mrCellPopulation);
            delete p_pde;
        }
    }
};

#endif /*TESTCELLBASEDELLIPTICPDES_HPP_*/
