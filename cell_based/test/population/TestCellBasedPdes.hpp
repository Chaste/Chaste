/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTCELLBASEDPDES_HPP_
#define TESTCELLBASEDPDES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"
#include "SimpleUniformSourcePde.hpp"
#include "CellwiseSourcePde.hpp"
#include "AveragedSourcePde.hpp"
#include "VolumeDependentAveragedSourcePde.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/**
 * This test suite covers any PDE classes defined in cell_based/src/population/pdes.
 */
class TestCellBasedPdes : public AbstractCellBasedTestSuite
{
public:

    void TestSimpleUniformSourcePdeMethods() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a PDE object
        SimpleUniformSourcePde<2> pde(0.05);

        // Test that the member variables have been initialised correctly
        TS_ASSERT_DELTA(pde.GetCoefficient(), 0.05, 1e-6);

        ChastePoint<2> point;

        TS_ASSERT_DELTA(pde.ComputeConstantInUSourceTerm(point, NULL), 0.0, 1e-6);
        TS_ASSERT_DELTA(pde.ComputeLinearInUCoeffInSourceTerm(point, NULL), 0.05, 1e-6);

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
    }

    void TestSimpleUniformSourcePdeArchiving() throw(Exception)
    {
        EXIT_IF_PARALLEL; // HoneycombMeshGenerator doesn't work in parallel

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "SimpleUniformSourcePde.arch";
        ArchiveLocationInfo::SetMeshFilename("SimpleUniformSourcePde");

        {
            // Create a PDE object
            AbstractLinearEllipticPde<2,2>* const p_pde = new SimpleUniformSourcePde<2>(0.05);

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
            TS_ASSERT(dynamic_cast<SimpleUniformSourcePde<2>*>(p_pde) != NULL);

            SimpleUniformSourcePde<2>* p_static_cast_pde = static_cast<SimpleUniformSourcePde<2>*>(p_pde);
            TS_ASSERT_DELTA(p_static_cast_pde->GetCoefficient(), 0.05, 1e-6);

            delete p_pde;
        }
    }

    void TestCellwiseSourcePdeMethods() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a PDE object
        CellwiseSourcePde<2> pde(cell_population, 0.05);

        // Test that the member variables have been initialised correctly
        TS_ASSERT_EQUALS(&(pde.rGetCellPopulation()), &cell_population);
        TS_ASSERT_DELTA(pde.GetCoefficient(), 0.05, 1e-6);

        // Test methods
        Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*(cell_population.Begin()));
        TS_ASSERT_DELTA(pde.ComputeLinearInUCoeffInSourceTermAtNode(*p_node), 0.05, 1e-6);

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
    }

    void TestCellwiseSourcePdeArchiving() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "CellwiseSourcePde.arch";
        ArchiveLocationInfo::SetMeshFilename("CellwiseSourcePde");

        {
            // Create a PDE object
            AbstractLinearEllipticPde<2,2>* const p_pde = new CellwiseSourcePde<2>(cell_population, 0.05);

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
            TS_ASSERT(dynamic_cast<CellwiseSourcePde<2>*>(p_pde) != NULL);

            CellwiseSourcePde<2>* p_static_cast_pde = static_cast<CellwiseSourcePde<2>*>(p_pde);
            TS_ASSERT_DELTA(p_static_cast_pde->GetCoefficient(), 0.05, 1e-6);
            TS_ASSERT_EQUALS(p_static_cast_pde->mrCellPopulation.GetNumRealCells(), 25u);

            // Avoid memory leaks
            delete &(p_static_cast_pde->mrCellPopulation);
            delete p_pde;
        }
    }

    void TestAveragedSourcePdeMethods() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a PDE object
        AveragedSourcePde<2> pde(cell_population, 0.05);

        // Test that the member variables have been initialised correctly
        TS_ASSERT_EQUALS(&(pde.rGetCellPopulation()), &cell_population);
        TS_ASSERT_DELTA(pde.GetCoefficient(), 0.05, 1e-6);

        // For simplicity we create a very large coarse mesh, so we know that all cells are contained in one element
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructFromMeshReader(mesh_reader);
        coarse_mesh.Scale(10.0, 10.0);

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

        // Test GetUptakeRateForElement()
        TS_ASSERT_DELTA(pde.GetUptakeRateForElement(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(pde.GetUptakeRateForElement(1), 0.0, 1e-6);

        // Test other methods
        ChastePoint<2> point;
        TS_ASSERT_DELTA(pde.ComputeLinearInUCoeffInSourceTerm(point, coarse_mesh.GetElement(0)), 0.05*0.5, 1e-6);

        TS_ASSERT_DELTA(pde.ComputeConstantInUSourceTerm(point, NULL), 0.0, 1e-6);

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
    }

    void TestAveragedSourcePdeArchiving() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "AveragedSourcePde.arch";
        ArchiveLocationInfo::SetMeshFilename("AveragedSourcePde");

        {
            // Create a PDE object
            AbstractLinearEllipticPde<2,2>* const p_pde = new AveragedSourcePde<2>(cell_population, 0.05);

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
            TS_ASSERT(dynamic_cast<AveragedSourcePde<2>*>(p_pde) != NULL);

            AveragedSourcePde<2>* p_static_cast_pde = static_cast<AveragedSourcePde<2>*>(p_pde);
            TS_ASSERT_DELTA(p_static_cast_pde->GetCoefficient(), 0.05, 1e-6);
            TS_ASSERT_EQUALS(p_static_cast_pde->mrCellPopulation.GetNumRealCells(), 25u);

            // Avoid memory leaks
            delete &(p_static_cast_pde->mrCellPopulation);
            delete p_pde;
        }
    }

    void TestVolumeDependentAveragedSourcePdeMethods() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create a PDE object
        VolumeDependentAveragedSourcePde<2> pde(cell_population, 0.05);

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
    }

    void TestVolumeDependentAveragedSourcePdeArchiving() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "VolumeDependentAveragedSourcePde.arch";
        ArchiveLocationInfo::SetMeshFilename("VolumeDependentAveragedSourcePde");

        {
            // Create a PDE object
            AbstractLinearEllipticPde<2,2>* const p_pde = new VolumeDependentAveragedSourcePde<2>(cell_population, 0.05);

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
            TS_ASSERT(dynamic_cast<VolumeDependentAveragedSourcePde<2>*>(p_pde) != NULL);

            VolumeDependentAveragedSourcePde<2>* p_static_cast_pde = static_cast<VolumeDependentAveragedSourcePde<2>*>(p_pde);
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

#endif /*TESTCELLBASEDPDES_HPP_*/
