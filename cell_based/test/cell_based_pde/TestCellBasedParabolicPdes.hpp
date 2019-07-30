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

#ifndef TESTCELLBASEDPARABOLICPDES_HPP_
#define TESTCELLBASEDPARABOLICPDES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AbstractCellBasedTestSuite.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * This test suite covers parabolic PDE classes.
 */
class TestCellBasedParabolicPdes : public AbstractCellBasedTestSuite
{
public:

    void TestCellwiseSourceParabolicPdeMethods()
    {
        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a PDE object
        CellwiseSourceParabolicPde<2> pde(cell_population, 0.1, 0.2, 0.3);

        // Test that the member variables have been initialised correctly
        // NOTE most member variables are tested below in the Methods section
        TS_ASSERT_EQUALS(&(pde.rGetCellPopulation()), &cell_population);

        // Test methods
        ChastePoint<2> point;
        Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*(cell_population.Begin()));
        TS_ASSERT_DELTA(pde.ComputeSourceTermAtNode(*p_node,2.0), 0.6, 1e-6);
        TS_ASSERT_DELTA(pde.ComputeDuDtCoefficientFunction(point), 0.1, 1e-6);
        c_matrix<double,2,2> diffusion_matrix = pde.ComputeDiffusionTerm(point);
        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<2; j++)
            {
                double value = 0.0;
                if (i == j)
                {
                    value = 0.2;
                }
                TS_ASSERT_DELTA(diffusion_matrix(i,j), value, 1e-6);
            }
        }
    }

    void TestCellwiseSourceParabolicPdeArchiving()
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
        std::string archive_file = "CellwiseSourceParabolicPde.arch";
        ArchiveLocationInfo::SetMeshFilename("CellwiseSourceParabolicPde");

        {
            // Create a PDE object
            AbstractLinearParabolicPde<2,2>* const p_pde = new CellwiseSourceParabolicPde<2>(cell_population, 0.1, 0.2, 0.3);

            // Create output archive and archive PDE object
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) << p_pde;

            delete p_pde;
        }

        {
            AbstractLinearParabolicPde<2,2>* p_pde;

            // Create an input archive and restore PDE object from archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) >> p_pde;

            // Test that the PDE and its member variables were archived correctly
            TS_ASSERT(dynamic_cast<CellwiseSourceParabolicPde<2>*>(p_pde) != NULL);

            CellwiseSourceParabolicPde<2>* p_static_cast_pde = static_cast<CellwiseSourceParabolicPde<2>*>(p_pde);
            TS_ASSERT_EQUALS(p_static_cast_pde->mrCellPopulation.GetNumRealCells(), 25u);

            ChastePoint<2> point;

            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*(cell_population.Begin()));
            TS_ASSERT_DELTA(p_static_cast_pde->ComputeSourceTermAtNode(*p_node,2.0), 0.6, 1e-6);
            TS_ASSERT_DELTA(p_static_cast_pde->ComputeDuDtCoefficientFunction(point), 0.1, 1e-6);
            c_matrix<double,2,2> diffusion_matrix = p_static_cast_pde->ComputeDiffusionTerm(point);
            for (unsigned i=0; i<2; i++)
            {
                for (unsigned j=0; j<2; j++)
                {
                    double value = 0.0;
                    if (i == j)
                    {
                        value = 0.2;
                    }
                    TS_ASSERT_DELTA(diffusion_matrix(i,j), value, 1e-6);
                }
            }

            // Avoid memory leaks
            delete &(p_static_cast_pde->mrCellPopulation);
            delete p_pde;
        }
    }

    void TestUniformSourceParabolicPdeMethods()
    {
        // Create a PDE object
        UniformSourceParabolicPde<2> pde(0.1);

        // Test that the member variables have been initialised correctly
        TS_ASSERT_EQUALS(pde.GetCoefficient(),0.1);

        // Test methods
        ChastePoint<2> point;
        TS_ASSERT_DELTA(pde.ComputeSourceTerm(point,DBL_MAX), 0.1, 1e-6);
        TS_ASSERT_DELTA(pde.ComputeDuDtCoefficientFunction(point), 1.0, 1e-6);
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

    void TestUniformSourceParabolicPdeArchiving()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "UniformSourceParabolicPde.arch";
        ArchiveLocationInfo::SetMeshFilename("UniformSourceParabolicPde");

        {
            // Create a PDE object
            AbstractLinearParabolicPde<2,2>* const p_pde = new UniformSourceParabolicPde<2>(0.1);

            // Create output archive and archive PDE object
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) << p_pde;

            delete p_pde;
        }

        {
            AbstractLinearParabolicPde<2,2>* p_pde;

            // Create an input archive and restore PDE object from archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) >> p_pde;

            // Test that the PDE and its member variables were archived correctly
            TS_ASSERT(dynamic_cast<UniformSourceParabolicPde<2>*>(p_pde) != NULL);

            UniformSourceParabolicPde<2>* p_static_cast_pde = static_cast<UniformSourceParabolicPde<2>*>(p_pde);
            TS_ASSERT_EQUALS(p_static_cast_pde->GetCoefficient(),0.1);

            // Avoid memory leaks
            delete p_pde;
        }
    }

    void TestAveragedSourceParabolicPdeMethods()
    {
        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a PDE object
        AveragedSourceParabolicPde<2> pde(cell_population, 0.1, 0.2, 0.3);

        // Test that the member variables have been initialised correctly
        // NOTE most member variables are tested below in the Methods section
        TS_ASSERT_EQUALS(&(pde.rGetCellPopulation()), &cell_population);

        // For simplicity we create a very large coarse mesh, so we know that all cells are contained in one element
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> fe_mesh;
        fe_mesh.ConstructFromMeshReader(mesh_reader);
        fe_mesh.Scale(10.0, 10.0);

        // Test SetupSourceTerms() when no map between cells and coarse mesh elements is supplied
        pde.SetupSourceTerms(fe_mesh);

        TS_ASSERT_EQUALS(pde.mCellDensityOnCoarseElements.size(), 2u);

        // The first element has area 0.5*10*10 = 50 and there are 5*5 = 25 cells, so the cell density is 25/50 = 0.5
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[0], 0.5, 1e-6);

        // The first element doesn't contain any cells, so the cell density is zero
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[1], 0.0, 1e-6);

        // Now test SetupSourceTerms() when a map between cells and FE mesh elements is supplied
        std::map<CellPtr, unsigned> cell_pde_element_map;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_pde_element_map[*cell_iter] = 0;
        }

        pde.SetupSourceTerms(fe_mesh, &cell_pde_element_map);
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(pde.mCellDensityOnCoarseElements[1], 0.0, 1e-6);

        // Test GetUptakeRateForElement()
        TS_ASSERT_DELTA(pde.GetUptakeRateForElement(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(pde.GetUptakeRateForElement(1), 0.0, 1e-6);

        // Test other methods
        ChastePoint<2> point;
        TS_ASSERT_DELTA(pde.ComputeSourceTerm(point, 1.0, fe_mesh.GetElement(0)), 0.3*0.5, 1e-6);
        TS_ASSERT_DELTA(pde.ComputeSourceTerm(point, 1.0, fe_mesh.GetElement(1)), 0.3*0.0, 1e-6);

        TS_ASSERT_DELTA(pde.ComputeDuDtCoefficientFunction(point), 0.1, 1e-6);

        c_matrix<double,2,2> diffusion_matrix = pde.ComputeDiffusionTerm(point);
        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<2; j++)
            {
                double value = 0.0;
                if (i == j)
                {
                    value = 0.2;
                }
                TS_ASSERT_DELTA(diffusion_matrix(i,j), value, 1e-6);
            }
        }
    }

    void TestAveragedSourceParabolicPdeArchiving()
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
        std::string archive_file = "AveragedSourceParabolicPde.arch";
        ArchiveLocationInfo::SetMeshFilename("AveragedSourceParabolicPde");

        {
            // Create a PDE object
            AbstractLinearParabolicPde<2,2>* const p_pde = new AveragedSourceParabolicPde<2>(cell_population, 0.1, 0.2, 0.3);

            // Create output archive and archive PDE object
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) << p_pde;

            delete p_pde;
        }

        {
            AbstractLinearParabolicPde<2,2>* p_pde;

            // Create an input archive and restore PDE object from archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) >> p_pde;

            // Test that the PDE and its member variables were archived correctly
            TS_ASSERT(dynamic_cast<AveragedSourceParabolicPde<2>*>(p_pde) != NULL);

            AveragedSourceParabolicPde<2>* p_static_cast_pde = static_cast<AveragedSourceParabolicPde<2>*>(p_pde);
            TS_ASSERT_EQUALS(p_static_cast_pde->mrCellPopulation.GetNumRealCells(), 25u);

            ChastePoint<2> point;

            TS_ASSERT_DELTA(p_static_cast_pde->mSourceCoefficient, 0.3 , 1e-6);
            TS_ASSERT_DELTA(p_static_cast_pde->ComputeDuDtCoefficientFunction(point), 0.1, 1e-6);
            c_matrix<double,2,2> diffusion_matrix = p_static_cast_pde->ComputeDiffusionTerm(point);
            for (unsigned i=0; i<2; i++)
            {
                for (unsigned j=0; j<2; j++)
                {
                    double value = 0.0;
                    if (i == j)
                    {
                        value = 0.2;
                    }
                    TS_ASSERT_DELTA(diffusion_matrix(i,j), value, 1e-6);
                }
            }

            // Avoid memory leaks
            delete &(p_static_cast_pde->mrCellPopulation);
            delete p_pde;
        }
    }
};

#endif /*TESTCELLBASEDPARABOLICPDES_HPP_*/
