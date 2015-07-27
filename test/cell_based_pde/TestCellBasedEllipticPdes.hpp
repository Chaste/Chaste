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

#ifndef TESTCELLBASEDELLIPTICPDES_HPP_
#define TESTCELLBASEDELLIPTICPDES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * This test suite covers any Elliptic PDE classes defined in cell_based_pdes/src.
 */
class TestCellBasedEllipticPdes : public AbstractCellBasedTestSuite
{
public:


    void TestCellwiseSourceEllipticPdeMethods() throw(Exception)
    {
        // Set up cell population
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a PDE object
        CellwiseSourceEllipticPde<2> pde(cell_population, 0.05);

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

    void TestCellwiseSourceEllipticPdeArchiving() throw(Exception)
    {
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
};

#endif /*TESTCELLBASEDELLIPTICPDES_HPP_*/
