/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef TESTCELLBETACATENINWRITER_HPP_
#define TESTCELLBETACATENINWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CryptCellsGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellBetaCateninWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCellBetaCateninWriter : public AbstractCellBasedTestSuite
{
public:

    void TestWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up a small MeshBasedCellPopulationWithGhostNodes using an appropriate cell-cycle model class
        CylindricalHoneycombMeshGenerator generator(5, 4, 1);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        double domain_length_for_wnt = 4.0*(sqrt(3.0)/2);
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisOne> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, false);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(domain_length_for_wnt);

        // Initialise the cell population to set up the ODE system associated with each cell-cycle model object
        cell_population.InitialiseCells();

        // Create an output directory for the writer
        std::string output_directory = "TestCellBetaCateninWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellBetaCateninWriter and test that the correct output is generated
        CellBetaCateninWriter<2,2> beta_catenin_writer(output_directory);
        beta_catenin_writer.OpenOutputFile();
        beta_catenin_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2,2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            beta_catenin_writer.VisitCell(*cell_iter, &cell_population);
        }
        beta_catenin_writer.CloseFile();

        /*
         * To verify this test, we can eyeball the results file to check that the correct default
         * initial conditions are output (see VanLeeuwen2009WntSwatCellCycleOdeSystem.hpp).
         */
        FileComparison(results_dir + "results.vizbetacatenin", "crypt/test/data/TestCellBetaCateninWriter/results.vizbetacatenin").CompareFiles();

        // Avoid memory leak
        WntConcentration<2>::Destroy();
    }
};

#endif /* TESTCELLBETACATENINWRITER_HPP_ */
