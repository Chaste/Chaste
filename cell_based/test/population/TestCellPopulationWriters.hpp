/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef TESTCELLPOPULATIONWRITERS_HPP_
#define TESTCELLPOPULATIONWRITERS_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"

#include "NodeLocationWriter.hpp"
#include "NodeBasedCellPopulation.hpp"

#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCellPopulationWriters : public AbstractCellBasedTestSuite
{
public:

    void TestNodeLocationWriter() throw(Exception)
    {
    	EXIT_IF_PARALLEL;

    	// Set up a node based cell population.
    	std::vector<Node<3>* > nodes;
    	nodes.push_back(new Node<3>(0, false));
    	nodes.push_back(new Node<3>(1, false, 1.0, 1.0, 1.0));

    	NodesOnlyMesh<3> mesh;
    	mesh.ConstructNodesWithoutMesh(nodes, 1.5);

    	std::vector<CellPtr> cells;
    	CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> generator;
    	generator.GenerateBasic(cells, mesh.GetNumNodes());

    	NodeBasedCellPopulation<3> cell_population(mesh, cells);

    	std::string output_directory = "TestWriteNodeLocations";
		OutputFileHandler output_file_handler(output_directory, false);

		std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

		// Create a node location writer and write the files.
    	NodeLocationWriter<3,3> location_writer(output_directory);

    	location_writer.OpenOutputFile();

    	location_writer.Visit(&cell_population);

    	location_writer.CloseFile();

    	FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestNodeLocationWriter/results.viznodes").CompareFiles();

    	// Make sure we can append to files.
    	location_writer.OpenOutputFileForAppend();

    	location_writer.Visit(&cell_population);

    	location_writer.CloseFile();

    	FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestNodeLocationWriter/results.viznodes_twice").CompareFiles();

    	// Tidy up.
    	delete nodes[0];
    	delete nodes[1];
    }
};

#endif /*TESTCELLPOPULATIONWRITERS_HPP_*/
