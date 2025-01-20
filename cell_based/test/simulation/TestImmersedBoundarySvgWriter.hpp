/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYSVGWRITER_HPP_
#define TESTIMMERSEDBOUNDARYSVGWRITER_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "DifferentiatedCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "SmartPointers.hpp"

// Includes from projects/ImmersedBoundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySvgWriter.hpp"
#include "FileComparison.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundarySvgWriter : public CxxTest::TestSuite
{
public:

    void TestContructor()
    {
      ImmersedBoundarySvgWriter<2> svg_writer;
      TS_ASSERT_EQUALS(svg_writer.GetSamplingMultiple(), 100u);
      TS_ASSERT_EQUALS(svg_writer.GetSvgSize(), 1600.0);
    }

    void TestGetSetMethods()
    {
      ImmersedBoundarySvgWriter<2> svg_writer;
      TS_ASSERT_EQUALS(svg_writer.GetSamplingMultiple(), 100u);
      TS_ASSERT_EQUALS(svg_writer.GetSvgSize(), 1600.0);

      svg_writer.SetSamplingMultiple(10u);
      TS_ASSERT_EQUALS(svg_writer.GetSamplingMultiple(), 10u);

      svg_writer.SetSvgSize(1200.0);
      TS_ASSERT_EQUALS(svg_writer.GetSvgSize(), 1200.0);
    }

    void TestOutputParametersWithImmersedBoundarySimulationModifier()
    {
        std::string output_directory = "TestOutputParametersWithImmersedBoundarySvgWriter";
        OutputFileHandler output_file_handler(output_directory, false);

        MAKE_PTR(ImmersedBoundarySvgWriter<2>, p_modifier);
        TS_ASSERT_EQUALS(p_modifier->GetIdentifier(), "ImmersedBoundarySvgWriter-2");

        out_stream modifier_parameter_file = output_file_handler.OpenOutputFile("ImmersedBoundarySvgWriter.parameters");
        p_modifier->OutputSimulationModifierParameters(modifier_parameter_file);
        modifier_parameter_file->close();

        // Compare the generated file in test output with a reference copy in the source code
        FileFinder generated = output_file_handler.FindFile("ImmersedBoundarySvgWriter.parameters");
        FileFinder reference("cell_based/test/data/TestImmersedBoundarySvgWriter/ImmersedBoundarySvgWriter.parameters",
                             RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated, reference);
        TS_ASSERT(comparer.CompareFiles());
    }

    void TestSetupSolve()
    {
        SimulationTime::Instance()->SetStartTime(0.0);

        // Create a single node, single element mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
        nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
        nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

        std::vector<ImmersedBoundaryElement<2, 2>*> elems;
        elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

        ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

        // Set up a cell population
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
        ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

        OutputFileHandler handler("TestImmersedBoundarySvgWriter", false);
        std::string svg_directory = handler.GetOutputDirectoryFullPath();

        ImmersedBoundarySvgWriter<2> modifier;
        modifier.SetupSolve(cell_population, svg_directory);

        TS_ASSERT_EQUALS(modifier.mOutputDirectory, svg_directory);
        TS_ASSERT_DIFFERS(modifier.mSvgHeader, "");
        TS_ASSERT_EQUALS(modifier.mSvgFooter, "</svg>\n");

        SimulationTime::Instance()->Destroy();
    }

    void TestWriting()
    {
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        // Create a single node, single element mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
        nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
        nodes.push_back(new Node<2>(2, true, 0.01, 0.2));
        nodes.push_back(new Node<2>(2, true, 0.99, 0.2));
        nodes.push_back(new Node<2>(2, true, 0.4, 0.01));
        nodes.push_back(new Node<2>(2, true, 0.7, 0.99));

        std::vector<ImmersedBoundaryElement<2, 2>*> elems;
        elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

        ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

        // Set up a cell population
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
        ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

        OutputFileHandler handler("TestImmersedBoundarySvgWriter", false);

        ImmersedBoundarySvgWriter<2> modifier;
        modifier.SetupSolve(cell_population, handler.GetRelativePath());
        modifier.SetSamplingMultiple(1);

        modifier.UpdateAtEndOfTimeStep(cell_population);

        FileComparison comparer(handler.GetOutputDirectoryFullPath() + "results_000000.svg", "cell_based/test/data/TestImmersedBoundarySvgWriter/results_000000.svg");
        comparer.CompareFiles();

    }

    void TestArchiving()
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ImmersedBoundarySvgWriter.arch";

        // Separate scope to write the archive
        {
            // Initialise a growth modifier and set a non-standard mature target area
            ImmersedBoundarySvgWriter<2> writer;
            writer.SetSamplingMultiple(10u);
            writer.SetSvgSize(1200.0);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize
            output_arch << writer;
        }

        // Separate scope to read the archive
        {
            ImmersedBoundarySvgWriter<2> writer;

            // Restore the modifier
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> writer;

            TS_ASSERT_EQUALS(writer.GetSamplingMultiple(), 10u);
            TS_ASSERT_EQUALS(writer.GetSvgSize(), 1200.0);

        }
    }

};

#endif /*TESTIMMERSEDBOUNDARYSVGWRITER_HPP_*/
