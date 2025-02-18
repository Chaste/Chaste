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

#ifndef TESTIMMERSEDBOUNDARYDIVISIONRULES_HPP_
#define TESTIMMERSEDBOUNDARYDIVISIONRULES_HPP_

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ShortAxisImmersedBoundaryDivisionRule.hpp"
#include "SmartPointers.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestImmersedBoundaryDivisionRules : public AbstractCellBasedTestSuite
{
public:
    void TestShortAxisImmersedBoundaryDivisionRule()
    {
        // Create an immersed boundary mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(1, true, 0.2, 0.3));
        nodes.push_back(new Node<2>(2, true, 0.4, 0.4));
        nodes.push_back(new Node<2>(3, true, 0.3, 0.2));

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);

        std::vector<ImmersedBoundaryElement<2, 2>*> elements;
        elements.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes_elem_1));

        ImmersedBoundaryMesh<2, 2> mesh(nodes, elements);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumElements());

        // Create cell population
        ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractImmersedBoundaryDivisionRule<2> > p_division_rule = cell_population.GetImmersedBoundaryDivisionRule();
        c_vector<double, 2> short_axis = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);

        TS_ASSERT_DELTA(short_axis[0], 1.0 / sqrt(2.0), 1e-9);
        TS_ASSERT_DELTA(short_axis[1], -1.0 / sqrt(2.0), 1e-9);
    }

    void TestArchiveShortAxisImmersedBoundaryDivisionRule()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "ShortAxisImmersedBoundaryDivisionRule.arch";

        // Create data structures to store variables to test for equality here
        {
            boost::shared_ptr<AbstractImmersedBoundaryDivisionRule<2> > p_division_rule(new ShortAxisImmersedBoundaryDivisionRule<2>());

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Record values to test into data structures
            (*p_arch) << p_division_rule;
        }

        {
            boost::shared_ptr<AbstractImmersedBoundaryDivisionRule<2> > p_division_rule;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_division_rule;

            TS_ASSERT(dynamic_cast<ShortAxisImmersedBoundaryDivisionRule<2>*>(p_division_rule.get()));
        }
    }
};

#endif /*TESTIMMERSEDBOUNDARYDIVISIONRULES_HPP_*/
