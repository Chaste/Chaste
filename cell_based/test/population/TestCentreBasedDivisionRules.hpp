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

#ifndef TESTCENTREBASEDDIVISIONRULES_HPP_
#define TESTCENTREBASEDDIVISIONRULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"
#include "FixedCentreBasedDivisionRule.hpp"
#include "SmartPointers.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCentreBasedDivisionRules : public AbstractCellBasedTestSuite
{
public:

    void TestRandomDirectionCentreBasedDivisionRule()
    {
        /**
         * This tests the RandomDirectionCentreBasedDivisionRule. We first create a centre-based cell population and check whether we can
         * give the division rule to the population and get it back. Then we create 10000 division vectors and check that they point
         * uniformly in random directions.
         */

        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);

        // Set the division rule for our population to be the random direction division rule
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule_to_set(new RandomDirectionCentreBasedDivisionRule<2,2>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule = cell_population.GetCentreBasedDivisionRule();

        // Get 10000 division vectors, check each length, their mean and their variance
        c_vector<double, 2> average_axis = zero_vector<double>(2);
        c_vector<double, 2> axis_variance = zero_vector<double>(2);
        double average_angle = 0.0;
        double angle_variance = 0.0;
        for (unsigned iteration = 0; iteration < 10000; iteration++)
        {
            std::pair<c_vector<double, 2>, c_vector<double, 2> > positions = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);
            c_vector<double, 2> random_axis = positions.second - positions.first;

            // Each random vector should have norm equal to 0.5*0.3 = 0.15
            TS_ASSERT_DELTA(norm_2(random_axis), 0.15,1e-6);

            average_axis(0) += random_axis(0);
            axis_variance(0) += random_axis(0)*random_axis(0);
            average_axis(1) += random_axis(1);
            axis_variance(1) += random_axis(1)*random_axis(1);
            average_angle += asin(random_axis(0)/norm_2(random_axis));
            angle_variance += asin(random_axis(0)/norm_2(random_axis))*asin(random_axis(0)/norm_2(random_axis));
        }
        average_axis(0) /= 10000.0;
        average_axis(1) /= 10000.0;
        axis_variance(0) /= 10000.0;
        axis_variance(1) /= 10000.0;
        average_angle /= 10000.0;
        angle_variance /= 10000.0;

        TS_ASSERT_DELTA(average_axis(0), 0.0, 1e-2);
        TS_ASSERT_DELTA(average_axis(1), 0.0, 1e-2);

        // Each component of the axis variance should equal 0.5*(0.15^2) = 0.01125
        TS_ASSERT_DELTA(axis_variance(0), 0.01125, 1e-2);
        TS_ASSERT_DELTA(axis_variance(1), 0.01125, 1e-2);
        TS_ASSERT_DELTA(average_angle, 0.0, 1e-2);
        TS_ASSERT_DELTA(angle_variance, M_PI*M_PI/12.0, 1e-2);
    }

    void TestFixedCentreBasedDivisionRule()
    {
        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);
        c_vector<double, 2> expected_parent_location;
        expected_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        c_vector<double, 2> expected_daughter_location;
        expected_daughter_location[0] = 1.2;
        expected_daughter_location[1] = 3.4;

        // Set the division rule for our population to be the random direction division rule
        typedef FixedCentreBasedDivisionRule<2,2> FixedRule;
        MAKE_PTR_ARGS(FixedRule, p_division_rule_to_set, (expected_daughter_location));

        TS_ASSERT_DELTA(p_division_rule_to_set->rGetDaughterLocation()[0], 1.2, 1e-6);
        TS_ASSERT_DELTA(p_division_rule_to_set->rGetDaughterLocation()[1], 3.4, 1e-6);

        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule = cell_population.GetCentreBasedDivisionRule();

        // Check that the division rule returns the correct pair of vectors
        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 2> parent_location;
        parent_location = positions.first;
        TS_ASSERT_DELTA(parent_location[0], expected_parent_location[0], 1e-6);
        TS_ASSERT_DELTA(parent_location[1], expected_parent_location[1], 1e-6);

        c_vector<double, 2> daughter_location;
        daughter_location = positions.second;
        TS_ASSERT_DELTA(daughter_location[0], expected_daughter_location[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location[1], expected_daughter_location[1], 1e-6);
    }

    void TestArchiveRandomDirectionCentreBasedDivisionRule()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "RandomDirectionCentreBasedDivisionRule.arch";

        {
            boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new RandomDirectionCentreBasedDivisionRule<2,2>());

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) << p_division_rule;
        }

        {
            boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule;

            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> p_division_rule;

            typedef RandomDirectionCentreBasedDivisionRule<2,2> RandomRule;
            TS_ASSERT(dynamic_cast<RandomRule*>(p_division_rule.get()));
        }
    }

    void TestArchiveFixedCentreBasedDivisionRule()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "FixedCentreBasedDivisionRule.arch";

        {
            c_vector<double, 2> location;
            location[0] = -0.73;
            location[1] = 5.82;

            boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new FixedCentreBasedDivisionRule<2,2>(location));

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) << p_division_rule;
        }

        {
            boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule;

            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> p_division_rule;

            typedef FixedCentreBasedDivisionRule<2,2> FixedRule;
            TS_ASSERT(dynamic_cast<FixedRule*>(p_division_rule.get()));

            c_vector<double, 2> location;
            location = (dynamic_cast<FixedRule*>(p_division_rule.get()))->rGetDaughterLocation();
            TS_ASSERT_DELTA(location[0], -0.73, 1e-6);
            TS_ASSERT_DELTA(location[1], 5.82, 1e-6);
        }
    }
};

#endif /*TESTCENTREBASEDDIVISIONRULES_HPP_*/
