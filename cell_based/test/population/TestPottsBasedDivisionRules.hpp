/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTPOTTSBASEDDIVISIONRULES_HPP_
#define TESTPOTTSBASEDDIVISIONRULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "FixedPottsBasedDivisionRule.hpp"
#include "RandomDirectionPottsBasedDivisionRule.hpp"
#include "SmartPointers.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestPottsBasedDivisionRules : public AbstractCellBasedTestSuite
{
public:

///\todo #1737 - test the default division rule is correct

    void TestFixedPottsBasedDivisionRule()
    {
        // Create a simple 2D PottsMesh with two cells
        PottsMeshGenerator<2> generator(4, 2, 2, 2, 1, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);

        c_vector<double, 2> expected_vector;
        expected_vector(0) = 6.0/5.0;
        expected_vector(1) = 3.0/5.0;

        TS_ASSERT_THROWS_THIS(new FixedPottsBasedDivisionRule<2>(expected_vector),
            "Input argument must be a unit vector");

        expected_vector(0) = 4.0/5.0;

        // Set the division rule
        MAKE_PTR_ARGS(FixedPottsBasedDivisionRule<2>, p_division_rule_to_set, (expected_vector));
        cell_population.SetPottsBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractPottsBasedDivisionRule<2> > p_division_rule = cell_population.GetPottsBasedDivisionRule();

        c_vector<double, 2> division_vector = boost::static_pointer_cast<FixedPottsBasedDivisionRule<2> >(p_division_rule)->rGetDivisionVector();
        TS_ASSERT_DELTA(division_vector(0), 0.8, 1e-6);
        TS_ASSERT_DELTA(division_vector(1), 0.6, 1e-6);

        c_vector<double, 2> division_vector_again = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);
        TS_ASSERT_DELTA(division_vector_again(0), 0.8, 1e-6);
        TS_ASSERT_DELTA(division_vector_again(1), 0.6, 1e-6);
    }

    void TestRandomDirectionPottsBasedDivisionRule()
    {
        /**
         * This tests the RandomDirectionPottsBasedDivisionRule. We first create a Potts-based cell population and check whether we can
         * give the division rule to the population and get it back. Then we create 10000 division vectors and check that they point
         * uniformly in random directions.
         */

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);

        // Set the division rule for our population to be the random direction division rule
        boost::shared_ptr<AbstractPottsBasedDivisionRule<2> > p_division_rule_to_set(new RandomDirectionPottsBasedDivisionRule<2>());
        cell_population.SetPottsBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractPottsBasedDivisionRule<2> > p_division_rule = cell_population.GetPottsBasedDivisionRule();

        // Get 10000 division vectors, check each length, their mean and their variance.
        c_vector<double, 2> average_axis = zero_vector<double>(2);
        c_vector<double, 2> axis_variance = zero_vector<double>(2);
        double average_angle = 0.0;
        double angle_variance = 0.0;
        for (unsigned iteration = 0; iteration < 10000; iteration++)
        {
            c_vector<double, 2> random_axis = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);
            TS_ASSERT_DELTA(norm_2(random_axis), 1.0,1e-6);
            average_axis(0) += random_axis(0);
            axis_variance(0) += random_axis(0)*random_axis(0);
            average_axis(1) += random_axis(1);
            axis_variance(1) += random_axis(1)*random_axis(1);
            average_angle += asin(random_axis(0));
            angle_variance += asin(random_axis(0))*asin(random_axis(0));
        }
        average_axis(0) /= 10000.0;
        average_axis(1) /= 10000.0;
        axis_variance(0) /= 10000.0;
        axis_variance(1) /= 10000.0;
        average_angle /= 10000.0;
        angle_variance /= 10000.0;

        TS_ASSERT_DELTA(average_axis(0), 0.0, 1e-2);
        TS_ASSERT_DELTA(average_axis(1), 0.0, 1e-2);
        TS_ASSERT_DELTA(axis_variance(0), 0.5, 1e-2);
        TS_ASSERT_DELTA(axis_variance(1), 0.5, 1e-2);
        TS_ASSERT_DELTA(average_angle, 0.0, 1e-2);
        TS_ASSERT_DELTA(angle_variance, M_PI*M_PI/12.0, 1e-2);
    }


    void TestArchiveFixedPottsBasedDivisionRule() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "FixedPottsBasedDivisionRule.arch";

        // Create data structures to store variables to test for equality here
        {
            c_vector<double, 2> vector;
            vector(0) = 5.0/13.0;
            vector(1) = 12.0/13.0;
            boost::shared_ptr<AbstractPottsBasedDivisionRule<2> > p_division_rule(new FixedPottsBasedDivisionRule<2>(vector));

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Record values to test into data structures
            (*p_arch) << p_division_rule;
        }

        {
            boost::shared_ptr<AbstractPottsBasedDivisionRule<2> > p_division_rule;

            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_division_rule;

            TS_ASSERT(dynamic_cast<FixedPottsBasedDivisionRule<2>*>(p_division_rule.get()));

            c_vector<double, 2> location = (dynamic_cast<FixedPottsBasedDivisionRule<2>*>(p_division_rule.get()))->rGetDivisionVector();
            TS_ASSERT_DELTA(location[0], 5.0/13.0, 1e-6);
            TS_ASSERT_DELTA(location[1], 12.0/13.0, 1e-6);
        }
    }

    void TestArchiveRandomDirectionPottsBasedDivisionRule() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "RandomDirectionPottsBasedDivisionRule.arch";

        // Create data structures to store variables to test for equality here
        {
            boost::shared_ptr<AbstractPottsBasedDivisionRule<2> > p_division_rule(new RandomDirectionPottsBasedDivisionRule<2>());

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Record values to test into data structures
            (*p_arch) << p_division_rule;
        }

        {
            boost::shared_ptr<AbstractPottsBasedDivisionRule<2> > p_division_rule;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_division_rule;

            TS_ASSERT(dynamic_cast<RandomDirectionPottsBasedDivisionRule<2>*>(p_division_rule.get()));
        }
    }

};

#endif /*TESTPOTTSBASEDDIVISIONRULES_HPP_*/
