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

#ifndef TESTNUMERICFILECOMPARISON_HPP_
#define TESTNUMERICFILECOMPARISON_HPP_

#include <cxxtest/TestSuite.h>
#include "OutputFileHandler.hpp"
#include "NumericFileComparison.hpp"

class TestNumericFileComparison : public CxxTest::TestSuite
{
public:

    void TestBasicFunctionality() throw(Exception)
    {
        std::string base_file = "./global/test/data/random_data.txt";
        std::string noised_file = "./global/test/data/same_random_data_with_1e-4_noise.txt";

        NumericFileComparison same_data(base_file, base_file);
        TS_ASSERT(same_data.CompareFiles());

        NumericFileComparison different_data(base_file, noised_file);
        TS_ASSERT(different_data.CompareFiles(1e-4));
    }

    void TestIgnoreHeader() throw(Exception)
    {
        std::string boost_33_file = "./global/test/data/fake_archive_boost_1_33.txt";
        std::string boost_34_file = "./global/test/data/fake_archive_boost_1_34.txt";

        NumericFileComparison same_data(boost_33_file, boost_34_file);
        TS_ASSERT(same_data.CompareFiles(5.1e-6, 1));
    }

    void TestIgnoreProvenanceComment() throw(Exception)
    {
        std::string v1_file = "./global/test/data/fake_v_day_one.txt";
        std::string v2_file = "./global/test/data/fake_v_day_two.txt";

        NumericFileComparison same_data(v1_file, v2_file);
        //TS_ASSERT(same_data.CompareFiles(5e-4)); // Fails due to difference after comment
        TS_ASSERT(same_data.CompareFiles(5e-3)); // Difference after comment is below this tolerance
    }

    void TestRelativeDifference() throw(Exception)
    {
        std::string base_file = "./global/test/data/random_data.txt";
        std::string noised_file = "./global/test/data/same_random_data_with_1e-4_noise.txt";

        // Lower bound on data is 1e-2 so 1e-4 absolute noise is within 1e-2 relative tolerance
        NumericFileComparison different_data(base_file, noised_file);
        TS_ASSERT(different_data.CompareFiles(1e-2, 0, false));
    }
};

#endif /*TESTNUMERICFILECOMPARISON_HPP_*/
