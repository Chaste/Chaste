/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTSIMPLEDATAWRITER_HPP_
#define TESTSIMPLEDATAWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include "SimpleDataWriter.hpp"
#include "OutputFileHandler.hpp"

class TestSimpleDataWriter : public CxxTest::TestSuite
{
public:

    void TestExceptions()
    {
        std::vector<std::vector<double> > empty_data;

        // No data
        TS_ASSERT_THROWS_THIS(SimpleDataWriter bad_writer("SimpleDataWriter", "bad1", empty_data),
                "Data vector is empty");

        std::vector<double> t;
        std::vector<double> x;

        t.push_back(0);
        x.push_back(0.1);
        x.push_back(0.2);

        // t and x have different sizes
        TS_ASSERT_THROWS_THIS(SimpleDataWriter bad_writer("SimpleDataWriter", "bad2", t,x),
                "Data vector sizes are not all equal");
    }

    void TestSimpleDataWriterWithStdVecs()
    {
        std::vector<double> t;
        std::vector<double> x;
        std::vector<double> y;

        for (unsigned i=0; i<4; i++)
        {
            t.push_back(i);
            x.push_back(2*i);
            y.push_back(i+10);
        }

        SimpleDataWriter writer1("SimpleDataWriter", "std_vecs1.dat", t, x);

        std::vector<std::vector<double> > full_data;
        full_data.push_back(t);
        full_data.push_back(x);
        full_data.push_back(y);

        SimpleDataWriter writer2("SimpleDataWriter", "std_vecs2.dat", full_data, false);

        SimpleDataWriter writer3("SimpleDataWriter", "std_vecs3.dat", t, false);

        // Do the testing now so that to also check the directory wasn't cleaned in the second write
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "SimpleDataWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "std_vecs1.dat  io/test/data/good_std_vec1.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "std_vecs2.dat  io/test/data/good_std_vec2.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "std_vecs3.dat  io/test/data/good_std_vec3.dat").c_str()), 0);
    }
};

#endif /*TESTSIMPLEDATAWRITER_HPP_*/
