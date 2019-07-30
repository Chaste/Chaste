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

#ifndef TESTSIMPLEDATAWRITER_HPP_
#define TESTSIMPLEDATAWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include "SimpleDataWriter.hpp"
#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

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
        FileComparison(results_dir + "std_vecs1.dat", "io/test/data/good_std_vec1.dat").CompareFiles();
        FileComparison(results_dir + "std_vecs2.dat", "io/test/data/good_std_vec2.dat").CompareFiles();
        FileComparison(results_dir + "std_vecs3.dat", "io/test/data/good_std_vec3.dat").CompareFiles();
    }
};

#endif /*TESTSIMPLEDATAWRITER_HPP_*/
