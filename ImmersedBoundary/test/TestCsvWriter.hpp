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

#include <cxxtest/TestSuite.h>
#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"

// Includes from projects/ImmersedBoundary
#include "CsvWriter.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestCsvWriter : public CxxTest::TestSuite
{
public:

    void TestDirectoryNameAndFileName() throw(Exception)
    {
        CsvWriter writer;
        TS_ASSERT_EQUALS(writer.GetDirectoryName(), "");
        TS_ASSERT_EQUALS(writer.GetFileName(), "");

        writer.SetDirectoryName("TestCsvWriter");
        writer.SetFileName("file");

        TS_ASSERT_EQUALS(writer.GetDirectoryName(), "TestCsvWriter");
        TS_ASSERT_EQUALS(writer.GetFileName(), "file.csv");
    }

    void TestOutput() throw(Exception)
    {
        CsvWriter writer;
        writer.SetDirectoryName("TestCsvWriter");
        writer.SetFileName("file.csv");

        std::vector<std::string> headers;
        headers.push_back("UnsignedHeader1");
        headers.push_back("UnsignedHeader2");
        headers.push_back("DoubleHeader1");
        headers.push_back("DoubleHeader2");
        headers.push_back("StringHeader1");
        headers.push_back("StringHeader2");
        writer.AddHeaders(headers);

        std::vector<double> double_data1;
        double_data1.push_back(1.1);
        double_data1.push_back(1.2);
        double_data1.push_back(1.3);
        writer.AddData(double_data1);

        std::vector<std::string> string_data1;
        string_data1.push_back("string11");
        string_data1.push_back("string12");
        string_data1.push_back("string13");
        writer.AddData(string_data1);

        std::vector<unsigned> unsigned_data1;
        unsigned_data1.push_back(11);
        unsigned_data1.push_back(12);
        unsigned_data1.push_back(13);
        writer.AddData(unsigned_data1);

        std::vector<std::string> string_data2;
        string_data2.push_back("string21");
        string_data2.push_back("string22");
        string_data2.push_back("string23");
        writer.AddData(string_data2);

        std::vector<unsigned> unsigned_data2;
        unsigned_data2.push_back(11);
        unsigned_data2.push_back(12);
        unsigned_data2.push_back(13);
        writer.AddData(unsigned_data2);

        std::vector<double> double_data2;
        double_data2.push_back(2.1);
        double_data2.push_back(2.2);
        double_data2.push_back(2.3);
        writer.AddData(double_data2);

        writer.WriteDataToFile();

        // Compare output with saved file of what it should look like
        OutputFileHandler output_file_handler("TestCsvWriter", false);
        std::string results_directory = output_file_handler.GetOutputDirectoryFullPath();

        FileComparison(results_directory + "file.csv",
                "projects/ImmersedBoundary/test/data/TestCsvWriter/file.csv").CompareFiles();
    }

    void TestExceptions() throw(Exception)
    {
        CsvWriter writer;

        std::vector<std::string> empty_headers;
        TS_ASSERT_THROWS_THIS(writer.AddHeaders(empty_headers),
                "No strings provided in vector of headers");

        std::vector<double> empty_data;
        TS_ASSERT_THROWS_THIS(writer.AddData(empty_data),
                "Vectors passed to this writer must contain at least one element");

        TS_ASSERT_THROWS_THIS(writer.WriteDataToFile(),
                "There is no data to write");

        std::vector<double> data;
        data.push_back(7.5);
        data.push_back(8.2);
        writer.AddData(data);

        std::vector<double> wrong_size_data;
        wrong_size_data.push_back(4.3);

        TS_ASSERT_THROWS_THIS(writer.AddData(wrong_size_data),
                "All data vectors added must be the same size");

        TS_ASSERT_THROWS_THIS(writer.WriteDataToFile(),
                "Output directory has not been specified");

        writer.SetDirectoryName("TestCsvWriter");

        TS_ASSERT_THROWS_THIS(writer.WriteDataToFile(),
                "File name has not been specified");

        TS_ASSERT_THROWS_THIS(writer.SetFileName("file.wrongformat"),
                "This class is designed to write a file with extension .csv");

        writer.SetFileName("file");

        std::vector<std::string> wrong_size_headers;
        wrong_size_headers.push_back("Header1");
        wrong_size_headers.push_back("Header2");
        wrong_size_headers.push_back("Header3");
        writer.AddHeaders(wrong_size_headers);

        TS_ASSERT_THROWS_THIS(writer.WriteDataToFile(),
                "Expecting to write a header row, but header length does not match number of data rows");
    }
};
