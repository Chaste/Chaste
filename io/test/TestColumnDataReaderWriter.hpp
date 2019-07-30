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


#ifndef _TESTCOLUMNDATAREADERWRITER_HPP_
#define _TESTCOLUMNDATAREADERWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"
#include "Exception.hpp"
#include <cassert>
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"


class TestColumnDataReaderWriter : public CxxTest::TestSuite
{
private:

    ColumnDataWriter* mpTestWriter;
    ColumnDataReader* mpTestReader;

    bool FilesMatch(std::string testfileName, std::string goodfileName)
    {
        bool matching = true;

        std::ifstream testfile(testfileName.c_str(), std::ios::in);
        std::ifstream goodfile(goodfileName.c_str(), std::ios::in);
        std::string teststring;
        std::string goodstring;

        if (!testfile.is_open() || !goodfile.is_open())
        {
            TS_FAIL("Files not present.");
        }

        while (getline(testfile, teststring))
        {
            getline(goodfile, goodstring);
            if (teststring != goodstring)
            {
                if (teststring.substr(0, 17) != "Created by Chaste")
                {
                    matching = false;
                }
            }
        }

        if (getline(goodfile, goodstring))
        {
            matching = false;
        }

        testfile.close();
        goodfile.close();
        return matching;
    }

    // Note: not using references so we can pass in temporary vectors
    void CompareVectors(std::vector<double> v1, std::vector<double> v2, double precision)
    {
        TS_ASSERT_EQUALS(v1.size(), v2.size());
        for (unsigned i=0; i<v1.size(); i++)
        {
            TS_ASSERT_DELTA(v1[i], v2[i], precision);
        }
    }

public:

    void TestCreateColumnWriter()
    {
        // Create a new writer
        std::string dirname("TestColumnDataReaderWriter");
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter(dirname, "test"));
        // Check that the output directory exists
        FileFinder dir(dirname, RelativeTo::ChasteTestOutput);
        TS_ASSERT(dir.Exists());
        TS_ASSERT(dir.IsDir());
        delete mpTestWriter;
    }

    void TestCreateColumnReader()
    {
        // File does not exist
        TS_ASSERT_THROWS_CONTAINS(mpTestReader = new ColumnDataReader("", "testdoesnotexist"), "Couldn't open info file: ");

        // File contains corrupt data
        TS_ASSERT_THROWS_THIS(mpTestReader = new ColumnDataReader("io/test/data", "testbad", false), "Couldn't read info file correctly");

        // .info file exists (unlimited) but _unlimited.dat file does not
        FileFinder dir("io/test/data", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_THIS(mpTestReader = new ColumnDataReader(dir, "UnlimitedMissing"), "Couldn't open ancillary data file");

        // .info file exists (fixed dim) but .dat file does not
        TS_ASSERT_THROWS_THIS(mpTestReader = new ColumnDataReader(dir, "DatMissing"), "Couldn't open data file");

        // Folder is missing
        FileFinder absent_dir("absent_dir", RelativeTo::ChasteTestOutput);
        TS_ASSERT_THROWS_CONTAINS(mpTestReader = new ColumnDataReader(absent_dir, "file"), "Directory does not exist: ");

        // Folder is defined with trailing '/' - check this runs OK.
        mpTestReader = new ColumnDataReader("io/test/data/", "testunlimitednegative2", false );
        delete mpTestReader;
    }

    void TestDetermineFieldWidth()
    {
        mpTestReader = new ColumnDataReader("io/test/data", "testfixed_good", false);
        TS_ASSERT_EQUALS(mpTestReader->GetFieldWidth(), 12u);

        delete mpTestReader;

        TS_ASSERT_THROWS_THIS(mpTestReader = new ColumnDataReader("io/test/data", "testenddefine_good", false),
                              "Unable to determine field width from file as cannot find any data entries");

        TS_ASSERT_THROWS_THIS(mpTestReader = new ColumnDataReader("io/test/data", "malformed_e_notation", false),
                              "Badly formatted scientific data field");
    }

    void TestDefineUnlimitedDimension()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testdefineunlimited", false));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"),
                "Unlimited dimension already set. Cannot be defined twice");

        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("Time", "m secs"),
                "Unlimited dimension already set. Cannot be defined twice");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("T,i,m,e", "msecs"),
                "Unlimited dimension already set. Cannot be defined twice");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("", "msecs"),
                "Unlimited dimension already set. Cannot be defined twice");

        delete mpTestWriter;
    }

    void TestDefineFixedDimension()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testdefinefixed", false));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension("Node","dimensionless", 5000));

        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension("Node ","dimensionless", 5000),
                "Variable name/units \'Node \' not allowed: may only contain alphanumeric characters or \'_\'.");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension("Node", "dimension.less", 5000),
                "Variable name/units \'dimension.less\' not allowed: may only contain alphanumeric characters or \'_\'.");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension("*Node*","dimensionless", 5000),
                "Variable name/units \'*Node*\' not allowed: may only contain alphanumeric characters or \'_\'.");

        delete mpTestWriter;
    }

    void TestDefineVariable()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testdefinevariable", false));
        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineVariable("Dummy", ""));

        // Bad variable names/units
        TS_ASSERT_THROWS_THIS(ik_var_id = mpTestWriter->DefineVariable("I_K", "milli amperes"), "Variable name/units \'milli amperes\' not allowed: may only contain alphanumeric characters or \'_\'.");
        TS_ASSERT_THROWS_THIS(ik_var_id = mpTestWriter->DefineVariable("I   K", "milliamperes"), "Variable name/units \'I   K\' not allowed: may only contain alphanumeric characters or \'_\'.");
        TS_ASSERT_THROWS_THIS(ik_var_id = mpTestWriter->DefineVariable("I.K", "milliamperes"), "Variable name/units \'I.K\' not allowed: may only contain alphanumeric characters or \'_\'.");
        TS_ASSERT_THROWS_THIS(ik_var_id = mpTestWriter->DefineVariable("", "milliamperes"), "Variable name not allowed: may not be blank.");

        TS_ASSERT_EQUALS(ina_var_id, 0);
        TS_ASSERT_EQUALS(ik_var_id, 1);

        delete mpTestWriter;
    }

    void TestEndDefineMode()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testenddefine", false));

        TS_ASSERT_THROWS_THIS(mpTestWriter->PutVariable(0, 0, 0), "Cannot put variables when in Define mode");

        // Ending define mode without having defined a dimension and a variable should raise an exception
        TS_ASSERT_THROWS_THIS(mpTestWriter->EndDefineMode(), "Cannot end define mode. No dimensions have been defined.");

        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_THIS(mpTestWriter->EndDefineMode(), "Cannot end define mode. No variables have been defined.");

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));

        TS_ASSERT_EQUALS(ina_var_id, 0);
        TS_ASSERT_EQUALS(ik_var_id, 1);

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineVariable("I_Ca", "milli amperes"), "Cannot define variables when not in Define mode");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"), "Unlimited dimension already set. Cannot be defined twice");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension("Node", "dimensionless", 5000), "Cannot define variables when not in Define mode");

        std::string output_dir = mpTestWriter->GetOutputDirectory();
        delete mpTestWriter;
#ifndef _MSC_VER
        TS_ASSERT(FilesMatch(output_dir + "testenddefine.dat", "io/test/data/testenddefine_good.dat"));
#else
        TS_ASSERT(FilesMatch(output_dir + "testenddefine.dat", "io/test/data/testenddefine_Windows_good.dat"));
#endif

        TS_ASSERT(FilesMatch(output_dir + "testenddefine.info", "io/test/data/testenddefine_good.info"));
    }

    void TestCantAddUnlimitedAfterEndDefine()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testdefine", false));
        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension("Node","dimensionless", 0), "Fixed dimension must be at least 1 long");
        mpTestWriter->DefineFixedDimension("Node","dimensionless", 5000);

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        TS_ASSERT_EQUALS(ina_var_id, 0);
        TS_ASSERT_EQUALS(ik_var_id,1);

        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("Time","msecs"), "Cannot define variables when not in Define mode");
        delete mpTestWriter;
    }

    void TestPutVariableInUnlimitedFile()
    {
        mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testunlimited", false);
        int time_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        for (int i=0; i<=10; i++)
        {
            mpTestWriter->PutVariable(time_var_id, (double)(i)/10 - 0.1);
            mpTestWriter->PutVariable(ina_var_id, 12.0);
            mpTestWriter->PutVariable(ica_var_id, ((double)((i+1)*(i+1)))*1e150); // NB: Last column
            mpTestWriter->PutVariable(ik_var_id, 7124.12355553*((double)(i+1))/12.0);
            mpTestWriter->AdvanceAlongUnlimitedDimension();
        }

        delete mpTestWriter;

        mpTestReader = new ColumnDataReader("TestColumnDataReaderWriter", "testunlimited");

        std::vector<double> values_ik = mpTestReader->GetValues("I_K");

        for (int i=0; i<11; i++)
        {
            TS_ASSERT_DELTA(values_ik[i]/(7124.12355553*((double)(i+1))/12.0), 1.0, 1e-3);
        }

        std::vector<double> time_values = mpTestReader->GetUnlimitedDimensionValues();
        for (int i=0; i < 10; i++)
        {
            TS_ASSERT_DELTA(time_values[i], i*0.1-0.1, 1e-3);
        }

        // Test for coverage
        TS_ASSERT_THROWS_THIS(values_ik = mpTestReader->GetValues("I_K", 3), "Data file has no fixed dimension");
        TS_ASSERT_THROWS_THIS(mpTestReader->GetValues("BadVar"), "'BadVar' is an unknown variable.");

        delete mpTestReader;

        // Compare with known good data using ColumnDataReader
        ColumnDataReader good_reader("io/test/data", "testunlimited", false);
        ColumnDataReader our_reader("TestColumnDataReaderWriter", "testunlimited");
        CompareVectors(our_reader.GetUnlimitedDimensionValues(), good_reader.GetUnlimitedDimensionValues(), 1e-6);
        CompareVectors(our_reader.GetValues("I_Na"), good_reader.GetValues("I_Na"), 1e-6);
        CompareVectors(our_reader.GetValues("I_K"), good_reader.GetValues("I_K"), 1e-6);
        CompareVectors(our_reader.GetValues("I_Ca"), good_reader.GetValues("I_Ca"), 1e-6);
    }

    void TestPutNegativeVariable()
    {
#ifdef  _MSC_VER
        _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testunlimitednegative", false));
        int time_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_THIS(time_var_id = mpTestWriter->DefineVariable("Time", "msecs"),
                              "Variable name: Time already in use as unlimited dimension");
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        int i = 12;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, -0.2));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) i));
        // Check very small values are OK now.
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, -1.1e-123));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, -3.3124e-123));

        mpTestWriter->AdvanceAlongUnlimitedDimension();
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.2));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) -i));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, 7124.12355553e99));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 3.124e123));

        // Check that an incorrect var id causes an exception:
        TS_ASSERT_THROWS_THIS(mpTestWriter->PutVariable(234, -33.124), "variableID unknown");

        std::string output_dir = mpTestWriter->GetOutputDirectory();
        delete mpTestWriter;

        TS_ASSERT(FilesMatch(output_dir + "testunlimitednegative.dat", "io/test/data/testunlimitednegative_good.dat"));

        // Compare with known good data using ColumnDataReader
        ColumnDataReader good_reader("io/test/data", "testunlimitednegative_good", false);
        ColumnDataReader our_reader("TestColumnDataReaderWriter", "testunlimitednegative");
        CompareVectors(our_reader.GetUnlimitedDimensionValues(), good_reader.GetUnlimitedDimensionValues(), 1e-6);
        CompareVectors(our_reader.GetValues("I_Na"), good_reader.GetValues("I_Na"), 1e-6);
        CompareVectors(our_reader.GetValues("I_K"), good_reader.GetValues("I_K"), 1e-6);
        CompareVectors(our_reader.GetValues("I_Ca"), good_reader.GetValues("I_Ca"), 1e-6);
    }

    void TestPutVariableInFixedFileAndPrecision()
    {
#ifdef  _MSC_VER
        _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testfixed", false, 3)); // precision = 3

        int node_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        int short_id = 0;

        TS_ASSERT_THROWS_NOTHING(node_var_id = mpTestWriter->DefineFixedDimension("Node", "dimensionless", 4));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_THIS(node_var_id = mpTestWriter->DefineVariable("Node", "dimensionless"),
                "Variable name: Node already in use as fixed dimension");
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(short_id = mpTestWriter->DefineVariable("Short_column", "dimensionless"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_THIS(mpTestWriter->PutVariable(node_var_id, 0, -1), "Dimension position not supplied");
        TS_ASSERT_THROWS_THIS(mpTestWriter->PutVariable(node_var_id, 0, -2), "Dimension position out of range");

        for (unsigned i=0; i<4; i++)
        {
            mpTestWriter->PutVariable(node_var_id, (double)(i+1), i);
            mpTestWriter->PutVariable(ina_var_id, 12.0, i);
            mpTestWriter->PutVariable(ica_var_id, ((double)((i+1)*(i+1)))/3.0, i);
            mpTestWriter->PutVariable(ik_var_id, 7124.12355553*((double)(i+1))/12.0, i);
        }

        for (unsigned i=0; i<2; i++)
        {
            mpTestWriter->PutVariable(short_id, (double)(i), i);
        }

        std::string output_dir = mpTestWriter->GetOutputDirectory();
        delete mpTestWriter;

        TS_ASSERT(FilesMatch(output_dir + "testfixed.dat", "io/test/data/testfixed_good.dat"));
        TS_ASSERT(FilesMatch(output_dir + "testfixed.info", "io/test/data/testfixed_good.info"));

        TS_ASSERT_THROWS_NOTHING(mpTestReader = new ColumnDataReader("TestColumnDataReaderWriter", "testfixed"));

        TS_ASSERT_THROWS_THIS(mpTestReader->GetValues("BadVar", 0), "Unknown variable");

        for (int i=0; i<4; i++)
        {
            std::vector<double> values_ik = mpTestReader->GetValues("I_K",  i);
            TS_ASSERT_DELTA(values_ik[0]/(7124.12355553*((double)(i+1))/12.0), 1.0, 1e-3);
        }

        for (int i=0; i<4; i++)
        {
            std::vector<double> values_short = mpTestReader->GetValues("Short_column", i);
            if (i<2)
            {
                TS_ASSERT_DELTA(values_short[0], (double) i, 1e-3);
            }
            else
            {
                //Missing data in short column
                TS_ASSERT_DELTA(values_short[0], DBL_MAX, 1e-3);
            }
        }

        TS_ASSERT_THROWS_THIS(mpTestReader->GetValues("non-existent_variable",1),
                "Unknown variable");

        // Check that get unlimited dimension values throws
        TS_ASSERT_THROWS_THIS(mpTestReader->GetUnlimitedDimensionValues(),
                "Data file has no unlimited dimension");

        delete mpTestReader;

        // Precision exceptions
        TS_ASSERT_THROWS_THIS(mpTestWriter = new ColumnDataWriter("","", false,1),
                "Precision must be between 2 and 20 (inclusive)");
        TS_ASSERT_THROWS_THIS(mpTestWriter = new ColumnDataWriter("","", false,21),
                "Precision must be between 2 and 20 (inclusive)");

        // Compare with known good data using ColumnDataReader
        ColumnDataReader good_reader("io/test/data", "testfixed_good", false);
        ColumnDataReader our_reader("TestColumnDataReaderWriter", "testfixed");
        for (int i=0; i<4; i++)
        {
            CompareVectors(our_reader.GetValues("I_Na", i), good_reader.GetValues("I_Na", i), 1e-6);
            CompareVectors(our_reader.GetValues("I_K", i), good_reader.GetValues("I_K", i), 1e-6);
            CompareVectors(our_reader.GetValues("I_Ca", i), good_reader.GetValues("I_Ca", i), 1e-6);
            CompareVectors(our_reader.GetValues("Short_column", i), good_reader.GetValues("Short_column", i), 1e-6);
        }
    }

    // This test is abstracted from TestIonicModels::TestSolverForLR91WithRegularStimulus where data from a recently run
    // cell model are checked against data from an ancient run of the same model
    void TestReadingAnOlderFile()
    {
        ColumnDataReader old_data_reader("io/test/data", "Lr91RegularStimValidData", false);
        std::vector<double> oldstyle_v_values = old_data_reader.GetValues("V");
        ColumnDataReader new_data_reader("io/test/data", "Lr91RegularStim", false);
        std::vector<double> newstyle_v_values = new_data_reader.GetValues("membrane_voltage");
        for (unsigned i=0; i<oldstyle_v_values.size(); i++)
        {
            TS_ASSERT_DELTA(oldstyle_v_values[i], newstyle_v_values[i], 1e-5);
        }
    }
    // This test is abstracted from TestPyCmlNightly
    void TestReadingAnotherOldFile()
    {
        ColumnDataReader old_data_reader("io/test/data", "earm_noble_model_1990", false);
        std::vector<double> oldstyle_v_values = old_data_reader.GetValues("membrane_voltage");
        for (unsigned i=10; i<oldstyle_v_values.size(); i++)
        {
            TS_ASSERT_DELTA(oldstyle_v_values[i], -75, 25);
        }
    }
    void TestPutNegativeVariableInFixedFile()
    {

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testfixed_negatives", false));
        int node_var_id = mpTestWriter->DefineFixedDimension("Node", "dimensionless", 4);
        int ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes");
        int ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes");
        int ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes");
        TS_ASSERT_THROWS_THIS(node_var_id = mpTestWriter->DefineVariable("Node", "dimensionless"), "Variable name: Node already in use as fixed dimension");

        mpTestWriter->EndDefineMode();
        int twelve = 12;

        mpTestWriter->PutVariable(ina_var_id, (double) twelve, 0);
        mpTestWriter->PutVariable(ina_var_id, (double) -twelve, 1);
        mpTestWriter->PutVariable(ica_var_id, -33.124,3);
        mpTestWriter->PutVariable(ik_var_id, 7124.12355553,3);
        mpTestWriter->AdvanceAlongUnlimitedDimension();
        TS_ASSERT_THROWS_THIS(mpTestWriter->PutVariable(ica_var_id, -63.124,2), "Cannot advance along unlimited dimension if it is not defined");

        // Note: the above call to PutVariable will, in effect, execute AdvanceAlongUnlimitedDimension and
        //       therefore throw an exception, hence we have to repeat that call to PutVariable below
        mpTestWriter->PutVariable(ica_var_id, -63.124,2);
        mpTestWriter->PutVariable(node_var_id, 1,0);
        mpTestWriter->PutVariable(node_var_id, -4,3);

        delete mpTestWriter;

        // This won't be true, as we use an old-format 'good' file, for coverage
        //TS_ASSERT(FilesMatch(output_dir + "testfixed_negatives.dat", "io/test/data/testfixed_negatives_good.dat"));

        // Compare with known good data using ColumnDataReader
        ColumnDataReader good_reader("io/test/data", "testfixed_negatives_good", false);
        ColumnDataReader our_reader("TestColumnDataReaderWriter", "testfixed_negatives");
        for (int i=0; i<4; i++)
        {
            CompareVectors(our_reader.GetValues("I_Na", i), good_reader.GetValues("I_Na", i), 1e-6);
            CompareVectors(our_reader.GetValues("I_K", i), good_reader.GetValues("I_K", i), 1e-6);
            CompareVectors(our_reader.GetValues("I_Ca", i), good_reader.GetValues("I_Ca", i), 1e-6);
        }
    }

    void TestPutVariableInFixedandUnlimitedFile()
    {
#ifdef _MSC_VER
        _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testfixedandunlimited", false));

        int time_var_id = 0;
        int node_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;

        TS_ASSERT_THROWS_NOTHING(node_var_id = mpTestWriter->DefineFixedDimension("Node", "dimensionless", 4));
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_THIS(time_var_id = mpTestWriter->DefineVariable("Time", "msecs"), "Variable name: Time already in use as unlimited dimension");
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.1));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(node_var_id, 0, 0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, 12.0, 0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, 12.0, 1));
        // Last column
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 1.1e130, 0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, -1.1e130, 1));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 1.1e-130, 2));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, -33.124, 3));
        // Penultimate column
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, -3.3e111, 0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, 3.3e-111, 1));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, -3.3e-111, 2));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, 7124.12355553, 3));

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->AdvanceAlongUnlimitedDimension());

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(node_var_id, 0, 0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 63.124,2));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, -35.124,3));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.2));
        TS_ASSERT_THROWS_THIS(mpTestWriter->PutVariable(time_var_id, 0.2,3), "Dimension position supplied, but not required");

        std::string output_dir = mpTestWriter->GetOutputDirectory();
        delete mpTestWriter;
        TS_ASSERT(FilesMatch(output_dir + "testfixedandunlimited_unlimited.dat",
                             "io/test/data/testfixedandunlimited_unlimited.dat"));

        TS_ASSERT(FilesMatch(output_dir + "testfixedandunlimited_000001.dat",
                             "io/test/data/testfixedandunlimited_000001.dat"));

        TS_ASSERT_THROWS_NOTHING(mpTestReader = new ColumnDataReader("TestColumnDataReaderWriter", "testfixedandunlimited"));

        std::vector<double> time_values = mpTestReader->GetUnlimitedDimensionValues();
        std::vector<double> ica_values = mpTestReader->GetValues("I_Ca", 3);
        for (int i=0; i<2; i++)
        {
            TS_ASSERT_DELTA(time_values[i],(i+1)*0.1,1e-3);
            TS_ASSERT_DELTA(ica_values[i],-33.124 - i * 2,1e-3);
        }

        // Check exception thrown if dimension is not given
        TS_ASSERT_THROWS_THIS(ica_values = mpTestReader->GetValues("I_Ca"), "Data file has fixed dimension which must be specified");

        delete mpTestReader;

        // Compare with known good data using ColumnDataReader
        ColumnDataReader good_reader("io/test/data", "testfixedandunlimited", false);
        ColumnDataReader our_reader("TestColumnDataReaderWriter", "testfixedandunlimited");
        for (int i=0; i<4; i++)
        {
            CompareVectors(our_reader.GetValues("I_Na", i), good_reader.GetValues("I_Na", i), 1e-6);
            CompareVectors(our_reader.GetValues("I_K", i), good_reader.GetValues("I_K", i), 1e-6);
            CompareVectors(our_reader.GetValues("I_Ca", i), good_reader.GetValues("I_Ca", i), 1e-6);
        }
        CompareVectors(our_reader.GetUnlimitedDimensionValues(), good_reader.GetUnlimitedDimensionValues(), 1e-6);
    }

    /*
     * This test is just to cover the line in ColumnDataWriter::PutVariable where
     * the fixed and unlimited dimensions are both set and the unlimited parameter
     * (i.e. time) is passed in negative.
     */
    void TestNegativeWithFixedAndUnlimitedDefined()
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testunlimitednegative2", false));

        int time_var_id = 0;
        int node_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_NOTHING(node_var_id = mpTestWriter->DefineFixedDimension("Node", "dimensionless", 4));
        TS_ASSERT_THROWS_NOTHING( ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, -0.2));

        // Make sure there's a data value to read the field width from
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 1.1e-130, 2));

        TS_ASSERT_EQUALS(node_var_id,-1);

        // Remember to delete - this closes the writer cleanly and means any data left
        // unwritten will be written to the datafile
        delete mpTestWriter;

        // Compare with known good data using ColumnDataReader
        ColumnDataReader good_reader("io/test/data", "testunlimitednegative2", false);
        ColumnDataReader our_reader("TestColumnDataReaderWriter", "testunlimitednegative2");
        for (int i=0; i<4; i++)
        {
            CompareVectors(our_reader.GetValues("I_Ca", i), good_reader.GetValues("I_Ca", i), 1e-6);
        }
        CompareVectors(our_reader.GetUnlimitedDimensionValues(), good_reader.GetUnlimitedDimensionValues(), 1e-6);
    }

    /** This test is to highlight portability issues with the column data reader.
     *  GNU/Linux and MacOSX write the same file when running a basic Luo-Rudy cell test, but MacOSX has
     *  issues in reading it back.
     *
     */
    void TestReadingFileWithVerySmallHighPrecisionNumbers()
    {
        ColumnDataReader reader("io/test/data", "lr91_chaste", false);

        //Column (index 0) has time
        std::vector<double> t = reader.GetValues("Time");
        //Column (index 2) is unaffected
        std::vector<double> m_gate = reader.GetValues("fast_sodium_current_m_gate__m");
        //Column (index 3) has 3-digit exponents - some less than 1e-300
        std::vector<double> h_gate = reader.GetValues("fast_sodium_current_h_gate__h");
        //Column (index 4) might be misread becuase of 3-digit exponents
        std::vector<double> j_gate = reader.GetValues("fast_sodium_current_j_gate__j");

        // The time row before the problem
        TS_ASSERT_DELTA(t[84], 84, 1e-15);
        TS_ASSERT_DELTA(m_gate[84], 9.97769687e-01, 1e-09);
        TS_ASSERT_DELTA(h_gate[84], 4.71914786e-98, 1e-106);
        TS_ASSERT_DELTA(j_gate[84], 6.56720224e-05, 1e-13);

        // This time row has a value below 10^-100
        TS_ASSERT_DELTA(t[85], 85, 1e-15);
        TS_ASSERT_DELTA(m_gate[85], 9.97792458e-01, 1e-09);
        TS_ASSERT_DELTA(h_gate[85], 5.87710052e-101, 1e-109);
        TS_ASSERT_DELTA(j_gate[85], 4.88974038e-05, 1e-13);

        // This time row has a value above DBL_MIN = 2.2250738585072014e-308
        TS_ASSERT_DELTA(t[157], 157, 1e-15);
        TS_ASSERT_DELTA(m_gate[157], 9.96673677e-01, 1e-09);
        TS_ASSERT_DELTA(h_gate[157], 3.56079375e-307, DBL_MIN);
        TS_ASSERT_DELTA(j_gate[157], 3.03889760e-14, 1e-23);

        // This time row has a value below DBL_MIN = 2.2250738585072014e-308
        // This means that it's not actually expressible in normalised double precision arithmetic
        // (This number mimicks what might happen in a Mac simulation)
        TS_ASSERT_DELTA(t[158], 158, 1e-15);
        TS_ASSERT_DELTA(m_gate[158], 9.96635012e-01, 1e-09);
        TS_ASSERT_DELTA(h_gate[158], 0.0, DBL_MIN);
        TS_ASSERT_DELTA(j_gate[158], 2.26824440e-14, 1e-23);

        // This number mimicks the lowest positive number produced in a Mac simulation
        // Number in file is 3.95252517e-323.  Note that gcc doesn't like positive numbers less than
        // "min positive subnormal number" = 4.9406564584124654e-324.
        TS_ASSERT_DELTA(h_gate[159], 0.0, DBL_MIN);

        // Number in file is 1.00000000e-330.  Note that this is not expressible in double precision arithmetic
        TS_ASSERT_DELTA(h_gate[160], 0.0, DBL_MIN);

        // Number in file is 1.00000000e-999.  Note that this is not expressible in double precision arithmetic
        TS_ASSERT_DELTA(h_gate[161], 0.0, DBL_MIN);
    }

    /**
     *
     *
     *  This test establishes that writing 3-digit exponents doesn't work properly (in Linux and Windows)
     *
     *  This test is also to highlight portability issues with the column data reader.
     *  MacOSX (clang) is able to write numbers which are smaller than DBL_MIN = 2.2250738585072014e-308
     *  but is unable to read them back without signalling iostream::fail.
     *
     *  Meanwhile, Gnu/Linux is able to read back all small numbers and round to zero if the number is inexpressible
     *
     */
    void TestWritingAndReadingWithThreeDigitExponents()
    {
        mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "widenumbers", false);

        int time_var_id = mpTestWriter->DefineUnlimitedDimension("Time", "aeons");
        int small_pos_id = mpTestWriter->DefineVariable("Small", "dimensionless");
        int large_pos_id = mpTestWriter->DefineVariable("Large", "dimensionless");
        int small_neg_id = mpTestWriter->DefineVariable("SmallNeg", "dimensionless");
        int large_neg_id = mpTestWriter->DefineVariable("LargeNeg", "dimensionless");
        mpTestWriter->EndDefineMode();

        double small = 100*DBL_MIN;
        // DBL_MIN = 2.2250738585072014e-308
        double large = 1e99;

        // Store what the numbers should be for comparison
        std::vector<double> small_values_check;
        std::vector<double> time_values_check;

        mpTestWriter->PutVariable(time_var_id, -1.0e123);
        time_values_check.push_back(-1.0e123);
        mpTestWriter->PutVariable(small_pos_id, 1.0);
        small_values_check.push_back(1.0);
        mpTestWriter->PutVariable(small_neg_id, 1.0);
        mpTestWriter->PutVariable(large_pos_id, 1.0);
        mpTestWriter->PutVariable(large_neg_id, 1.0);
        mpTestWriter->AdvanceAlongUnlimitedDimension();

        mpTestWriter->PutVariable(time_var_id, -0.0);
        time_values_check.push_back(-0.0);
        mpTestWriter->PutVariable(small_pos_id, -1.0);
        small_values_check.push_back(-1.0);
        mpTestWriter->PutVariable(small_neg_id, -1.0);
        mpTestWriter->PutVariable(large_pos_id, -1.0);
        mpTestWriter->PutVariable(large_neg_id, -1.0);
        mpTestWriter->AdvanceAlongUnlimitedDimension();

        for (unsigned i= 1; i<6; i++, small/=10.0, large*=10.0)
        {
            mpTestWriter->PutVariable(time_var_id, i);
            time_values_check.push_back(i);
            mpTestWriter->PutVariable(small_pos_id, small);
            small_values_check.push_back(small);
            mpTestWriter->PutVariable(small_neg_id, -small);
            mpTestWriter->PutVariable(large_pos_id, large);
            mpTestWriter->PutVariable(large_neg_id, -large);
            mpTestWriter->AdvanceAlongUnlimitedDimension();
        }
        delete mpTestWriter;

        ColumnDataReader reader("TestColumnDataReaderWriter", "widenumbers");
        /*
         * Check our small numbers are within DBL_MIN of 0. Sometimes they get
         * rounded to 0, sometimes not, depending on the compiler. Also check
         * our "time" numbers.
         */
        std::vector<double> small_values = reader.GetValues("Small");
        std::vector<double> time_values  = reader.GetValues("Time");
        for (unsigned i=0; i<6; i++)
        {
            TS_ASSERT_DELTA(small_values[i], small_values_check[i], DBL_MIN);
            TS_ASSERT_DELTA(time_values[i], time_values_check[i], 0.1);
        }
    }
};

#endif //_TESTCOLUMNDATAREADERWRITER_HPP_
