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


#ifndef _TESTCOLUMNDATAREADERWRITER_HPP_
#define _TESTCOLUMNDATAREADERWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"
#include "Exception.hpp"
#include <cassert>

using namespace std;

class TestColumnDataReaderWriter : public CxxTest::TestSuite
{
private:

    ColumnDataWriter* mpTestWriter;
    ColumnDataReader* mpTestReader;

    bool FilesMatch(std::string testfileName, std::string goodfileName)
    {
        bool matching = true;

        ifstream testfile(testfileName.c_str(),ios::in);
        ifstream goodfile(goodfileName.c_str(),ios::in);
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

    void TestCreateColumnWriter() throw(Exception)
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

    void TestCreateColumnReader() throw(Exception)
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
    }

    void TestDetermineFieldWidth() throw(Exception)
    {
        mpTestReader = new ColumnDataReader("io/test/data", "testfixed_good", false);
        TS_ASSERT_EQUALS(mpTestReader->GetFieldWidth(), 10u);

        delete mpTestReader;

        TS_ASSERT_THROWS_THIS(mpTestReader = new ColumnDataReader("io/test/data", "testenddefine_good", false),
                              "Unable to determine field width from file as cannot find any data entries");

        TS_ASSERT_THROWS_THIS(mpTestReader = new ColumnDataReader("io/test/data", "malformed_e_notation", false),
                              "Badly formatted scientific data field");
    }

    void TestDefineUnlimitedDimension() throw(Exception)
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

    void TestDefineFixedDimension() throw(Exception)
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

    void TestDefineVariable() throw(Exception)
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

    void TestEndDefineMode() throw(Exception)
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

        TS_ASSERT(FilesMatch(output_dir + "testenddefine.dat", "io/test/data/testenddefine_good.dat"));

        TS_ASSERT(FilesMatch(output_dir + "testenddefine.info", "io/test/data/testenddefine_good.info"));
    }

    void TestCantAddUnlimitedAfterEndDefine() throw(Exception)
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

    void TestPutVariableInUnlimitedFile() throw(Exception)
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

        std::string output_dir = mpTestWriter->GetOutputDirectory();
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

    void TestPutNegativeVariable() throw(Exception)
    {
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

    void TestPutVariableInFixedFileAndPrecision() throw(Exception)
    {
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

        TS_ASSERT_THROWS_THIS(std::vector<double> values_dodgy = mpTestReader->GetValues("non-existent_variable",1),
                "Unknown variable");

        // Check that get unlimited dimension values throws
        TS_ASSERT_THROWS_THIS(std::vector<double> unlimited_values = mpTestReader->GetUnlimitedDimensionValues(),
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

    void TestPutNegativeVariableInFixedFile() throw(Exception)
    {

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("TestColumnDataReaderWriter", "testfixed_negatives", false));
        int node_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(node_var_id = mpTestWriter->DefineFixedDimension("Node", "dimensionless", 4));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_THIS(node_var_id = mpTestWriter->DefineVariable("Node", "dimensionless"), "Variable name: Node already in use as fixed dimension");
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        int i = 12;

        mpTestWriter->PutVariable(ina_var_id, (double) i,0);
        mpTestWriter->PutVariable(ina_var_id, (double) -i,1);
        mpTestWriter->PutVariable(ica_var_id, -33.124,3);
        mpTestWriter->PutVariable(ik_var_id, 7124.12355553,3);
        mpTestWriter->AdvanceAlongUnlimitedDimension();
        TS_ASSERT_THROWS_THIS(mpTestWriter->PutVariable(ica_var_id, -63.124,2), "Cannot advance along unlimited dimension if it is not defined");

        // Note: the above call to PutVariable will, in effect, execute AdvanceAlongUnlimitedDimension and
        //       therefore throw an exception, hence we have to repeat that call to PutVariable below
        mpTestWriter->PutVariable(ica_var_id, -63.124,2);
        mpTestWriter->PutVariable(node_var_id, 1,0);
        mpTestWriter->PutVariable(node_var_id, -4,3);

        std::string output_dir = mpTestWriter->GetOutputDirectory();
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

    void TestPutVariableInFixedandUnlimitedFile() throw(Exception)
    {
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

        int i = 12;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.1));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(node_var_id, 0, 0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) i,0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) i,1));
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
    void TestNegativeWithFixedAndUnlimitedDefined() throw(Exception)
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
};

#endif //_TESTCOLUMNDATAREADERWRITER_HPP_
