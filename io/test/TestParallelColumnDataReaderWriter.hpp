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

#ifndef _TESTPARALLELCOLUMNDATAREADERWRITER_HPP_
#define _TESTPARALLELCOLUMNDATAREADERWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include "ColumnDataWriter.hpp"
#include "ParallelColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"
#include "Exception.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include <petsc.h>
#include "FileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestParallelColumnDataReaderWriter : public CxxTest::TestSuite
{
private:

    ParallelColumnDataWriter* mpParallelWriter;
    ColumnDataReader* mpReader;

    static const int num_nodes = 10;

public:

    void TestParallelColumnWriter()
    {
        int time_var_id=-1, var1_id=-1, var2_id=-1;

#ifdef _MSC_VER
        _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

        // Make a parallel data writer
        mpParallelWriter = new ParallelColumnDataWriter("TestParallelColumnDataWriter", "ParallelColumnWriter");

        TS_ASSERT_THROWS_NOTHING(time_var_id = mpParallelWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_NOTHING(mpParallelWriter->DefineFixedDimension("Node", "dimensionless", num_nodes));

        TS_ASSERT_THROWS_NOTHING(var1_id = mpParallelWriter->DefineVariable("Var1", "LightYears"));
        TS_ASSERT_THROWS_NOTHING(var2_id = mpParallelWriter->DefineVariable("Var2", "Angstroms"));

        TS_ASSERT_THROWS_NOTHING(mpParallelWriter->EndDefineMode());

        PetscTools::Barrier("TestParallelColumnWriter1");

        std::string output_dir = mpParallelWriter->GetOutputDirectory();

        // Test that the output directory and the .info file was created
        FileFinder info_file(output_dir+"ParallelColumnWriter.info", RelativeTo::Absolute);

        TS_ASSERT(info_file.Exists());

        // Set up some data in PETSc vectors
        Vec var1=PetscTools::CreateVec(num_nodes);
        Vec var2=PetscTools::CreateVec(num_nodes);
        Vec var3=PetscTools::CreateVec(num_nodes+1);

        double* var1_array,* var2_array;
        VecGetArray(var1, &var1_array);
        VecGetArray(var2, &var2_array);
        int lo, hi;
        VecGetOwnershipRange(var1,&lo,&hi);

        for (int global_index=lo; global_index<hi; global_index++)
        {
            var1_array[global_index-lo] = global_index;
            var2_array[global_index-lo] = -global_index * 1e100;
        }

        VecRestoreArray(var1, &var1_array);
        VecAssemblyBegin(var1);
        VecAssemblyEnd(var1);
        VecRestoreArray(var2, &var2_array);
        VecAssemblyBegin(var2);
        VecAssemblyEnd(var2);

        // Write out the data (Conventional)
        /*VecGetArray(var1, &var1_array);
        VecGetArray(var2, &var2_array);
        mpWriter->PutVariable(time_var_id, 0);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            mpWriter->PutVariable(var1_id, var1_array[global_index - lo], global_index);
            mpWriter->PutVariable(var2_id, var2_array[global_index - lo], global_index);
        }
        mpWriter->AdvanceAlongUnlimitedDimension();
        VecRestoreArray(var1, &var1_array);
        VecRestoreArray(var2, &var2_array);
        */

        // Write out the data (Parallel)
        mpParallelWriter->PutVariable(time_var_id, 0.1);
        TS_ASSERT_THROWS_NOTHING(mpParallelWriter->PutVector(var1_id, var1));
        TS_ASSERT_THROWS_NOTHING(mpParallelWriter->PutVector(var2_id, var2));

        // Throws since var3 is the wrong size
        TS_ASSERT_THROWS_THIS(mpParallelWriter->PutVector(var1_id, var3), "Size of vector does not match FixedDimensionSize.");

        // No-op if not master, writes anyway if we are
        TS_ASSERT_THROWS_NOTHING(mpParallelWriter->PutVariable(var1_id, 0.0, 0));

        TS_ASSERT_THROWS_NOTHING(mpParallelWriter->AdvanceAlongUnlimitedDimension());

        // Change the data
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
        VecSqrtAbs(var1);
#else
        VecSqrt(var1);
#endif
        VecAbs(var2);
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
        VecSqrtAbs(var2);
#else
        VecSqrt(var2);
#endif

        // Write out the data again (Conventional)
        /*
        VecGetArray(var1, &var1_array);
        VecGetArray(var2, &var2_array);
        mpWriter->PutVariable(time_var_id, 1);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            mpWriter->PutVariable(var1_id, var1_array[global_index - lo], global_index);
            mpWriter->PutVariable(var2_id, var2_array[global_index - lo], global_index);
        }
        */

        // Write out the data (Parallel)
        mpParallelWriter->PutVariable(time_var_id, 0.2);
        mpParallelWriter->PutVector(var1_id, var1);
        mpParallelWriter->PutVector(var2_id, var2);
        //mpParallelWriter->AdvanceAlongUnlimitedDimension();

        //delete mpParallelWriter;
        delete mpParallelWriter;

        PetscTools::Barrier("TestParallelColumnWriter2");

        std::vector<std::string> files_to_compare;
        files_to_compare.push_back("ColumnWriter.info");
        files_to_compare.push_back("ColumnWriter_000000.dat");
        files_to_compare.push_back("ColumnWriter_000001.dat");
        files_to_compare.push_back("ColumnWriter_unlimited.dat");

        for (unsigned i=0; i<files_to_compare.size(); i++)
        {
            FileComparison comparer(output_dir+ "Parallel" + files_to_compare[i],
                                    "io/test/data/" + files_to_compare[i]);
            TS_ASSERT(comparer.CompareFiles());
        }

        PetscTools::Destroy(var1);
        PetscTools::Destroy(var2);
        PetscTools::Destroy(var3);
    }

    void TestPutSlice()
    {
        // Create a vector slice

        const unsigned problem_size = 10;
        DistributedVectorFactory factory(problem_size);
        Vec striped = factory.CreateVec(2);

        DistributedVector distributed_vector = factory.CreateDistributedVector(striped);
        DistributedVector::Stripe zeros(distributed_vector, 0);
        DistributedVector::Stripe ones(distributed_vector, 1);

        // Write some values
        for (DistributedVector::Iterator index = distributed_vector.Begin();
             index!= distributed_vector.End();
             ++index)
        {
            zeros[index] = 0;
            ones[index] = 1;
        }

        // Write to file with parallel data writer
#ifdef _MSC_VER
        _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

        ParallelColumnDataWriter* p_parallel_writer = new ParallelColumnDataWriter("TestParallelColumnDataWriterStripe","Stripe");
        unsigned time_var_id = p_parallel_writer->DefineUnlimitedDimension("Time","msecs");
        unsigned var1_id = p_parallel_writer->DefineVariable("Var1","LightYears");
        p_parallel_writer->DefineFixedDimension("Node","dimensionless", problem_size);
        p_parallel_writer->EndDefineMode();
        PetscTools::Barrier("TestPutSlice1");
        std::string output_dir = p_parallel_writer->GetOutputDirectory();

        p_parallel_writer->PutVariable(time_var_id, 0.1);
        p_parallel_writer->PutVectorStripe(var1_id, ones);
        p_parallel_writer->AdvanceAlongUnlimitedDimension();

        // Check file
        PetscTools::Barrier("TestPutSlice2");

        std::vector<std::string> files_to_compare;
        files_to_compare.push_back("Stripe.info");
        files_to_compare.push_back("Stripe_000000.dat");
        files_to_compare.push_back("Stripe_unlimited.dat");

        p_parallel_writer->Close();//ensure files flushed before comparison

        for (unsigned i=0; i<files_to_compare.size(); i++)
        {
            FileComparison comparer(output_dir+files_to_compare[i], "io/test/data/" + files_to_compare[i]);
            TS_ASSERT(comparer.CompareFiles());
        }

        PetscTools::Destroy(striped);
        delete p_parallel_writer;
    }

    // Read back the data written in the test TestParallelColumnWriter above
    void TestColumnReader()
    {
        /*
         * There is no *Parallel* ColumnDataReader.  Since everyone might
         * need to know everything there's no point in only one processor
         * opening the file
         */

        // Make a parallel data writer
        mpReader = new ColumnDataReader("TestParallelColumnDataWriter", "ParallelColumnWriter");

        // Check that there's the correct number of files
        std::vector<double> time_stamps;
        time_stamps = mpReader->GetUnlimitedDimensionValues();
        TS_ASSERT_EQUALS(time_stamps[0],0.1);
        TS_ASSERT_EQUALS(time_stamps[1],0.2);
        TS_ASSERT_EQUALS(time_stamps.size(),2u);

        // Check that some of the data is correct
        std::vector<double> var1_node4;
        var1_node4 = mpReader->GetValues("Var1",4);
        TS_ASSERT_EQUALS(var1_node4[0], 4.0); // First time step
        TS_ASSERT_EQUALS(var1_node4[1], 2.0); // Second time step
        std::vector<double> var2_node4;
        var2_node4 = mpReader->GetValues("Var2",4);
        TS_ASSERT_EQUALS(var2_node4[0], -4 * 1e100); // First time step
        TS_ASSERT_DELTA(var2_node4[1], sqrt(4 * 1e100), 1e-4); // Second time step

        TS_ASSERT_THROWS_THIS(mpReader->GetValues("LifeSigns",4), "Unknown variable");
        TS_ASSERT_THROWS(mpReader->GetValues("Var1",10), const std::out_of_range &);

        // Delete the reader: makes sure that files are closed
        delete mpReader;
    }
};

#endif //_TESTPARALLELCOLUMNDATAREADERWRITER_HPP_
