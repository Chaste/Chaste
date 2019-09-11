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

#ifndef TESTHDF5DATAWRITER_HPP_
#define TESTHDF5DATAWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <cstring> // For strcpy

#include "ChasteSyscalls.hpp"
#include "DistributedVectorFactory.hpp"
#include "Hdf5DataReader.hpp"
#include "Hdf5DataWriter.hpp"
#include "MathsCustomFunctions.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "Warnings.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CompareHdf5ResultsFiles.hpp"

class TestHdf5DataWriter : public CxxTest::TestSuite
{
private:
    Hdf5DataWriter* mpTestWriter;

public:
    void TestSimpleParallelWriteDirectlyWithHdf5()
    {
        // File to write
        OutputFileHandler oh("TestHdf5DataWriter");
        std::string results_dir = oh.GetOutputDirectoryFullPath();
        std::string file_name = results_dir + "test.h5";

        // Set up a property list saying how we'll open the file
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL);

        // Create a file (collectively) and free the property list
        hid_t file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);

        // Initialise the data for this process
        const unsigned DIMS = 2;
        const unsigned X = 2;
        const unsigned Y = 5;
        int data[X][Y];
        for (unsigned i = 0; i < X; i++)
        {
            for (unsigned j = 0; j < Y; j++)
            {
                data[i][j] = 100 * PetscTools::GetMyRank() + 10 * i + j;
            }
        }

        // Create the dataspace for the dataset.
        hsize_t dimsf[DIMS]; // dataset dimensions
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        dimsf[0] = X * num_procs;
        dimsf[1] = Y;
        hid_t filespace = H5Screate_simple(DIMS, dimsf, NULL);

        // Create the dataset with default properties and close filespace.
        hid_t dset_id = H5Dcreate(file_id, "IntArray", H5T_NATIVE_INT, filespace,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        // Define a dataset in memory for this process
        hsize_t count[DIMS] = { X, Y };
        hid_t memspace = H5Screate_simple(DIMS, count, NULL);

        // Select hyperslab in the file.
        hsize_t offset[DIMS] = { PetscTools::GetMyRank() * X, 0 };
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        // Create property list for collective dataset write, and write!  Finally.
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        herr_t status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
        TS_ASSERT_EQUALS(status, 0);

        // Create dataspace for the name, unit attribute
        hsize_t columns[2] = { Y, 21 };
        hid_t colspace = H5Screate_simple(1, columns, NULL);

        // Create attribute
        char col_data[5][21];
        strcpy(col_data[0], "Noughth");
        strcpy(col_data[1], "First");
        strcpy(col_data[2], "Second");
        strcpy(col_data[3], "Third");
        strcpy(col_data[4], "Fourth");

        // Create the type 'char'
        hid_t char_type = H5Tcopy(H5T_C_S1);
        // H5Tset_strpad(char_type, H5T_STR_NULLPAD);
        H5Tset_size(char_type, 21);
        hid_t attr = H5Acreate(dset_id, "Name", char_type, colspace,
                               H5P_DEFAULT, H5P_DEFAULT);

        // Write to the attribute
        status = H5Awrite(attr, char_type, col_data);

        // Close dataspace & attribute
        H5Sclose(colspace);
        H5Aclose(attr);

        // Release resources and close the file
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Pclose(plist_id);
        H5Fclose(file_id);
    }

    static const unsigned data_size = 17;
    void TestPetscWriteDirectlyWithHdf5()
    {
        // Initialise a PETSc vector
        Vec a_vec = PetscTools::CreateVec(data_size);
        double* p_a_vec;
        VecGetArray(a_vec, &p_a_vec);
        int lo, hi;
        VecGetOwnershipRange(a_vec, &lo, &hi);
        for (int global_index = lo; global_index < hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            p_a_vec[local_index] = global_index + 100 * PetscTools::GetMyRank();
        }
        VecRestoreArray(a_vec, &p_a_vec);
        VecAssemblyBegin(a_vec);
        VecAssemblyEnd(a_vec);

        //VecView(a_vec, PETSC_VIEWER_STDOUT_WORLD);

        // File to write
        OutputFileHandler oh("TestHdf5DataWriter", false);
        std::string results_dir = oh.GetOutputDirectoryFullPath();
        std::string file_name = results_dir + "vec.h5";

        // Set up a property list saying how we'll open the file
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL);

        // Create a file (collectively) and free the property list
        hid_t file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);

        const unsigned DIMS = 1;

        // Create the dataspace for the dataset.
        //TS_ASSERT_EQUALS(data_size, hi-lo);
        hsize_t dimsf[DIMS] = { data_size }; // dataset dimensions

        hid_t filespace = H5Screate_simple(DIMS, dimsf, NULL);

        // Create the dataset with default properties and close filespace.
        hid_t dset_id = H5Dcreate(file_id, "TheVector", H5T_NATIVE_DOUBLE, filespace,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        // Define a dataset in memory for this process
        hsize_t count[DIMS] = { (unsigned)(hi - lo) };
        hid_t memspace = H5Screate_simple(DIMS, count, NULL);

        // Select hyperslab in the file.
        hsize_t offset[DIMS] = { (unsigned)(lo) };
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        // Create property list for collective dataset write, and write!  Finally.
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        VecGetArray(a_vec, &p_a_vec);
        herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, p_a_vec);
        VecRestoreArray(a_vec, &p_a_vec);

        TS_ASSERT_EQUALS(status, 0);

        // Release resources and close the file
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Pclose(plist_id);
        H5Fclose(file_id);

        PetscTools::Destroy(a_vec);
    }

    void TestReadAndChecksumHdf5()
    {
        // File to read
        OutputFileHandler oh("TestHdf5DataWriter", false);
        double data[data_size];
        std::string results_dir = oh.GetOutputDirectoryFullPath();
        std::string file_name = results_dir + "vec.h5";

        hsize_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        //H5Tget_nmembers(file_id);
        hsize_t dataset_id = H5Dopen(file_id, "TheVector", H5P_DEFAULT);
        hsize_t dxpl = H5Pcreate(H5P_DATASET_XFER);
        hsize_t edc = H5Pget_edc_check(dxpl);
        TS_ASSERT_EQUALS(edc, (hsize_t)1) //Checksum is enabled

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl, data);

        TS_ASSERT_EQUALS(status, 0);

        // Check the index
        for (unsigned i = 0; i < data_size; i++)
        {
            TS_ASSERT_EQUALS(((unsigned)data[i] % 100), i);
        }

        // Check the final component
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        // The last component was owned by processor "num_procs-1"
        TS_ASSERT_EQUALS(((int)data[data_size - 1] / 100), num_procs - 1);

        H5Pclose(dxpl);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
    }

    void TestHdf5DataWriterMultipleColumns()
    {
        int number_nodes = 100;

        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_test_multi_column", false);

        // Coverage
        TS_ASSERT_EQUALS(writer.GetUsingCache(), false);

        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        int ik_id = writer.DefineVariable("I_K", "milliamperes");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");

        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2 = factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        // Write some values
        for (DistributedVector::Iterator index = distributed_vector_1.Begin();
             index != distributed_vector_1.End();
             ++index)
        {
            distributed_vector_1[index] = index.Global;
            distributed_vector_2[index] = 1000 + index.Global;
        }
        distributed_vector_1.Restore();
        distributed_vector_2.Restore();

        // Write the vector
        writer.PutVector(node_id, petsc_data_1);
        writer.PutVector(ik_id, petsc_data_1);
        writer.PutVector(ina_id, petsc_data_2);

        writer.Close();

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_multi_column", true,
                                                "io/test/data", "hdf5_test_multi_column", false));

        // We make sure that the comparison returns the same result with the weaker method
        TS_ASSERT(CompareFilesViaHdf5DataReaderGlobalNorm("TestHdf5DataWriter", "hdf5_test_multi_column", true,
                                                          "io/test/data", "hdf5_test_multi_column", false));

        std::cout << "The next comparison should fail...\n";
        // We make that the comparison fails as the files are different
        TS_ASSERT(!CompareFilesViaHdf5DataReaderGlobalNorm("TestHdf5DataWriter", "hdf5_test_multi_column", true,
                                                           "io/test/data", "hdf5_test_full_format_extended", false));

        PetscTools::Destroy(petsc_data_1);
        PetscTools::Destroy(petsc_data_2);
    }

    void TestHdf5DataWriterSingleColumnNoTimeCachedFails()
    {
        int number_nodes = 100;

        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory,
                              "TestHdf5DataWriter",
                              "hdf5_test_single_column_cached_fails",
                              false,
                              false,
                              "Data",
                              true); // cache
        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node", "dimensionless");

        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        for (DistributedVector::Iterator index = distributed_vector_1.Begin();
             index != distributed_vector_1.End();
             ++index)
        {
            distributed_vector_1[index] = index.Global;
        }
        distributed_vector_1.Restore();

        TS_ASSERT_THROWS_THIS(writer.PutVector(node_id, petsc_data_1),
                              "Cached writes require an unlimited dimension.");

        PetscTools::Destroy(petsc_data_1);
    }

    void TestHdf5DataWriterMultipleColumnsCachedFails()
    {
        int number_nodes = 100;

        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory,
                              "TestHdf5DataWriter",
                              "hdf5_test_multi_column_cached_fails",
                              false,
                              false,
                              "Data",
                              true); // cache
        writer.DefineFixedDimension(number_nodes);

        // Define TWO variables
        int node_id = writer.DefineVariable("Node", "dimensionless");
        writer.DefineVariable("I_K", "milliamperes");

        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        for (DistributedVector::Iterator index = distributed_vector_1.Begin();
             index != distributed_vector_1.End();
             ++index)
        {
            distributed_vector_1[index] = index.Global;
        }
        distributed_vector_1.Restore();

        // Try and write ONE -> exception
        TS_ASSERT_THROWS_THIS(writer.PutVector(node_id, petsc_data_1),
                              "Cached writes must write all variables at once.");

        PetscTools::Destroy(petsc_data_1);
    }

    void TestHdf5DataWriterSingleColumnCached()
    {
        int number_nodes = 100;
        DistributedVectorFactory factory(number_nodes);
        std::string folder("TestHdf5DataWriter");
        std::string filename("hdf5_test_single_column_cached");

        // Write the original file
        {
            Hdf5DataWriter writer(factory,
                                  folder,
                                  filename,
                                  false,
                                  false, // extend
                                  "Data",
                                  true); // cache
            writer.DefineFixedDimension(number_nodes);

            int node_id = writer.DefineVariable("Node", "dimensionless");
            writer.DefineUnlimitedDimension("Time", "msec");
            writer.EndDefineMode();

            Vec petsc_data_1 = factory.CreateVec();
            DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index != distributed_vector_1.End();
                 ++index)
            {
                distributed_vector_1[index] = index.Global;
            }
            distributed_vector_1.Restore();
            writer.PutVector(node_id, petsc_data_1);
            writer.PutUnlimitedVariable(0.0);
            writer.AdvanceAlongUnlimitedDimension();
            writer.Close();
            PetscTools::Destroy(petsc_data_1);
        }

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", filename, true,
                                                "io/test/data", filename, false));

        // Now extend it
        {
            Hdf5DataWriter writer(factory,
                                  folder,
                                  filename,
                                  false,
                                  true, // extend
                                  "Data",
                                  true); // cache

            // Check chunk info was read correctly
            hsize_t expected_chunk_size[3] = { 1, 100, 1 };
            for (int i = 0; i < 3; ++i)
            {
                TS_ASSERT_EQUALS(writer.mChunkSize[i], expected_chunk_size[i]);
            }
            // Check the cache has reserved the right amount of space
            TS_ASSERT_EQUALS(writer.mDataCache.capacity(), writer.mNumberOwned);

            // Get IDs for the variables in the file
            int node_id = writer.GetVariableByName("Node");

            Vec petsc_data_1 = factory.CreateVec();
            DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index != distributed_vector_1.End();
                 ++index)
            {
                distributed_vector_1[index] = number_nodes * 2 + index.Global;
            }
            distributed_vector_1.Restore();
            writer.PutVector(node_id, petsc_data_1);
            writer.PutUnlimitedVariable(1.0);
            writer.AdvanceAlongUnlimitedDimension();
            writer.Close();
            PetscTools::Destroy(petsc_data_1);
        }

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", filename, true,
                                                "io/test/data", filename + "_extended", false));
    }

    void TestHdf5DataWriterNonEvenRowDistribution()
    {
        int number_nodes = 100;

        PetscInt local_number_of_nodes;

        if (PetscTools::AmMaster())
        {
            local_number_of_nodes = number_nodes - PetscTools::GetNumProcs() + 1;
        }
        else
        {
            local_number_of_nodes = 1;
        }

        DistributedVectorFactory factory(number_nodes, local_number_of_nodes);

        Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_non_even_row_dist", false);
        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        int ik_id = writer.DefineVariable("I_K", "milliamperes");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");

        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2 = factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        // Write some values
        for (DistributedVector::Iterator index = distributed_vector_1.Begin();
             index != distributed_vector_1.End();
             ++index)
        {
            distributed_vector_1[index] = index.Global;
            distributed_vector_2[index] = 1000 + index.Global;
        }
        distributed_vector_1.Restore();
        distributed_vector_2.Restore();

        // Write the vector
        writer.PutVector(node_id, petsc_data_1);
        writer.PutVector(ik_id, petsc_data_1);
        writer.PutVector(ina_id, petsc_data_2);

        writer.Close();

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_non_even_row_dist", true,
                                                "io/test/data", "hdf5_test_multi_column", false));

        PetscTools::Destroy(petsc_data_1);
        PetscTools::Destroy(petsc_data_2);
    }

    void TestHdf5DataWriterFullFormatIncomplete()
    {
        int number_nodes = 100;

        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_test_full_format_incomplete", false);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        int ik_id = writer.DefineVariable("I_K", "milliamperes");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");
        writer.DefineUnlimitedDimension("Time", "msec");

        std::vector<unsigned> node_numbers;
        node_numbers.push_back(21);
        node_numbers.push_back(47);
        node_numbers.push_back(60);
        writer.DefineFixedDimension(node_numbers, node_numbers, number_nodes);

        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2 = factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        Vec petsc_data_3 = factory.CreateVec();
        DistributedVector distributed_vector_3 = factory.CreateDistributedVector(petsc_data_3);

        for (unsigned time_step = 0; time_step < 10; time_step++)
        {
            // Write some values
            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index != distributed_vector_1.End();
                 ++index)
            {
                distributed_vector_1[index] = index.Global;
                distributed_vector_2[index] = time_step * 1000 + 100 + index.Global;
                distributed_vector_3[index] = time_step * 1000 + 200 + index.Global;
            }
            distributed_vector_1.Restore();
            distributed_vector_2.Restore();
            distributed_vector_3.Restore();

            // Write the vector

            writer.PutVector(node_id, petsc_data_1);
            writer.PutVector(ik_id, petsc_data_2);
            writer.PutVector(ina_id, petsc_data_3);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        writer.Close();

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_full_format_incomplete", true,
                                                "io/test/data", "hdf5_test_full_format_incomplete", false));

        // Test whether one with big-endian datatypes looks the same:
        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_full_format_incomplete", true,
                                                "io/test/data", "hdf5_test_full_format_incomplete_bigendian", false));

        PetscTools::Destroy(petsc_data_1);
        PetscTools::Destroy(petsc_data_2);
        PetscTools::Destroy(petsc_data_3);
    }

    void TestHdf5DataWriterFullFormatIncompleteCached()
    {
        int number_nodes = 100;

        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory,
                              "TestHdf5DataWriter",
                              "hdf5_test_full_format_incomplete_cached",
                              false,
                              false,
                              "Data",
                              true); // cache

        int node_id = writer.DefineVariable("Node", "dimensionless");
        writer.DefineUnlimitedDimension("Time", "msec");

        std::vector<unsigned> node_numbers;
        node_numbers.push_back(21);
        node_numbers.push_back(47);
        node_numbers.push_back(48);
        writer.DefineFixedDimension(node_numbers, node_numbers, number_nodes);

        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        for (unsigned time_step = 0; time_step < 10; time_step++)
        {
            // Write some values
            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index != distributed_vector_1.End();
                 ++index)
            {
                distributed_vector_1[index] = index.Global;
            }
            distributed_vector_1.Restore();

            // Write the vector

            writer.PutVector(node_id, petsc_data_1);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        writer.Close();

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_full_format_incomplete_cached", true,
                                                "io/test/data", "hdf5_test_full_format_incomplete_cached", false));

        PetscTools::Destroy(petsc_data_1);
    }

    void TestHdf5DataWriterFullFormat()
    {
        int number_nodes = 100;

        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_test_full_format", false);
        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        int ik_id = writer.DefineVariable("I_K", "milliamperes");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");
        writer.DefineUnlimitedDimension("Time", "msec", 10);

        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2 = factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        Vec petsc_data_3 = factory.CreateVec();
        DistributedVector distributed_vector_3 = factory.CreateDistributedVector(petsc_data_3);

        for (unsigned time_step = 0; time_step < 10; time_step++)
        {
            // Write some values
            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index != distributed_vector_1.End();
                 ++index)
            {
                distributed_vector_1[index] = index.Global;
                distributed_vector_2[index] = time_step * 1000 + 100 + index.Global;
                distributed_vector_3[index] = time_step * 1000 + 200 + index.Global;
            }
            distributed_vector_1.Restore();
            distributed_vector_2.Restore();
            distributed_vector_3.Restore();

            // Write the vector
            writer.PutVector(node_id, petsc_data_1);
            writer.PutVector(ik_id, petsc_data_2);
            writer.PutVector(ina_id, petsc_data_3);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        writer.Close();

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_full_format", true,
                                                "io/test/data", "hdf5_test_full_format", false));

        PetscTools::Destroy(petsc_data_1);
        PetscTools::Destroy(petsc_data_2);
        PetscTools::Destroy(petsc_data_3);
    }

    void TestHdf5DataWriterFullFormatStriped()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        Hdf5DataWriter writer(vec_factory, "TestHdf5DataWriter", "hdf5_test_full_format_striped", false);
        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        int vm_id = writer.DefineVariable("V_m", "millivolts");
        int phi_e_id = writer.DefineVariable("Phi_e", "millivolts");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");

        std::vector<int> striped_variable_IDs;
        striped_variable_IDs.push_back(vm_id);
        striped_variable_IDs.push_back(phi_e_id);

        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        DistributedVectorFactory factory(number_nodes);

        Vec petsc_data_short = vec_factory.CreateVec();
        DistributedVector distributed_vector_short = vec_factory.CreateDistributedVector(petsc_data_short);

        Vec node_number = vec_factory.CreateVec();
        DistributedVector distributed_node_number = vec_factory.CreateDistributedVector(node_number);

        for (DistributedVector::Iterator index = distributed_vector_short.Begin();
             index != distributed_vector_short.End();
             ++index)
        {
            distributed_node_number[index] = index.Global;
            distributed_vector_short[index] = -0.5;
        }
        distributed_node_number.Restore();
        distributed_vector_short.Restore();

        Vec petsc_data_long = factory.CreateVec(2);
        DistributedVector distributed_vector_long = factory.CreateDistributedVector(petsc_data_long);
        DistributedVector::Stripe vm_stripe(distributed_vector_long, 0);
        DistributedVector::Stripe phi_e_stripe(distributed_vector_long, 1);

        for (unsigned time_step = 0; time_step < 10; time_step++)
        {
            for (DistributedVector::Iterator index = distributed_vector_long.Begin();
                 index != distributed_vector_long.End();
                 ++index)
            {
                vm_stripe[index] = time_step * 1000 + index.Global * 2;
                phi_e_stripe[index] = time_step * 1000 + index.Global * 2 + 1;
            }
            distributed_vector_long.Restore();

            writer.PutVector(node_id, node_number);
            writer.PutVector(ina_id, petsc_data_short);
            writer.PutStripedVector(striped_variable_IDs, petsc_data_long);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        writer.Close();

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_full_format_striped", true,
                                                "io/test/data", "hdf5_test_full_format_striped", false));

        PetscTools::Destroy(node_number);
        PetscTools::Destroy(petsc_data_long);
        PetscTools::Destroy(petsc_data_short);
    }

    void TestHdf5DataWriterStripedCached()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        Hdf5DataWriter writer(vec_factory,
                              "TestHdf5DataWriter",
                              "hdf5_test_striped_with_cache",
                              false,
                              false,
                              "Data",
                              true); // use cache
        writer.DefineFixedDimension(number_nodes);

        /* Set specific chunk dims for coverage.
         * We expect that the writer will flush whole chunks (every 3 steps in
         * this case) automatically, but we have 10 entries, so Close() will do
         * for the final flush. */
        writer.SetFixedChunkSize(3, 10, 2);

        int vm_id = writer.DefineVariable("V_m", "millivolts");
        int phi_e_id = writer.DefineVariable("Phi_e", "millivolts");

        std::vector<int> striped_variable_IDs;
        striped_variable_IDs.push_back(vm_id);
        striped_variable_IDs.push_back(phi_e_id);

        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        // Check chunk info was read correctly
        hsize_t expected_chunk_size[3] = { 3, 10, 2 };
        for (int i = 0; i < 3; ++i)
        {
            TS_ASSERT_EQUALS(writer.mChunkSize[i], expected_chunk_size[i]);
        }
        // Check the cache has reserved the right amount of space
        unsigned expected_capacity = 3 * writer.mNumberOwned * 2;
        TS_ASSERT_EQUALS(writer.mDataCache.capacity(), expected_capacity);

        DistributedVectorFactory factory(number_nodes);

        Vec petsc_data_long = factory.CreateVec(2);
        DistributedVector distributed_vector_long = factory.CreateDistributedVector(petsc_data_long);
        DistributedVector::Stripe vm_stripe(distributed_vector_long, 0);
        DistributedVector::Stripe phi_e_stripe(distributed_vector_long, 1);

        for (unsigned time_step = 0; time_step < 10; time_step++)
        {
            for (DistributedVector::Iterator index = distributed_vector_long.Begin();
                 index != distributed_vector_long.End();
                 ++index)
            {
                vm_stripe[index] = time_step * 1000 + index.Global * 2;
                phi_e_stripe[index] = time_step * 1000 + index.Global * 2 + 1;
            }
            distributed_vector_long.Restore();

            writer.PutStripedVector(striped_variable_IDs, petsc_data_long);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();

            // Check that the cache is emptied on whole chunks. Size of cache
            // should go 0, 200, 400, 0, ... when run with one process.
            unsigned expected_cache_size = ((time_step + 1) % 3) * writer.mNumberOwned * 2;
            TS_ASSERT_EQUALS(writer.mDataCache.size(), expected_cache_size);
        }

        // Final flush happens here
        writer.Close();

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_striped_with_cache", true,
                                                "io/test/data", "hdf5_test_striped_with_cache", false));

        PetscTools::Destroy(petsc_data_long);
    }

    void TestHdf5DataWriterStripedNoTimeCachedFails()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        Hdf5DataWriter writer(vec_factory,
                              "TestHdf5DataWriter",
                              "hdf5_test_striped_no_time_cache",
                              false,
                              false,
                              "Data",
                              true); // use cache
        writer.DefineFixedDimension(number_nodes);

        int vm_id = writer.DefineVariable("V_m", "millivolts");
        int phi_e_id = writer.DefineVariable("Phi_e", "millivolts");

        std::vector<int> striped_variable_IDs;
        striped_variable_IDs.push_back(vm_id);
        striped_variable_IDs.push_back(phi_e_id);

        writer.EndDefineMode();

        DistributedVectorFactory factory(number_nodes);

        Vec petsc_data_long = factory.CreateVec(2);
        DistributedVector distributed_vector_long = factory.CreateDistributedVector(petsc_data_long);
        DistributedVector::Stripe vm_stripe(distributed_vector_long, 0);
        DistributedVector::Stripe phi_e_stripe(distributed_vector_long, 1);

        for (DistributedVector::Iterator index = distributed_vector_long.Begin();
             index != distributed_vector_long.End();
             ++index)
        {
            vm_stripe[index] = 1000 + index.Global * 2;
            phi_e_stripe[index] = 1000 + index.Global * 2 + 1;
        }
        distributed_vector_long.Restore();

        TS_ASSERT_THROWS_THIS(writer.PutStripedVector(striped_variable_IDs, petsc_data_long),
                              "Cached writes require an unlimited dimension.")

        PetscTools::Destroy(petsc_data_long);
    }

    void TestHdf5DataWriterFullFormatStripedWith3Variables()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        Hdf5DataWriter writer(vec_factory, "TestHdf5DataWriter", "hdf5_test_full_format_striped_3vars", false);
        writer.DefineFixedDimension(number_nodes);

        int vm_id = writer.DefineVariable("V_m", "millivolts");
        int phi_e_id = writer.DefineVariable("Phi_e", "millivolts");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");

        std::vector<int> striped_variable_IDs;
        striped_variable_IDs.push_back(vm_id);
        striped_variable_IDs.push_back(phi_e_id);
        striped_variable_IDs.push_back(ina_id);

        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        DistributedVectorFactory factory(number_nodes);

        Vec petsc_data = factory.CreateVec(3);
        DistributedVector distributed_vector = factory.CreateDistributedVector(petsc_data);
        DistributedVector::Stripe vm_stripe(distributed_vector, 0);
        DistributedVector::Stripe phi_e_stripe(distributed_vector, 1);
        DistributedVector::Stripe ina_stripe(distributed_vector, 2);

        for (unsigned time_step = 0; time_step < 10; time_step++)
        {
            for (DistributedVector::Iterator index = distributed_vector.Begin();
                 index != distributed_vector.End();
                 ++index)
            {
                vm_stripe[index] = time_step * 1000 + index.Global * 2;
                phi_e_stripe[index] = time_step * 1000 + index.Global * 2 + 1;
                ina_stripe[index] = -56.0;
            }
            distributed_vector.Restore();

            writer.PutStripedVector(striped_variable_IDs, petsc_data);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }
        writer.Close();

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_full_format_striped_3vars", true,
                                                "io/test/data", "hdf5_test_full_format_striped_3vars", false));
        PetscTools::Destroy(petsc_data);
    }

    void Hdf5DataWriterFullFormatStripedIncomplete(bool useCache, std::string outputFile, std::string expectedException)
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        Hdf5DataWriter writer(vec_factory,
                              "TestHdf5DataWriter",
                              outputFile,
                              false,
                              false,
                              "Data",
                              useCache);

        std::vector<unsigned> node_numbers;
        node_numbers.push_back(4);
        node_numbers.push_back(8);
        node_numbers.push_back(15);
        node_numbers.push_back(16);
        node_numbers.push_back(23);
        node_numbers.push_back(42);
        writer.DefineFixedDimension(node_numbers, node_numbers, number_nodes);

        int vm_id = writer.DefineVariable("V_m", "millivolts");
        int phi_e_id = writer.DefineVariable("Phi_e", "millivolts");

        std::vector<int> variable_IDs;
        variable_IDs.push_back(vm_id);
        variable_IDs.push_back(phi_e_id);

        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        DistributedVectorFactory factory(number_nodes);
        Vec petsc_data_long = factory.CreateVec(2);
        DistributedVector distributed_vector_long = factory.CreateDistributedVector(petsc_data_long);
        DistributedVector::Stripe vm_stripe(distributed_vector_long, 0);
        DistributedVector::Stripe phi_e_stripe(distributed_vector_long, 1);

        for (unsigned time_step = 0; time_step < 2; time_step++)
        {
            for (DistributedVector::Iterator index = distributed_vector_long.Begin();
                 index != distributed_vector_long.End();
                 ++index)
            {
                vm_stripe[index] = (time_step + 1) * 1000 + index.Global;
                phi_e_stripe[index] = index.Global;
            }
            distributed_vector_long.Restore();

            writer.PutStripedVector(variable_IDs, petsc_data_long);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }
        writer.Close();

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", outputFile, true,
                                                "io/test/data", "hdf5_test_full_format_striped_incomplete", false));

        PetscTools::Destroy(petsc_data_long);

        // Now cover two exceptions: one is the unsupported PutStripedVector for incomplete data and 3 vars...
        int first = writer.DefineVariable("first", "millivolts");
        int second = writer.DefineVariable("second", "millivolts");
        int third = writer.DefineVariable("third", "milliAmps");

        std::vector<int> three_variable_IDs;
        three_variable_IDs.push_back(first);
        three_variable_IDs.push_back(second);
        three_variable_IDs.push_back(third);
        writer.EndDefineMode();

        Vec petsc_data_3vars = factory.CreateVec(3);
        TS_ASSERT_THROWS_THIS(writer.PutStripedVector(three_variable_IDs, petsc_data_3vars), expectedException);

        PetscTools::Destroy(petsc_data_3vars);
        writer.Close();

        // ...and one is the case when we pass in a short vector to the PutStripedVector method
        int single_var = writer.DefineVariable("only", "one");
        std::vector<int> one_ID;
        one_ID.push_back(single_var);
        writer.EndDefineMode();
        Vec petsc_data_1var = factory.CreateVec(1);
        TS_ASSERT_THROWS_THIS(writer.PutStripedVector(one_ID, petsc_data_1var),
                              "The PutStripedVector method requires at least two variables ID. If only one is needed, use PutVector method instead");
        PetscTools::Destroy(petsc_data_1var);
    }

    void TestHdf5DataWriterFullFormatStripedIncomplete()
    {
        Hdf5DataWriterFullFormatStripedIncomplete(false,
                                                  "hdf5_test_full_format_striped_incomplete",
                                                  "The PutStripedVector functionality for incomplete data is supported for only 2 stripes");
    }

    void TestHdf5DataWriterFullFormatStripedIncompleteCached()
    {
        Hdf5DataWriterFullFormatStripedIncomplete(true,
                                                  "hdf5_test_full_format_striped_incomplete_cached",
                                                  "Cached writes must write all variables at once.");
    }

    void TestNonImplementedFeatures()
    {
        int number_nodes = 100;
        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_test_non_implemented", false);
        writer.DefineFixedDimension(number_nodes);

        int vm_id = writer.DefineVariable("V_m", "millivolts");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");
        int phi_e_id = writer.DefineVariable("Phi_e", "millivolts");

        std::vector<int> unoredred_variable_IDs;
        unoredred_variable_IDs.push_back(vm_id);
        unoredred_variable_IDs.push_back(phi_e_id);

        std::vector<int> ordered_variable_IDs;
        ordered_variable_IDs.push_back(vm_id);
        ordered_variable_IDs.push_back(ina_id);

        writer.EndDefineMode();

        Vec petsc_data_short = factory.CreateVec();
        DistributedVector distributed_vector_short = factory.CreateDistributedVector(petsc_data_short);

        for (DistributedVector::Iterator index = distributed_vector_short.Begin();
             index != distributed_vector_short.End();
             ++index)
        {
            distributed_vector_short[index] = -0.5;
        }
        distributed_vector_short.Restore();

        DistributedVectorFactory factory2(2 * number_nodes);
        Vec petsc_data_long = factory2.CreateVec();
        DistributedVector distributed_vector_long = factory.CreateDistributedVector(petsc_data_long);
        DistributedVector::Stripe vm_stripe(distributed_vector_long, 0);
        for (DistributedVector::Iterator index = distributed_vector_long.Begin();
             index != distributed_vector_long.End();
             ++index)
        {
            vm_stripe[index] = index.Global;
        }
        distributed_vector_long.Restore();

        writer.PutVector(ina_id, petsc_data_short);
        //Try to write striped data in the wrong columns
        TS_ASSERT_THROWS_THIS(writer.PutStripedVector(unoredred_variable_IDs, petsc_data_long),
                              "Columns should be consecutive. Try reordering them.");
        //Try to write data of wrong size
        TS_ASSERT_THROWS_THIS(writer.PutVector(ina_id, petsc_data_long),
                              "Vector size doesn\'t match fixed dimension");
        TS_ASSERT_THROWS_THIS(writer.PutStripedVector(ordered_variable_IDs, petsc_data_short),
                              "Vector size doesn\'t match fixed dimension");

        writer.Close();

        PetscTools::Destroy(petsc_data_long);
        PetscTools::Destroy(petsc_data_short);
    }

    /**
     * Tests copied (with some minor modifications) from TestColumnDataReaderWriter:
     * to be refactored at some point.
     */
    void TestDefineThings()
    {
        DistributedVectorFactory vec_factory(100);

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new Hdf5DataWriter(vec_factory, "", "test"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"),
                              "Unlimited dimension already set. Cannot be defined twice");

        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("Time", "m secs"),
                              "Unlimited dimension already set. Cannot be defined twice");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("T,i,m,e", "msecs"),
                              "Unlimited dimension already set. Cannot be defined twice");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("", "msecs"),
                              "Unlimited dimension already set. Cannot be defined twice");

        std::vector<unsigned> node_numbers;
        node_numbers.push_back(21);
        node_numbers.push_back(47);
        node_numbers.push_back(6);

        // Data not increasing
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension(node_numbers, node_numbers, 100), "Input should be monotonic increasing");
        node_numbers[2] = 100;
        // Data is increasing but the last number is too large
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension(node_numbers, node_numbers, 100), "Vector size doesn\'t match nodes to output");

        mpTestWriter->DefineFixedDimension(5000);
        // Can't set fixed dimension more than once
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension(5000), "Fixed dimension already set");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension(node_numbers, node_numbers, 100), "Vector size doesn\'t match nodes to output");

        int ina_var_id = INT_UNSET;
        int ik_var_id = INT_UNSET;
        int ik2_var_id = INT_UNSET;

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineVariable("Dummy", ""));

        // Defined twice
        TS_ASSERT_THROWS_THIS(ik2_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"),
                              "Variable name already exists");
        TS_ASSERT_EQUALS(ik2_var_id, INT_UNSET);

        // Bad variable names/units
        TS_ASSERT_THROWS_THIS(ik_var_id = mpTestWriter->DefineVariable("I_K", "milli amperes"),
                              "Variable name/units \'milli amperes\' not allowed: may only contain alphanumeric characters or \'_\'.");
        TS_ASSERT_THROWS_THIS(ik_var_id = mpTestWriter->DefineVariable("I   K", "milliamperes"),
                              "Variable name/units \'I   K\' not allowed: may only contain alphanumeric characters or \'_\'.");
        TS_ASSERT_THROWS_THIS(ik_var_id = mpTestWriter->DefineVariable("I.K", "milliamperes"),
                              "Variable name/units \'I.K\' not allowed: may only contain alphanumeric characters or \'_\'.");
        TS_ASSERT_THROWS_THIS(ik_var_id = mpTestWriter->DefineVariable("", "milliamperes"),
                              "Variable name not allowed: may not be blank.");

        TS_ASSERT_EQUALS(ina_var_id, 0);
        TS_ASSERT_EQUALS(ik_var_id, 1);

        delete mpTestWriter;
    }

    void TestEndDefineMode()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new Hdf5DataWriter(vec_factory, "", "testdefine"));

        // Ending define mode without having defined at least a variable and a fixed dimension should raise an exception
        TS_ASSERT_THROWS_THIS(mpTestWriter->EndDefineMode(), "Cannot end define mode. No variables have been defined.");

        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"));
        TS_ASSERT_THROWS_THIS(mpTestWriter->EndDefineMode(),
                              "Cannot end define mode. No variables have been defined.");

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_EQUALS(ina_var_id, 0);
        TS_ASSERT_EQUALS(ik_var_id, 1);

        // In Hdf5 a fixed dimension should be defined always
        TS_ASSERT_THROWS_THIS(mpTestWriter->EndDefineMode(),
                              "Cannot end define mode. One fixed dimension should be defined.");

        std::vector<unsigned> node_numbers;
        node_numbers.push_back(21);
        node_numbers.push_back(47);
        node_numbers.push_back(60);
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension(node_numbers, node_numbers, 100));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineVariable("I_Ca", "milli amperes"),
                              "Cannot define variables when not in Define mode");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"),
                              "Unlimited dimension already set. Cannot be defined twice");
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension(5000),
                              "Cannot define variables when not in Define mode");

        // Can't call define fixed dimension again
        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension(node_numbers, node_numbers, 100),
                              "Cannot define variables when not in Define mode");

        // Test that we can't write incomplete data from a vector that doesn't have the right entries (0 to 59)
        DistributedVectorFactory factory(60);
        Vec petsc_data_short = factory.CreateVec();
        TS_ASSERT_THROWS_THIS(mpTestWriter->PutVector(0, petsc_data_short),
                              "Vector size doesn\'t match fixed dimension");
        PetscTools::Destroy(petsc_data_short);

        mpTestWriter->Close();
        delete mpTestWriter;
    }

    void TestCantAddUnlimitedAfterEndDefine()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new Hdf5DataWriter(vec_factory, "", "testdefine"));
        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineFixedDimension(0), "Fixed dimension must be at least 1 long");
        mpTestWriter->DefineFixedDimension(5000);

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K", "milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        TS_ASSERT_EQUALS(ina_var_id, 0);
        TS_ASSERT_EQUALS(ik_var_id, 1);

        TS_ASSERT_THROWS_THIS(mpTestWriter->DefineUnlimitedDimension("Time", "msecs"),
                              "Cannot define variables when not in Define mode");
        mpTestWriter->Close();
        delete mpTestWriter;
    }

    void TestAdvanceAlongUnlimitedDimension()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new Hdf5DataWriter(vec_factory, "", "testdefine"));

        int ina_var_id = 9999;
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension(5000));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na", "milliamperes"));
        TS_ASSERT_EQUALS(ina_var_id, 0);

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_THIS(mpTestWriter->PutUnlimitedVariable(0.0),
                              "PutUnlimitedVariable() called but no unlimited dimension has been set");
        TS_ASSERT_THROWS_THIS(mpTestWriter->AdvanceAlongUnlimitedDimension(),
                              "Trying to advance along an unlimited dimension without having defined any");

        mpTestWriter->Close();
        delete mpTestWriter;
    }

    void TestCantWriteDataWhileInDefineMode()
    {
        int number_nodes = 100;
        DistributedVectorFactory vec_factory(number_nodes);

        Hdf5DataWriter writer(vec_factory, "", "testdefine", false);
        //        writer.DefineFixedDimension(number_nodes);
        //
        //        int node_id = writer.DefineVariable("Node","dimensionless");
        //        int vm_id = writer.DefineVariable("V_m","millivolts");
        //        int phi_e_id = writer.DefineVariable("Phi_e","millivolts");
        //        int ina_id = writer.DefineVariable("I_Na","milliamperes");
        //
        //        writer.DefineUnlimitedDimension("Time", "msec");

        //        writer.EndDefineMode();

        int node_id = 1;
        int vm_id = 2;
        int phi_e_id = 3;
        int ina_id = 4;

        std::vector<int> variable_IDs;
        variable_IDs.push_back(vm_id);
        variable_IDs.push_back(phi_e_id);

        DistributedVectorFactory factory(number_nodes);

        Vec petsc_data_short = vec_factory.CreateVec();
        DistributedVector distributed_vector_short = vec_factory.CreateDistributedVector(petsc_data_short);

        Vec node_number = vec_factory.CreateVec();
        DistributedVector distributed_node_number = vec_factory.CreateDistributedVector(node_number);

        for (DistributedVector::Iterator index = distributed_node_number.Begin();
             index != distributed_node_number.End();
             ++index)
        {
            distributed_node_number[index] = index.Global;
            distributed_vector_short[index] = -0.5;
        }
        distributed_node_number.Restore();
        distributed_vector_short.Restore();

        DistributedVectorFactory factory2(2 * number_nodes);
        Vec petsc_data_long = factory2.CreateVec();
        DistributedVector distributed_vector_long = factory2.CreateDistributedVector(petsc_data_long);

        for (DistributedVector::Iterator index = distributed_vector_long.Begin();
             index != distributed_vector_long.End();
             ++index)
        {
            distributed_vector_long[index] = 1000 + index.Global;
        }
        distributed_vector_long.Restore();

        TS_ASSERT_THROWS_THIS(writer.PutVector(node_id, node_number),
                              "Cannot write data while in define mode.");
        TS_ASSERT_THROWS_THIS(writer.PutVector(ina_id, petsc_data_short),
                              "Cannot write data while in define mode.");
        TS_ASSERT_THROWS_THIS(writer.PutStripedVector(variable_IDs, petsc_data_long),
                              "Cannot write data while in define mode.");
        TS_ASSERT_THROWS_THIS(writer.PutUnlimitedVariable(0.0),
                              "Cannot write data while in define mode.");
        TS_ASSERT_THROWS_THIS(writer.AdvanceAlongUnlimitedDimension(),
                              "Trying to advance along an unlimited dimension without having defined any");

        writer.Close();
        PetscTools::Destroy(petsc_data_short);
        PetscTools::Destroy(node_number);
        PetscTools::Destroy(petsc_data_long);
    }

    /**
     * Test the functionality for adding further data to an existing file.
     */
    void TestFailCreateFile(void)
    {
        // This test causes memory leak within MPIO
        int number_nodes = 100;
        DistributedVectorFactory factory(number_nodes);

        OutputFileHandler handler("TestHdf5DataWriter", false);

        // Create file and remove permission to overwrite file
        handler.OpenOutputFile("empty.h5")->close();
        PetscTools::Barrier("TestFailCreateFile"); //Wait for all processes to leave file before checking that it's there
        FileFinder empty = handler.FindFile("empty.h5");
        TS_ASSERT(empty.Exists());
        chmod(empty.GetAbsolutePath().c_str(), CHASTE_READONLY);

        Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "empty", false);
        writer.DefineVariable("Node", "dimensionless");
        writer.DefineFixedDimension(number_nodes);
        H5E_BEGIN_TRY //Suppress HDF5 error in this test
        {
            TS_ASSERT_THROWS_CONTAINS(writer.EndDefineMode(), "Hdf5DataWriter could not create");
        }
        H5E_END_TRY;
        writer.Close();

        // Re-instate permission to overwrite file
        chmod(empty.GetAbsolutePath().c_str(), CHASTE_READ_WRITE);
    }

    /**
     * Test the functionality for adding a new dataset ("Postprocessing") to an existing file.
     *
     * This test must come after TestHdf5DataWriterFullFormat and TestHdf5DataWriterFullFormatStripedIncomplete,
     * as we extend their files.
     *
     * NB And this test must come before TestWriteToExistingFile which test that
     * there is a Postprocessing block!
     */
    void TestWriteNewDatasetToExistingFile(void)
    {
        int number_nodes = 100;
        DistributedVectorFactory factory(number_nodes);

        // Open the real file
        Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_test_full_format", false, true, "Postprocessing");

        // Can't set alignment on existing file (with pre-existing dataset or not).
        TS_ASSERT_THROWS_THIS(writer.SetAlignment(123), "Alignment parameter can only be set for new HDF5 files.");

        // CAN set chunk size target on new dataset in existing file
        writer.SetTargetChunkSize(0x800); // 2 K

        // Define what the new dataset is going to look like.
        writer.DefineFixedDimension(number_nodes);
        writer.DefineVariable("Phase", "dimensionless");
        writer.DefineUnlimitedDimension("Time", "msec", 10);
        writer.EndDefineMode();

        // Get IDs for the variables in the file
        int phase_id = writer.GetVariableByName("Phase");

        // Create some extra test data
        Vec phase_petsc = factory.CreateVec();
        DistributedVector phase_data = factory.CreateDistributedVector(phase_petsc);

        for (unsigned time_step = 0; time_step < 10; time_step++)
        {
            // Fill in data
            for (DistributedVector::Iterator index = phase_data.Begin();
                 index != phase_data.End();
                 ++index)
            {
                phase_data[index] = index.Global + 0.5;
            }
            phase_data.Restore();

            // Write to file
            writer.PutVector(phase_id, phase_petsc);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        // Close and test
        writer.Close();
        PetscTools::Destroy(phase_petsc);

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_full_format", true,
                                                "io/test/data", "hdf5_test_full_format_extended", false, 1e-10, "Postprocessing"));

        // Check chunk dimensions are as expected
        OutputFileHandler file_handler("TestHdf5DataWriter", false);
        FileFinder file = file_handler.FindFile("hdf5_test_full_format.h5");
        hid_t h5_file = H5Fopen(file.GetAbsolutePath().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t dset = H5Dopen(h5_file, "Postprocessing", H5P_DEFAULT); // open dataset
        hid_t dcpl = H5Dget_create_plist(dset); // get dataset creation property list
        hsize_t expected_dims[3] = { 10, 25, 1 };
        hsize_t chunk_dims[3];
        H5Pget_chunk(dcpl, 3, chunk_dims);
        for (int i = 0; i < 3; ++i)
        {
            TS_ASSERT_EQUALS(chunk_dims[i], expected_dims[i]);
        }
        H5Pclose(dcpl);
        H5Dclose(dset);
        H5Fclose(h5_file);
    }

    /**
     *  Test for adding a new dataset ("Extra stuff") to an existing HDF5 file.
     */
    void TestHdf5DataWriterAddNewVariable()
    {
        int number_nodes = 100;

        DistributedVectorFactory factory(number_nodes);

        /* Make a simple file so that we can add extra data later.  This block is so that this
         * test doesn't have to depend on a previous test.
         */
        {
            Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_test_adding_variables", false);

            writer.DefineFixedDimension(number_nodes);

            int node_id = writer.DefineVariable("Node", "dimensionless");
            int ik_id = writer.DefineVariable("I_K", "milliamperes");
            writer.DefineUnlimitedDimension("Time", "msec", 10);

            writer.EndDefineMode();

            Vec petsc_data_1 = factory.CreateVec();
            DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

            Vec petsc_data_2 = factory.CreateVec();
            DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

            for (unsigned time_step = 0; time_step < 10; time_step++)
            {
                // Write some values
                for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                     index != distributed_vector_1.End();
                     ++index)
                {
                    distributed_vector_1[index] = index.Global;
                    distributed_vector_2[index] = time_step * 1000 + 100 + index.Global;
                }
                distributed_vector_1.Restore();
                distributed_vector_2.Restore();

                // Write the vector
                writer.PutVector(node_id, petsc_data_1);
                writer.PutVector(ik_id, petsc_data_2);
                writer.PutUnlimitedVariable(time_step);
                writer.AdvanceAlongUnlimitedDimension();
            }

            writer.Close();
            PetscTools::Destroy(petsc_data_1);
            PetscTools::Destroy(petsc_data_2);
        }

        /* Re-open the file and add the new data */

        TS_ASSERT_THROWS_THIS(Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_test_adding_variables", false,
                                                    false /*(do not extend)*/, "Extra stuff"),
                              "Adding new data only makes sense when extending an existing file");
        Hdf5DataWriter annotating_writer(factory, "TestHdf5DataWriter", "hdf5_test_adding_variables", false, true, "Extra stuff");

        int phase_id = annotating_writer.DefineVariable("Phase", "radians");
        int plasma_id = annotating_writer.DefineVariable("Plasma", "gloops");
        int node_id_again = annotating_writer.DefineVariable("Node", "dimensionless");
        annotating_writer.DefineFixedDimension(number_nodes);
        annotating_writer.DefineUnlimitedDimension("Distance", "LightYears", 2);
        annotating_writer.EndDefineMode();
        {
            Vec petsc_data_1 = factory.CreateVec();
            DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

            Vec petsc_data_2 = factory.CreateVec();
            DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

            Vec petsc_data_3 = factory.CreateVec();
            DistributedVector distributed_vector_3 = factory.CreateDistributedVector(petsc_data_2);

            for (unsigned time_step = 0; time_step < 2; time_step++)
            {
                // Write some values
                for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                     index != distributed_vector_1.End();
                     ++index)
                {
                    distributed_vector_3[index] = index.Global;
                    distributed_vector_1[index] = index.Global;
                    distributed_vector_2[index] = time_step * 1000 + 100 + index.Global;
                }
                distributed_vector_1.Restore();
                distributed_vector_2.Restore();
                distributed_vector_3.Restore();

                // Write the vector
                annotating_writer.PutVector(phase_id, petsc_data_1);
                annotating_writer.PutVector(plasma_id, petsc_data_2);
                annotating_writer.PutVector(node_id_again, petsc_data_3);
                annotating_writer.PutUnlimitedVariable(time_step);
                annotating_writer.AdvanceAlongUnlimitedDimension();
            }

            annotating_writer.Close();
            PetscTools::Destroy(petsc_data_1);
            PetscTools::Destroy(petsc_data_2);
            PetscTools::Destroy(petsc_data_3);
        }

        // This one has the wrong name on the unlimited variable of the Extra Stuff dataset.
        TS_ASSERT_EQUALS(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_adding_variables", true,
                                                       "io/test/data", "hdf5_test_adding_variables_bad", false, 1e-10, "Extra stuff"),
                         false);

        // The 'bad' file isn't actually bad for the original Data (voltage etc. just the units of the new 'Extra Stuff').
        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_adding_variables", true,
                                                "io/test/data", "hdf5_test_adding_variables_bad", false));

        // This one is correct for both original and extra stuff datasets.
        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_adding_variables", true,
                                                "io/test/data", "hdf5_test_adding_variables", false));

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_adding_variables", true,
                                                "io/test/data", "hdf5_test_adding_variables", false, 1e-10, "Extra stuff"));
    }

    /**
     * Test the functionality for adding further data to an existing dataset in an existing file.
     *
     * This test must come after TestHdf5DataWriterFullFormat and TestHdf5DataWriterFullFormatStripedIncomplete,
     * as we extend their files.
     */
    void TestWriteToExistingFile(void)
    {
        int number_nodes = 100;
        DistributedVectorFactory factory(number_nodes);

        // Test some exceptions
        {
            TS_ASSERT_THROWS_THIS(Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_test_full_format", true, true),
                                  "You are asking to delete a file and then extend it, change arguments to constructor.");

            // Check for a wrong file format exception
            OutputFileHandler file_handler("TestHdf5DataWriter", false);
            out_stream p_wrong_format = file_handler.OpenOutputFile("hdf5_wrong_format.h5");
            *p_wrong_format << "gobbledegook" << std::endl;
            p_wrong_format->close();
            H5E_BEGIN_TRY //Supress HDF5 error in this test
            {
                TS_ASSERT_THROWS_CONTAINS(Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_wrong_format", false, true),
                                          "H5Fopen error code");
            }
            H5E_END_TRY;
            TS_ASSERT_THROWS_CONTAINS(Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "absent_file", false, true),
                                      "as it does not exist");

            TS_ASSERT_THROWS_THIS(Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_test_multi_column", false, true),
                                  "Tried to open a datafile for extending which doesn't have an unlimited dimension.");
        }

        // Open the real file
        Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_test_full_format", false, true);

        // Get IDs for the variables in the file
        int node_id = writer.GetVariableByName("Node");
        int ik_id = writer.GetVariableByName("I_K");
        int ina_id = writer.GetVariableByName("I_Na");

        TS_ASSERT_THROWS_THIS(writer.GetVariableByName("bob"),
                              "Variable does not exist in hdf5 definitions.");

        // Can't set chunk size on existing dataset
        TS_ASSERT_THROWS_THIS(writer.SetTargetChunkSize(123),
                              "Cannot set chunk target size when not in define mode.");

        // Can't set alignment on existing file (with pre-existing dataset or not).
        TS_ASSERT_THROWS_THIS(writer.SetAlignment(456), "Alignment parameter can only be set for new HDF5 files.");

        // Create some extra test data
        Vec node_petsc = factory.CreateVec();
        Vec ik_petsc = factory.CreateVec();
        Vec ina_petsc = factory.CreateVec();
        DistributedVector node_data = factory.CreateDistributedVector(node_petsc);
        DistributedVector ik_data = factory.CreateDistributedVector(ik_petsc);
        DistributedVector ina_data = factory.CreateDistributedVector(ina_petsc);

        for (unsigned time_step = 10; time_step < 15; time_step++)
        {
            // Fill in data
            for (DistributedVector::Iterator index = node_data.Begin();
                 index != node_data.End();
                 ++index)
            {
                node_data[index] = index.Global;
                ik_data[index] = time_step * 1000 + 100 + index.Global;
                ina_data[index] = time_step * 1000 + 200 + index.Global;
            }
            node_data.Restore();
            ik_data.Restore();
            ina_data.Restore();

            // Write to file
            writer.PutVector(node_id, node_petsc);
            writer.PutVector(ina_id, ina_petsc);
            writer.PutVector(ik_id, ik_petsc);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        // Close and test
        writer.Close();
        PetscTools::Destroy(node_petsc);
        PetscTools::Destroy(ik_petsc);
        PetscTools::Destroy(ina_petsc);

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_full_format", true,
                                                "io/test/data", "hdf5_test_full_format_extended", false));

        TS_ASSERT_THROWS_THIS(Hdf5DataWriter another_writer(factory, "TestHdf5DataWriter", "hdf5_test_full_format_striped_incomplete", false, true),
                              "Unable to extend an incomplete data file at present.");
    }

    /**
     * This test must come after TestWriteToExistingFile as it extends even further.
     */
    void TestWriteToExistingFileWithCache(void)
    {
        int number_nodes = 100;
        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory,
                              "TestHdf5DataWriter",
                              "hdf5_test_full_format",
                              false,
                              true, // extend
                              "Data",
                              true); // cache

        // Check chunk info was read correctly
        hsize_t expected_chunk_size[3] = { 10, 100, 3 };
        for (int i = 0; i < 3; ++i)
        {
            TS_ASSERT_EQUALS(writer.mChunkSize[i], expected_chunk_size[i]);
        }
        // Check the cache has reserved the right amount of space
        unsigned expected_capacity = 10 * writer.mNumberOwned * 3;
        TS_ASSERT_EQUALS(writer.mDataCache.capacity(), expected_capacity);

        // Get IDs for the variables in the file
        int node_id = writer.GetVariableByName("Node");
        int ik_id = writer.GetVariableByName("I_K");
        int ina_id = writer.GetVariableByName("I_Na");

        std::vector<int> variable_IDs;
        variable_IDs.push_back(node_id);
        variable_IDs.push_back(ik_id);
        variable_IDs.push_back(ina_id);

        Vec petsc_data_long = factory.CreateVec(3);
        DistributedVector distributed_vector_long = factory.CreateDistributedVector(petsc_data_long);
        DistributedVector::Stripe node_stripe(distributed_vector_long, 0);
        DistributedVector::Stripe ik_stripe(distributed_vector_long, 1);
        DistributedVector::Stripe ina_stripe(distributed_vector_long, 2);

        for (unsigned time_step = 15; time_step < 20; time_step++)
        {
            // Fill in data
            for (DistributedVector::Iterator index = distributed_vector_long.Begin();
                 index != distributed_vector_long.End();
                 ++index)
            {
                node_stripe[index] = index.Global;
                ik_stripe[index] = time_step * 1000 + 100 + index.Global;
                ina_stripe[index] = time_step * 1000 + 200 + index.Global;
            }
            distributed_vector_long.Restore();

            // Write to file
            writer.PutStripedVector(variable_IDs, petsc_data_long);
            writer.PutUnlimitedVariable(time_step);

            // Check that the cache is growing at the expected rate.
            // Should go 300, 600, 900, ... when run with one process.
            unsigned expected_cache_size = ((time_step - 15 + 1) % 10) * writer.mNumberOwned * 3;
            TS_ASSERT_EQUALS(writer.mDataCache.size(), expected_cache_size);

            writer.AdvanceAlongUnlimitedDimension();

            /*
             * We started halfway through a chunk, so after the final
             * iteration there should be a flush (despite only having a half-
             * full cache).
             */
            if (time_step == 19)
            {
                TS_ASSERT_EQUALS(writer.mDataCache.size(), 0u);
            }
        }

        // Close and test
        writer.Close();
        PetscTools::Destroy(petsc_data_long);

        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_test_full_format", true,
                                                "io/test/data", "hdf5_test_full_format_extended_twice", false));
    }

    void TestPermutation()
    {
        int number_nodes = 10;
        //Note that we need a 10x10 permutation matrix which is difficult to share over more processes
        if (PetscTools::GetNumProcs() > 10u)
        {
            TS_TRACE("This test is designed for fewer processes!");
            return;
        }
        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory, "TestHdf5DataWriter", "hdf5_permuted", false);
        writer.DefineFixedDimension(number_nodes);

        int index_id = writer.DefineVariable("index", "dimensionless");
        int vm_id = writer.DefineVariable("V_m", "millivolts");
        int phi_e_id = writer.DefineVariable("Phi_e", "millivolts");

        std::vector<int> variable_IDs;
        variable_IDs.push_back(vm_id);
        variable_IDs.push_back(phi_e_id);

        std::vector<unsigned> rotation_perm;
        std::vector<unsigned> identity_perm;
        std::vector<unsigned> short_perm;

        //Can't apply empty permutation - nothing happens
        TS_ASSERT_EQUALS(writer.ApplyPermutation(rotation_perm), false);
        //Can't apply a permutation of the wrong length
        short_perm.push_back(0u);
        TS_ASSERT_THROWS_THIS(writer.ApplyPermutation(short_perm), "Permutation doesn't match the expected problem size of 10");

        for (unsigned index = 0; index < (unsigned)number_nodes; index++)
        {
            rotation_perm.push_back((index + 3) % number_nodes); // 3, 4, ... 0, 1, 2
            identity_perm.push_back(index);
        }

        // Make the permutation incorrect
        TS_ASSERT_EQUALS(rotation_perm[0], 3u);
        rotation_perm[0] = 0;

        TS_ASSERT_THROWS_THIS(writer.ApplyPermutation(rotation_perm), "Permutation vector doesn't contain a valid permutation");

        // Correct the mistake imposed above
        rotation_perm[0] = 3;

        TS_ASSERT_EQUALS(writer.ApplyPermutation(identity_perm), false); //Does nothing

        // +++ This is where the permutation is really applied +++
        TS_ASSERT(writer.ApplyPermutation(rotation_perm));

        writer.EndDefineMode();

        // Can't apply permutation after define mode
        TS_ASSERT_THROWS_THIS(writer.ApplyPermutation(rotation_perm), "Cannot define permutation when not in Define mode");

        // However, we can force the permutation to be applied by using the unsafe/extending flag.
        Hdf5DataWriter writer_unsafe_perm(factory, "TestHdf5DataWriter", "hdf5_permuted_unsafe", false);
        writer_unsafe_perm.DefineFixedDimension(number_nodes);
        writer_unsafe_perm.DefineVariable("index", "dimensionless");
        writer_unsafe_perm.DefineVariable("V_m", "millivolts");
        writer_unsafe_perm.DefineVariable("Phi_e", "millivolts");
        writer_unsafe_perm.EndDefineMode();
        writer_unsafe_perm.ApplyPermutation(rotation_perm, /*unsafe*/ true);

        Vec petsc_data_short = factory.CreateVec();
        DistributedVector distributed_vector_short = factory.CreateDistributedVector(petsc_data_short);
        for (DistributedVector::Iterator index = distributed_vector_short.Begin();
             index != distributed_vector_short.End();
             ++index)
        {
            distributed_vector_short[index] = index.Global;
        }
        distributed_vector_short.Restore();
        writer.PutVector(index_id, petsc_data_short);
        writer_unsafe_perm.PutVector(index_id, petsc_data_short);

        Vec petsc_data_long = factory.CreateVec(2);
        DistributedVector distributed_vector_long = factory.CreateDistributedVector(petsc_data_long);
        DistributedVector::Stripe vm_stripe(distributed_vector_long, 0);
        DistributedVector::Stripe phi_e_stripe(distributed_vector_long, 1);
        for (DistributedVector::Iterator index = distributed_vector_long.Begin();
             index != distributed_vector_long.End();
             ++index)
        {
            vm_stripe[index] = 100 + index.Global;
            phi_e_stripe[index] = 1000 + index.Global;
        }
        distributed_vector_long.Restore();
        writer.PutStripedVector(variable_IDs, petsc_data_long);
        writer_unsafe_perm.PutStripedVector(variable_IDs, petsc_data_long);

        writer.Close();
        writer_unsafe_perm.Close();

        PetscTools::Destroy(petsc_data_short);
        PetscTools::Destroy(petsc_data_long);
        TS_ASSERT(CompareFilesViaHdf5DataReaderGlobalNorm("TestHdf5DataWriter", "hdf5_permuted", true,
                                                          "io/test/data", "hdf5_unpermuted", false));
        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_permuted", true,
                                                "io/test/data", "hdf5_permuted", false));
        TS_ASSERT(CompareFilesViaHdf5DataReader("TestHdf5DataWriter", "hdf5_permuted_unsafe", true,
                                                "io/test/data", "hdf5_permuted", false));
    }

    void TestHdf5DataWriterManualChunkSizeAndAlignment()
    {
        std::string folder("TestHdf5DataWriter");
        std::string filename("hdf5_test_manual_chunk_size_and_alignment");

        {
            int number_nodes = 2000;
            DistributedVectorFactory factory(number_nodes);
            Hdf5DataWriter writer(factory, folder, filename, false);

            // Define some dimensions and variables
            writer.DefineUnlimitedDimension("Time", "msec", 101);
            writer.DefineFixedDimension(number_nodes);
            // Odd number of variables for
            writer.DefineVariable("Node", "dimensionless");
            writer.DefineVariable("I_K", "milliamperes");
            writer.DefineVariable("I_Na", "milliamperes");

            /* Set the target chunk size to 8 K (smaller than normal) and the
             * alignment to 16 K.
             * Note: this is a stupid example and just for testing. Because
             * every 8 K chunk will be aligned to 16 K boundaries the file
             * will be about twice the size it needs to be on disk!
             */
            writer.SetTargetChunkSize(0x2000); // 8 K
            writer.SetAlignment(0x4000); // 16 K
            writer.EndDefineMode();

            // Test assertions
            TS_ASSERT_THROWS_THIS(writer.SetTargetChunkSize(123),
                                  "Cannot set chunk target size when not in define mode.");
            TS_ASSERT_THROWS_THIS(writer.SetAlignment(456),
                                  "Cannot set alignment parameter when not in define mode.");

            // Don't bother actually writing anything, that's tested elsewhere
            writer.Close();
        }

        // Open file
        OutputFileHandler file_handler(folder, false);
        FileFinder file = file_handler.FindFile(filename + ".h5");
        hid_t h5_file = H5Fopen(file.GetAbsolutePath().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t dset = H5Dopen(h5_file, "Data", H5P_DEFAULT); // open dataset

        /* Check chunk dimensions are as expected for 8 K chunks with these
         * dataset dimensions. */
        hid_t dcpl = H5Dget_create_plist(dset); // get dataset creation property list
        hsize_t expected_dims[3] = { 17, 20, 3 };
        hsize_t chunk_dims[3];
        H5Pget_chunk(dcpl, 3, chunk_dims);
        for (int i = 0; i < 3; ++i)
        {
            TS_ASSERT_EQUALS(chunk_dims[i], expected_dims[i]);
        }
        H5Pclose(dcpl);

        /*
         * Check the "location" of the datasets (the offset from the start of
         * file, a bit like a pointer to the start of the dataset) to confirm
         * alignment was switched on.
         * With alignment switched off, Data is usually located at 800 B. With
         * alignment = 16 K it is at 64 K. For the Data_Unlimited dataset the
         * numbers are 4935472 B (about 4.7 MB) and 10125312 B (about 9.7 MB),
         * respectively.
         * (These numbers might be machine-dependent!)
         */
        H5O_info_t data_info;
        H5Oget_info(dset, &data_info);
        TS_ASSERT_EQUALS(data_info.addr, 0x10000u); // 64 KB
        H5Dclose(dset);

        dset = H5Dopen(h5_file, "Data_Unlimited", H5P_DEFAULT);
        H5Oget_info(dset, &data_info);
        TS_ASSERT_EQUALS(data_info.addr, 0x9A8000u); // About 9.7 MB

        // Tidy up
        H5Dclose(dset);
        H5Fclose(h5_file);
    }
};

#endif /*TESTHDF5DATAWRITER_HPP_*/
