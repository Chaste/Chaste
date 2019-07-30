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

#ifndef TESTHDF5READER_HPP_
#define TESTHDF5READER_HPP_

#include <cxxtest/TestSuite.h>

#include "Hdf5DataWriter.hpp"
#include "Hdf5DataReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "DistributedVectorFactory.hpp"


/*
 * Operator function to be called by H5Literate in TestListingDatasetsInAnHdf5File.
 */
herr_t op_func (hid_t loc_id,
                const char *name,
                const H5L_info_t *info,
                void *operator_data);


class TestHdf5DataReader : public CxxTest::TestSuite
{
private:

    std::string h5file_name;
    #define DATASETNAME "IntArray"

    void WriteDataTestSimpleReadDirectlyWithHdf5()
    {
        int const NX = 5;                      /* dataset dimensions */
        int const NY = 6;
        int const RANK = 2;

        h5file_name = OutputFileHandler::GetChasteTestOutputDirectory() + "SDS.h5";

        hid_t       file, dataset;         /* file and dataset handles */
        hid_t       datatype, dataspace;   /* handles */
        hsize_t     dimsf[2];              /* dataset dimensions */
        herr_t      status;
        int         data[NX][NY];          /* data to write */
        int         i, j;

        /*
         * Data  and output buffer initialization.
         */
        for (j = 0; j < NX; j++) {
        for (i = 0; i < NY; i++)
            data[j][i] = i + j;
        }
        /*
         * 0 1 2 3 4 5
         * 1 2 3 4 5 6
         * 2 3 4 5 6 7
         * 3 4 5 6 7 8
         * 4 5 6 7 8 9
         */

        /*
         * Create a new file using H5F_ACC_TRUNC access,
         * default file creation properties, and default file
         * access properties.
         */
        file = H5Fcreate(h5file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        /*
         * Describe the size of the array and create the data space for fixed
         * size dataset.
         */
        dimsf[0] = NX;
        dimsf[1] = NY;
        dataspace = H5Screate_simple(RANK, dimsf, NULL);

        /*
         * Define datatype for the data in the file.
         * We will store little endian INT numbers.
         */
        datatype = H5Tcopy(H5T_NATIVE_INT);
        status = H5Tset_order(datatype, H5T_ORDER_LE);

        /*
         * Create a new dataset within the file using defined dataspace and
         * datatype and default dataset creation properties.
         */
        dataset = H5Dcreate(file, DATASETNAME, datatype, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        /*
         * Write the data to the dataset using default transfer properties.
         */
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        UNUSED_OPT(status);
        assert(status==0);

        /*
         * Close/release resources.
         */
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Fclose(file);
        PetscTools::Barrier("WriteDataTestSimpleReadDirectlyWithHdf5");
    }

public:

    void TestSimpleReadDirectlyWithHdf5()
    {
        int const NX_SUB = 3;           /* hyperslab dimensions */
        int const NY_SUB = 4;
        int const NX = 7;           /* output buffer dimensions */
        int const NY = 7;
        int const NZ = 3;
        //int const RANK      =   2;
        int const RANK_OUT  =   3;

        hid_t       file, dataset;         /* handles */
        hid_t       datatype, dataspace;
        hid_t       memspace;
        H5T_class_t datatype_class;        /* datatype class */
        H5T_order_t order;                 /* data order */
        size_t      size;                  /* size of the data element stored in file */
        hsize_t     dimsm[3];              /* memory space dimensions */
        hsize_t     dims_out[2];           /* dataset dimensions */
        herr_t      status;

        int         data_out[NX][NY][NZ ]; /* output buffer */

        hsize_t      count[2];              /* size of the hyperslab in the file */
        hsize_t      offset[2];             /* hyperslab offset in the file */
        hsize_t      count_out[3];          /* size of the hyperslab in memory */
        hsize_t      offset_out[3];         /* hyperslab offset in memory */
        int          i, j, k, status_n, rank;

        // Create the file it's gonna be read
        WriteDataTestSimpleReadDirectlyWithHdf5();

        for (j = 0; j < NX; j++)
        {
            for (i = 0; i < NY; i++)
            {
                for (k = 0; k < NZ; k++)
                    data_out[j][i][k] = 0;
            }
        }

        // Open the file and the dataset
        file = H5Fopen(h5file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        dataset = H5Dopen(file, DATASETNAME, H5P_DEFAULT);

        /*
         * Get datatype and dataspace handles and then query
         * dataset class, order, size, rank and dimensions.
         */
        datatype  = H5Dget_type(dataset);     /* datatype handle */
        datatype_class     = H5Tget_class(datatype);
        TS_ASSERT(datatype_class == H5T_INTEGER);
        //if (datatype_class == H5T_INTEGER) printf("Data set has INTEGER type \n");
        order     = H5Tget_order(datatype);
        TS_ASSERT(order == H5T_ORDER_LE);
        //if (order == H5T_ORDER_LE) printf("Little endian order \n");

        size  = H5Tget_size(datatype);
        TS_ASSERT_EQUALS(size,4u);
        //printf(" Data size is %d \n", size);

        dataspace = H5Dget_space(dataset);    /* dataspace handle */
        rank      = H5Sget_simple_extent_ndims(dataspace);
        status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
        TS_ASSERT_EQUALS(rank,2);
        TS_ASSERT_EQUALS(status_n,2);
        //printf("rank %d, dimensions %lu x %lu \n", rank,
        //   (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));

        // Define hyperslab in the dataset
        offset[0] = 1;
        offset[1] = 2;
        count[0]  = NX_SUB;
        count[1]  = NY_SUB;
        status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,
                     count, NULL);
        TS_ASSERT_EQUALS(status, 0);

        // Define the memory dataspace
        dimsm[0] = NX;
        dimsm[1] = NY;
        dimsm[2] = NZ;
        memspace = H5Screate_simple(RANK_OUT,dimsm,NULL);

        // Define memory hyperslab
        offset_out[0] = 3;
        offset_out[1] = 0;
        offset_out[2] = 0;
        count_out[0]  = NX_SUB;
        count_out[1]  = NY_SUB;
        count_out[2]  = 1;
        status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL,
                     count_out, NULL);

        /*
         * Read data from hyperslab in the file into the hyperslab in
         * memory and display.
         */
        status = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace,
                 H5P_DEFAULT, data_out);

        /*
         * 0 0 0 0 0 0 0
         * 0 0 0 0 0 0 0
         * 0 0 0 0 0 0 0
         * 3 4 5 6 0 0 0
         * 4 5 6 7 0 0 0
         * 5 6 7 8 0 0 0
         * 0 0 0 0 0 0 0
         */
        for (int row=0; row<NX; row++)
        {
            for (int column=0; column<NY; column++)
            {
                if (row>2 && row<6 && column<4)
                {
                    TS_ASSERT_EQUALS(data_out[row][column][0], row+column);
                }
                else
                {
                    TS_ASSERT_EQUALS(data_out[row][column][0], 0);
                }
            }
        }

        // Close/release resources
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Sclose(memspace);
        H5Fclose(file);
    }

private:

    static const unsigned NUMBER_NODES = 100;

    void WriteMultiStepData()
    {
        DistributedVectorFactory factory(NUMBER_NODES);

        Hdf5DataWriter writer(factory, "hdf5_reader", "hdf5_test_complete_format", false);
        writer.DefineFixedDimension(NUMBER_NODES);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        int ik_id = writer.DefineVariable("I_K", "milliamperes");
        int ina_id = writer.DefineVariable("I_Na", "milliamperes");
        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        Vec petsc_data_1=factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2=factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        Vec petsc_data_3=factory.CreateVec();
        DistributedVector distributed_vector_3 = factory.CreateDistributedVector(petsc_data_3);

        for (unsigned time_step=0; time_step<10; time_step++)
        {
            // Write some values
            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index!= distributed_vector_1.End();
                 ++index)
            {
                distributed_vector_1[index] =  index.Global;
                distributed_vector_2[index] =  time_step*1000 + 100 + index.Global;
                distributed_vector_3[index] =  time_step*1000 + 200 + index.Global;
            }
            distributed_vector_1.Restore();
            distributed_vector_2.Restore();
            distributed_vector_3.Restore();

            // Write the vector
            writer.PutVector(node_id, petsc_data_1);
            writer.PutVector(ik_id, petsc_data_2);
            writer.PutVector(ina_id, petsc_data_3);
            writer.PutUnlimitedVariable(time_step);

            // Don't advance in the last iteration
            if (time_step < 9)
            {
                writer.AdvanceAlongUnlimitedDimension();
            }
        }
        PetscTools::Destroy(petsc_data_1);
        PetscTools::Destroy(petsc_data_2);
        PetscTools::Destroy(petsc_data_3);
        writer.Close();
    }

public:

    void TestMultiStepReader()
    {
        WriteMultiStepData();

        Hdf5DataReader reader("hdf5_reader", "hdf5_test_complete_format");

        std::vector<std::string> variable_names=reader.GetVariableNames();
        TS_ASSERT_EQUALS(variable_names.size(), 3U);
        TS_ASSERT_EQUALS(variable_names[0], "Node");
        TS_ASSERT_EQUALS(reader.GetUnit("Node"), "dimensionless");
        TS_ASSERT_EQUALS(variable_names[1], "I_K");
        TS_ASSERT_EQUALS(reader.GetUnit("I_K"), "milliamperes");
        TS_ASSERT_EQUALS(variable_names[2], "I_Na");
        TS_ASSERT_EQUALS(reader.GetUnit("I_Na"), "milliamperes");

        for (unsigned node_index=0; node_index<NUMBER_NODES; node_index++)
        {
            std::vector<double> node_values = reader.GetVariableOverTime("Node", node_index);
            std::vector<double> i_k_values = reader.GetVariableOverTime("I_K", node_index);
            std::vector<double> i_na_values = reader.GetVariableOverTime("I_Na", node_index);

            TS_ASSERT_EQUALS(node_values.size(), 10u);
            TS_ASSERT_EQUALS(i_k_values.size(), 10u);
            TS_ASSERT_EQUALS(i_na_values.size(), 10u);

            for (unsigned i=0; i<node_values.size(); i++)
            {
                TS_ASSERT_DELTA( node_values[i], node_index, 1e-9);
                TS_ASSERT_DELTA( i_k_values[i], i*1000 + 100 + node_index, 1e-9);
                TS_ASSERT_DELTA( i_na_values[i], i*1000 + 200 + node_index, 1e-9);
            }
        }

        std::vector<double> i_k_values = reader.GetVariableOverTime("I_K", 15);
        std::vector<double> i_k_values_over_multiple = reader.GetVariableOverTimeOverMultipleNodes("I_K", 10, 19)[5];

        // If 19 is 'one past the end' of the requested nodes, shouldn't this be of size 9 ??
        TS_ASSERT_EQUALS(i_k_values_over_multiple.size(),10u);

        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTimeOverMultipleNodes("I_K", 0, NUMBER_NODES+5),
                              "The dataset 'Data' doesn't contain info for node 104");

        TS_ASSERT_EQUALS(i_k_values.size(), i_k_values_over_multiple.size());

        for (unsigned index=0; index< i_k_values.size(); index++)
        {
            TS_ASSERT_EQUALS(i_k_values[index], i_k_values_over_multiple[index]);
        }

        unsigned NUMBER_NODES=100;
        DistributedVectorFactory factory(NUMBER_NODES);

        Vec petsc_data_1=factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2=factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        Vec petsc_data_3=factory.CreateVec();
        DistributedVector distributed_vector_3 = factory.CreateDistributedVector(petsc_data_3);

        TS_ASSERT_EQUALS(reader.GetNumberOfRows(), NUMBER_NODES);
        for (unsigned time_step=0; time_step<10; time_step++)
        {
            reader.GetVariableOverNodes(petsc_data_1, "Node", time_step);
            reader.GetVariableOverNodes(petsc_data_2, "I_K", time_step);
            reader.GetVariableOverNodes(petsc_data_3, "I_Na", time_step);

            distributed_vector_1.Restore();
            distributed_vector_2.Restore();
            distributed_vector_3.Restore();

            // Check values
            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index!= distributed_vector_1.End();
                 ++index)
            {
                TS_ASSERT_EQUALS(distributed_vector_1[index], index.Global);
                TS_ASSERT_EQUALS(distributed_vector_2[index], time_step*1000 + 100 + index.Global);
                TS_ASSERT_EQUALS(distributed_vector_3[index], time_step*1000 + 200 + index.Global);
            }
        }

        std::vector<double> unlimited_values = reader.GetUnlimitedDimensionValues();

        for (unsigned i=0; i< unlimited_values.size(); i++)
        {
            TS_ASSERT_EQUALS(unlimited_values[i], i);
        }

        PetscTools::Destroy(petsc_data_1);
        PetscTools::Destroy(petsc_data_2);
        PetscTools::Destroy(petsc_data_3);
        reader.Close();
    }

    void TestNonMultiStepExceptions()
    {
        DistributedVectorFactory factory(NUMBER_NODES);

        Hdf5DataWriter writer(factory, "hdf5_reader", "hdf5_test_overtime_exceptions", false);
        writer.DefineFixedDimension(NUMBER_NODES);

        writer.DefineVariable("Node", "dimensionless");
        writer.DefineVariable("I_K", "milliamperes");
        writer.DefineVariable("I_Na", "milliamperes");

        writer.EndDefineMode();
        writer.Close();

        Hdf5DataReader reader("hdf5_reader", "hdf5_test_overtime_exceptions");

        // Unlimited dimension has a default return value
        std::vector<double> times = reader.GetUnlimitedDimensionValues();
        TS_ASSERT_EQUALS(times.size(),1U);
        TS_ASSERT_EQUALS(times[0],0.0);

        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTime("Node", 99/*node*/),
                "The dataset 'Data' does not contain time dependent data");

        Vec data = factory.CreateVec();
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverNodes(data, "Node", 1/*timestep*/),
                "The dataset 'Data' does not contain time dependent data");

        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTimeOverMultipleNodes("Node",0,10),
                "The dataset 'Data' does not contain time dependent data");

        PetscTools::Destroy(data);
        reader.Close();

        // Missing dir
        FileFinder absent_dir("absent_dir", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS(Hdf5DataReader(absent_dir, "base"), "Directory does not exist: ");

        H5E_BEGIN_TRY //Supress HDF5 error in this test
        {
            TS_ASSERT_THROWS_CONTAINS(Hdf5DataReader("hdf5_reader", "hdf5_test_overtime_exceptions", true, "Postprocessing"),
                    "but could not find the dataset 'Postprocessing',");
        }
        H5E_END_TRY;
    }

    void TestMultiStepExceptions()
    {
        DistributedVectorFactory factory(NUMBER_NODES);

        Hdf5DataWriter writer(factory, "hdf5_reader", "hdf5_test_overtime_exceptions", false);
        DistributedVectorFactory vec_factor(NUMBER_NODES);
        writer.DefineFixedDimension(NUMBER_NODES);

        writer.DefineVariable("Node", "dimensionless");
        writer.DefineVariable("I_K", "milliamperes");
        writer.DefineVariable("I_Na", "milliamperes");
        writer.DefineUnlimitedDimension("Time", "msec");

        writer.EndDefineMode();

        writer.AdvanceAlongUnlimitedDimension();

        writer.Close();

        // Test some exceptions
        {
            // Check for a missing file exception
            TS_ASSERT_THROWS_CONTAINS(Hdf5DataReader reader2("hdf5_reader", "hdf5_wrong_name"),
                                      "as it does not exist");

            // Check for a wrong file format exception
            OutputFileHandler file_handler("hdf5_reader",false);
            out_stream p_wrong_format = file_handler.OpenOutputFile("hdf5_wrong_format.h5");
            p_wrong_format->close();
            H5E_BEGIN_TRY //Supress HDF5 error in this test
            {
                TS_ASSERT_THROWS_CONTAINS(Hdf5DataReader reader2("hdf5_reader", "hdf5_wrong_format"),
                                          "H5Fopen error code");
            }
            H5E_END_TRY;
        }

        Hdf5DataReader reader("hdf5_reader", "hdf5_test_overtime_exceptions");

        TS_ASSERT_THROWS_NOTHING(reader.GetVariableOverTime("Node", 99/*node*/));
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTime("WrongName", 99/*node*/),
                "The dataset 'Data' doesn\'t contain data for variable WrongName");
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTime("Node", 100/*node*/),
                "The dataset 'Data' doesn\'t contain info of node 100");
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTimeOverMultipleNodes("WrongName", 90, 99),
                "The dataset 'Data' doesn\'t contain data for variable WrongName");

        Vec data = factory.CreateVec();
        reader.GetVariableOverNodes(data, "Node", 0/*timestep*/);
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverNodes(data, "WrongName"),
                "The dataset 'Data' does not contain data for variable WrongName"); //Wrong name
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverNodes(data, "I_K", 1/*timestep*/),
                "The dataset 'Data' does not contain data for timestep number 1"); //Time step doesn't exist

        DistributedVectorFactory factory2(NUMBER_NODES+1);
        Vec data_too_big = factory2.CreateVec();
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverNodes(data_too_big, "Node", 0/*timestep*/),
                "Could not read data because Vec is the wrong size"); //Data too big

        PetscTools::Destroy(data);
        PetscTools::Destroy(data_too_big);
        reader.Close();
    }

    void TestIncompleteData()
    {
        DistributedVectorFactory factory(NUMBER_NODES);

        Hdf5DataReader reader("io/test/data","hdf5_test_full_format_incomplete", false);

        std::vector<std::string> variable_names = reader.GetVariableNames();
        TS_ASSERT_EQUALS(variable_names.size(), 3u);
        TS_ASSERT_EQUALS(variable_names[0], "Node");
        TS_ASSERT_EQUALS(reader.GetUnit("Node"), "dimensionless");
        TS_ASSERT_EQUALS(variable_names[1], "I_K");
        TS_ASSERT_EQUALS(reader.GetUnit("I_K"), "milliamperes");
        TS_ASSERT_EQUALS(variable_names[2], "I_Na");
        TS_ASSERT_EQUALS(reader.GetUnit("I_Na"), "milliamperes");

        // Can't read into a PETSc Vec
        Vec data = factory.CreateVec();
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverNodes(data, "Node", 1/*timestep*/),
                "You can only get a vector for complete data");
        PetscTools::Destroy(data);

        std::vector<unsigned> nodes=reader.GetIncompleteNodeMap();
        TS_ASSERT_EQUALS(nodes.size(), 3U);
        TS_ASSERT_EQUALS(nodes[0], 21U);
        TS_ASSERT_EQUALS(nodes[1], 47U);
        TS_ASSERT_EQUALS(nodes[2], 60U);

        // Can read one of the nodes that was written
        std::vector<double> twenty_one = reader.GetVariableOverTime("Node", 21);
        TS_ASSERT_EQUALS(twenty_one.size(), 10u);
        for (unsigned i=0; i<twenty_one.size(); i++)
        {
            TS_ASSERT_EQUALS(twenty_one[i], 21u);
        }
        // Can read more of the data
        std::vector<double> na_47 = reader.GetVariableOverTime("I_Na", 47);
        TS_ASSERT_EQUALS(na_47.size(), 10u);
        for (unsigned i=0; i<na_47.size(); i++)
        {
            TS_ASSERT_EQUALS(na_47[i], i*1000u + 200u + 47u);
        }

        // Data not included
        TS_ASSERT_THROWS_THIS(reader.GetVariableOverTime("Node", 22),
                              "The incomplete dataset 'Data' does not contain info of node 22");

        // another exception
        TS_ASSERT_THROWS_CONTAINS(reader.GetVariableOverTimeOverMultipleNodes("I_Na", 0, 1),
                                  "GetVariableOverTimeOverMultipleNodes() cannot be called using incomplete data sets");
    }

    void TestReadingExtraData()
    {
        // In this test we read data from a file in the source
        // (which is re-created and compared with in TestHdf5DataWriter.hpp)

        { // Original / Normal data in a datastructure called "Data"
            Hdf5DataReader reader("io/test/data","hdf5_test_adding_variables", false);
            std::vector<std::string> variable_names = reader.GetVariableNames();
            TS_ASSERT_EQUALS(variable_names[0], "Node");
            TS_ASSERT_EQUALS(variable_names[1], "I_K");
        }

        { // Added data in a datastructure called "Extra stuff"
            Hdf5DataReader reader("io/test/data","hdf5_test_adding_variables", false, "Extra stuff");
            std::vector<std::string> variable_names = reader.GetVariableNames();
            TS_ASSERT_EQUALS(variable_names[0], "Phase");
            TS_ASSERT_EQUALS(variable_names[1], "Plasma");

            TS_ASSERT_EQUALS(reader.GetUnlimitedDimensionName(), "Distance");
            TS_ASSERT_EQUALS(reader.GetUnlimitedDimensionUnit(), "LightYears");

            DistributedVectorFactory factory(reader.GetNumberOfRows());
            Vec data = factory.CreateVec();
            reader.GetVariableOverNodes(data, "Plasma", 1/*timestep*/);
            DistributedVector distributed_vector_1 = factory.CreateDistributedVector(data);

            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                                 index!= distributed_vector_1.End();
                                 ++index)
            {
                // These magic numbers come from TestHdf5DataWriter - where they are written into the HDF5 file.
                TS_ASSERT_DELTA(distributed_vector_1[index], 1*1000 + 100 + index.Global, 1e-9);
            }

            PetscTools::Destroy(data);
        }
    }

    void TestListingDatasetsInAnHdf5File()
    {
        /*
         * Open file.
         */
        hid_t file = H5Fopen ("mesh/test/data/vtk_extending/SimulationResults.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

        /*
         * Make an object for us to store the dataset names in.
         */
        std::vector<std::string> dataset_names;

        /*
         * Begin iteration.
         */
        herr_t status = H5Literate(file, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, op_func, &dataset_names);

        /*
         * Close and release resources.
         */
        TS_ASSERT_EQUALS(status, 0);
        status = H5Fclose (file);
        TS_ASSERT_EQUALS(status, 0);

        /*
         * Check that all the datasets have been listed correctly.
         */
        TS_ASSERT_EQUALS(dataset_names.size(), 4u);
        TS_ASSERT_EQUALS(dataset_names[0], "Apd_60_minus_30_Map");
        TS_ASSERT_EQUALS(dataset_names[1], "Apd_60_minus_30_Map_Unlimited" );
        TS_ASSERT_EQUALS(dataset_names[2], "Data");
        TS_ASSERT_EQUALS(dataset_names[3], "Data_Unlimited");
    }
};

/************************************************************

  Operator function.  Prints the name and type of the object
  being examined.

 ************************************************************/
herr_t op_func (hid_t loc_id, const char *name,
                const H5L_info_t *info,
                void *operator_data)
{
    std::vector<std::string>* p_dataset_names = static_cast<std::vector< std::string > * >(operator_data);

    /*
     * Get type of the object and display its name and type.
     * The name of the object is passed to this function by
     * the Library.
     */
    H5O_info_t infobuf;
    H5Oget_info_by_name (loc_id, name, &infobuf, H5P_DEFAULT);
    switch (infobuf.type)
    {
//          case H5O_TYPE_GROUP:
//              printf ("  Group: %s\n", name);
//              break;
        case H5O_TYPE_DATASET:
            p_dataset_names->push_back(name);
            break;
//          case H5O_TYPE_NAMED_DATATYPE:
//              printf ("  Datatype: %s\n", name);
//              break;
        default:
            EXCEPTION("File includes HDF5 object that it shouldn't.");
    }

    return 0;
}

#endif /*TESTHDF5READER_HPP_*/
