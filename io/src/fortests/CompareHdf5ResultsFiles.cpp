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

#include "CompareHdf5ResultsFiles.hpp"

#include <iostream>
#include <vector>
#include "DistributedVectorFactory.hpp"
#include "Hdf5DataReader.hpp"
#include "petscvec.h"

bool CompareFilesViaHdf5DataReader(std::string pathname1, std::string filename1, bool makeAbsolute1,
                                   std::string pathname2, std::string filename2, bool makeAbsolute2,
                                   double tol, std::string datasetName)
{
    Hdf5DataReader reader1(pathname1, filename1, makeAbsolute1, datasetName);
    Hdf5DataReader reader2(pathname2, filename2, makeAbsolute2, datasetName);

    unsigned number_nodes1 = reader1.GetNumberOfRows();
    unsigned number_nodes2 = reader2.GetNumberOfRows();
    if (number_nodes1 != number_nodes2)
    {
        std::cout << "Number of nodes " << number_nodes1 << " and " << number_nodes2 << " don't match\n";
        return false;
    }

    // Check the variable names and units
    std::vector<std::string> variable_names1 = reader1.GetVariableNames();
    std::vector<std::string> variable_names2 = reader2.GetVariableNames();
    std::string unlimited_variable1 = reader1.GetUnlimitedDimensionName();
    std::string unlimited_variable2 = reader2.GetUnlimitedDimensionName();
    std::string unlimited_variable_unit1 = reader1.GetUnlimitedDimensionUnit();
    std::string unlimited_variable_unit2 = reader2.GetUnlimitedDimensionUnit();

    unsigned num_vars = variable_names1.size();
    if (num_vars != variable_names2.size())
    {
        std::cout << "Number of variables " << variable_names1.size()
                  << " and " << variable_names2.size() << " don't match\n";
        return false;
    }
    if (unlimited_variable1 != unlimited_variable2)
    {
        std::cout << "Unlimited variable names " << unlimited_variable1 << " and "
                  << unlimited_variable2 << " don't match\n";
        return false;
    }
    if (unlimited_variable_unit1 != unlimited_variable_unit2)
    {
        std::cout << "Unlimited variable units " << unlimited_variable_unit1 << " and "
                  << unlimited_variable_unit2 << " don't match\n";
        return false;
    }
    for (unsigned var = 0; var < num_vars; var++)
    {
        std::string var_name = variable_names1[var];
        if (var_name != variable_names2[var])
        {
            std::cout << "Variable names " << var_name << " and "
                      << variable_names2[var] << " don't match\n";
            return false;
        }
        if (reader1.GetUnit(var_name) != reader2.GetUnit(var_name))
        {
            std::cout << "Units names " << reader1.GetUnit(var_name)
                      << " and " << reader2.GetUnit(var_name) << " don't match\n";
            return false;
        }
    }

    // Check the timestep vectors
    std::vector<double> times1 = reader1.GetUnlimitedDimensionValues();
    std::vector<double> times2 = reader2.GetUnlimitedDimensionValues();

    if (times1.size() != times2.size())
    {
        std::cout << "Number of time steps " << times1.size()
                  << " and " << times2.size() << " don't match\n";
        return false;
    }

    for (unsigned timestep = 0; timestep < times1.size(); timestep++)
    {
        ///\todo remove magic number? (#1884)
        if (fabs(times1[timestep] - times2[timestep]) > 1e-8)
        {
            std::cout << "Time step " << timestep << " doesn't match. First file says " << times1[timestep]
                      << " and second says " << times2[timestep] << std::endl;
            return false;
        }
    }

    bool is_complete1 = reader1.IsDataComplete();
    bool is_complete2 = reader2.IsDataComplete();

    if (is_complete1 != is_complete2)
    {
        std::cout << "One of the readers has incomplete data and the other doesn't\n";
        return false;
    }

    if (is_complete1)
    {
        DistributedVectorFactory factory(number_nodes1);

        Vec data1 = factory.CreateVec();
        Vec data2 = factory.CreateVec();

        for (unsigned timestep = 0; timestep < times1.size(); timestep++)
        {
            for (unsigned var = 0; var < num_vars; var++)
            {
                reader1.GetVariableOverNodes(data1, variable_names1[var], timestep);
                reader2.GetVariableOverNodes(data2, variable_names2[var], timestep);

#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
                double minus_one = -1.0;
                VecAXPY(&minus_one, data2, data1);
#else
                //[note: VecAXPY(y,a,x) computes y = ax+y]
                VecAXPY(data1, -1.0, data2);
#endif

                PetscReal difference_norm;
                VecNorm(data1, NORM_2, &difference_norm);

                if (difference_norm > tol)
                {
                    std::cout << "Time " << times1[timestep] << ": vectors differ in NORM_2 by " << difference_norm << std::endl;
                    return false;
                }
            }
        }
        PetscTools::Destroy(data1);
        PetscTools::Destroy(data2);
    }
    else
    {
        // Incomplete data

        // Check the index vectors
        std::vector<unsigned> indices1 = reader1.GetIncompleteNodeMap();
        std::vector<unsigned> indices2 = reader2.GetIncompleteNodeMap();

        if (indices1.size() != indices2.size())
        {
            std::cout << "Index map sizes " << indices1.size() << " and " << indices2.size() << " don't match\n";
            return false;
        }

        for (unsigned index = 0; index < indices1.size(); index++)
        {
            if (indices1[index] != indices2[index])
            {
                std::cout << "Node indices " << indices1[index] << " and " << indices2[index] << " don't match\n";
                return false;
            }
        }

        // Check all the data
        for (unsigned index = 0; index < indices1.size(); index++)
        {
            unsigned node_index = indices1[index];
            for (unsigned var = 0; var < num_vars; var++)
            {
                std::vector<double> var_over_time1 = reader1.GetVariableOverTime(variable_names1[var], node_index);
                std::vector<double> var_over_time2 = reader2.GetVariableOverTime(variable_names1[var], node_index);
                for (unsigned time_step = 0; time_step < var_over_time1.size(); time_step++)
                {
                    if (fabs(var_over_time1[time_step] - var_over_time2[time_step]) > tol)
                    {
                        std::cout << "Node " << node_index << " at time step " << time_step << " variable " << variable_names1[var]
                                  << " differs (" << var_over_time1[time_step] << " != " << var_over_time2[time_step] << ")" << std::endl;
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

bool CompareFilesViaHdf5DataReaderGlobalNorm(std::string pathname1, std::string filename1, bool makeAbsolute1,
                                             std::string pathname2, std::string filename2, bool makeAbsolute2,
                                             double tol, std::string datasetName)
{
    Hdf5DataReader reader1(pathname1, filename1, makeAbsolute1, datasetName);
    Hdf5DataReader reader2(pathname2, filename2, makeAbsolute2, datasetName);

    unsigned number_nodes1 = reader1.GetNumberOfRows();
    bool is_the_same = true;
    std::vector<std::string> variable_names1 = reader1.GetVariableNames();
    std::vector<std::string> variable_names2 = reader2.GetVariableNames();
    std::vector<double> times1 = reader1.GetUnlimitedDimensionValues();
    unsigned num_vars = variable_names1.size();
    DistributedVectorFactory factory(number_nodes1);

    Vec data1 = factory.CreateVec();
    Vec data2 = factory.CreateVec();

    for (unsigned timestep = 0; timestep < times1.size(); timestep++)
    {
        for (unsigned var = 0; var < num_vars; var++)
        {
            reader1.GetVariableOverNodes(data1, variable_names1[var], timestep);
            reader2.GetVariableOverNodes(data2, variable_names2[var], timestep);

            PetscReal data1_norm;
            PetscReal data2_norm;
            VecNorm(data1, NORM_2, &data1_norm);
            VecNorm(data2, NORM_2, &data2_norm);
            PetscReal norm = fabs(data1_norm - data2_norm);
            if (norm > tol)
            {
                is_the_same = false;
                std::cout << "Vectors differ in global NORM_2 by " << norm << std::endl;
            }
        }
    }

    PetscTools::Destroy(data1);
    PetscTools::Destroy(data2);

    return is_the_same;
}
