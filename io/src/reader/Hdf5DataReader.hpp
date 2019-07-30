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

#ifndef HDF5DATAREADER_HPP_
#define HDF5DATAREADER_HPP_

#include <petscvec.h>
#include <vector>
#include <map>

#include "AbstractHdf5Access.hpp"

/**
 * A concrete HDF5 data reader class.
 *
 * \todo Add to documentation whether we can call this in serial/parallel, whether we can call it when file is open already etc...
 */
class Hdf5DataReader : public AbstractHdf5Access
{
private:

    unsigned mVariablesDatasetRank;                         /**< The rank of the variables data set. */

    hsize_t mNumberTimesteps;                               /**< The number of time steps recorded in the data file. */

    std::vector<std::string> mVariableNames;                /**< The variable names. */
    std::map<std::string, unsigned> mVariableToColumnIndex; /**< Map between variable names and data column numbers. */
    std::map<std::string, std::string> mVariableToUnit;     /**< Map between variable names and variable units. */

    bool mClosed;                                           /**< Whether we've already closed the file. */

    /**
     * Contains functionality common to both constructors.
     */
    void CommonConstructor();

public:

    /**
     * Read data from the given files into memory.
     *
     * @param rDirectory  The directory the files are stored in
     * @param rBaseName  The base name of the files to read (i.e. without the extensions)
     * @param makeAbsolute  Whether the h5 file should be treated as relative to Chaste test output,
     *                      and converted to absolute, the file is otherwise treated as relative to current working directory.
     * @param datasetName  The name of the HDF5 dataset to read, defaults to "Data".
     */
    Hdf5DataReader(const std::string& rDirectory,
                   const std::string& rBaseName,
                   bool makeAbsolute=true,
                   std::string datasetName="Data");

    /**
     * Alternative constructor taking a FileFinder to specify the directory.
     *
     * @param rDirectory  The directory the files are stored in
     * @param rBaseName  The base name of the files to read (i.e. without the extensions)
     * @param datasetName The name of the HDF5 dataset to read, defaults to "Data".
     */
    Hdf5DataReader(const FileFinder& rDirectory,
                   const std::string& rBaseName,
                   std::string datasetName="Data");

    /**
     * @return the values of a given variable at each time step at a given node.
     *
     * @param rVariableName  name of a variable in the data file
     * @param nodeIndex the index of the node for which the data is obtained
     */
    std::vector<double> GetVariableOverTime(const std::string& rVariableName,
                                            unsigned nodeIndex);

    /**
     * @return the values of a given variable at each time step over multiple nodes.
     *
     * @param rVariableName  name of a variable in the data file
     * @param lowerIndex the index of the lower node for which the data is obtained
     * @param upperIndex one past the index of the upper node for which the data is obtained
     */
    std::vector<std::vector<double> > GetVariableOverTimeOverMultipleNodes(const std::string& rVariableName,
                                                                           unsigned lowerIndex,
                                                                           unsigned upperIndex);

    /**
     * @return the values of a given variable at each node at a given time step.
     *
     * @param data  PETSc vec to hold the data
     * @param rVariableName  name of a variable in the data file
     * @param timestep the time step for which the data is obtained (defaults to 0)
     */
    void GetVariableOverNodes(Vec data, const std::string& rVariableName, unsigned timestep=0);

    /**
     * @return the unlimited dimension values.
     */
    std::vector<double> GetUnlimitedDimensionValues();

    /**
     * @return the number of rows in the data file.
     */
    unsigned GetNumberOfRows();

    /**
     * @return the variable names.
     */
    std::vector<std::string> GetVariableNames();

    /**
     * @return the units in which a given variable is measured.
     *
     * @param rVariableName  name of a variable in the data file
     */
    std::string GetUnit(const std::string& rVariableName);

    /**
     * Close any open files.
     */
    void Close();

    /**
     * Destructor just calls Close.
     */
    ~Hdf5DataReader();
};

#endif /*HDF5DATAREADER_HPP_*/
