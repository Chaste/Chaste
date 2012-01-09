/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef HDF5DATAREADER_HPP_
#define HDF5DATAREADER_HPP_

#ifndef H5_USE_16_API
#define H5_USE_16_API 1
#endif

#include <hdf5.h>
#include <petscvec.h>
#include <string>
#include <vector>
#include <map>

#include "FileFinder.hpp"

const unsigned MAX_STRING_SIZE = 100; /// \todo: magic number

/**
 * A concrete HDF5 data reader class.
 */
class Hdf5DataReader
{
private:

    static const unsigned MAX_DATASET_RANK = 3;             /**< Defined in HDF5 writer too. \todo: define it once */

    std::string mDirectory;                                 /**< Directory output files will be stored in (absolute path). */
    std::string mBaseName;                                  /**< The base name for the output data files. */

    hid_t mFileId;                                          /**< The data file ID. */

    hid_t mVariablesDatasetId;                              /**< The variables data set ID. */
    unsigned mVariablesDatasetRank;                         /**< The rank of the variables data set. */
    hsize_t mVariablesDatasetSizes[MAX_DATASET_RANK];       /**< The sizes of each variable data set. */

    bool mIsUnlimitedDimensionSet;                          /**< Is the unlimited dimension set */
    hid_t mTimeDatasetId;                                   /**< The time data set ID. */
    hsize_t mNumberTimesteps;                               /**< The number of time steps recorded in the data file. */

    std::vector<std::string> mVariableNames;                /**< The variable names. */
    std::map<std::string, unsigned> mVariableToColumnIndex; /**< Map between variable names and data column numbers. */
    std::map<std::string, std::string> mVariableToUnit;     /**< Map between variable names and variable units. */

    bool mIsDataComplete;                                   /**< Whether the data file is complete. */
    std::vector<unsigned> mIncompleteNodeIndices;           /**< Vector of node indices for which the data file does not contain data. */

    bool mClosed;                                           /**< Whether we've already closed the file. */

    /**
     * Contains functionality common to both constructors.
     *
     * @param rDirectory  The directory the files are stored in
     * @param rBaseName  The base name of the files to read (i.e. without the extensions)
     */
    void CommonConstructor(const FileFinder& rDirectory,
                           const std::string& rBaseName);

public:

    /**
     * Read data from the given files into memory.
     *
     * @param rDirectory  The directory the files are stored in
     * @param rBaseName  The base name of the files to read (i.e. without the extensions)
     * @param makeAbsolute  Whether to convert directory to an absolute path using the
     *                      OutputFileHandler (defaults to true)
     */
    Hdf5DataReader(const std::string& rDirectory,
                   const std::string& rBaseName,
                   bool makeAbsolute=true);

    /**
     * Alternative constructor taking a FileFinder to specify the directory.
     *
     * @param rDirectory  The directory the files are stored in
     * @param rBaseName  The base name of the files to read (i.e. without the extensions)
     */
    Hdf5DataReader(const FileFinder& rDirectory,
                   const std::string& rBaseName);

    /**
     * Get the values of a given variable at each time step at a given node.
     *
     * @param rVariableName  name of a variable in the data file
     * @param nodeIndex the index of the node for which the data is obtained
     */
    std::vector<double> GetVariableOverTime(const std::string& rVariableName,
                                            unsigned nodeIndex);

    /**
     * Get the values of a given variable at each time step over multiple nodes.
     *
     * @param rVariableName  name of a variable in the data file
     * @param lowerIndex the index of the lower node for which the data is obtained
     * @param upperIndex one past the index of the upper node for which the data is obtained
     */
    std::vector<std::vector<double> > GetVariableOverTimeOverMultipleNodes(const std::string& rVariableName,
                                                                           unsigned lowerIndex,
                                                                           unsigned upperIndex);

    /**
     * Get the values of a given variable at each node at a given time step.
     *
     * @param data  PETSc vec to hold the data
     * @param rVariableName  name of a variable in the data file
     * @param timestep the time step for which the data is obtained (defaults to 0)
     */
    void GetVariableOverNodes(Vec data, const std::string& rVariableName, unsigned timestep=0);

    /**
     * Get the unlimited dimension values.
     */
    std::vector<double> GetUnlimitedDimensionValues();

    /**
     * Get the number of rows in the data file.
     */
    unsigned GetNumberOfRows();

    /**
     * Get the variable names.
     */
    std::vector<std::string> GetVariableNames();

    /**
     * Get the units in which a given variable is measured.
     *
     * @param rVariableName  name of a variable in the data file
     */
    std::string GetUnit(const std::string& rVariableName);

    /**
     * Get method for mIsDataComplete.
     */
    bool IsDataComplete();

    /**
     * Get method for mIncompleteNodeIndices.
     */
    std::vector<unsigned> GetIncompleteNodeMap();

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
