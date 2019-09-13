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

#ifndef ABSTRACTHDF5ACCESS_HPP_
#define ABSTRACTHDF5ACCESS_HPP_

#include <hdf5.h>
#include <string>
#include "FileFinder.hpp"

const unsigned MAX_STRING_SIZE = 100; /// \todo: magic number

/**
 * An abstract class to get common code for reading and writing HDF5 files into one place.
 *
 * It doesn't do very much, but provides common member variables.
 */
class AbstractHdf5Access
{

protected:
    std::string mBaseName;                        /**< The base name for the data files. */
    std::string mDatasetName;                     /**< The base name for the dataset we are reading/writing. */
    FileFinder mDirectory;                        /**< Directory HDF5 file will be, or is, stored in. */
    static const unsigned DATASET_DIMS = 3;       /**< The dimensions of each dataset (variable, node, time). */

    bool mIsDataComplete;                         /**< Whether the data file contains entries for all nodes. */
    std::string mUnlimitedDimensionName;          /**< The name of the unlimited dimension. */
    std::string mUnlimitedDimensionUnit;          /**< The physical units of the unlimited dimension. */
    bool mIsUnlimitedDimensionSet;                /**< Is the unlimited dimension set up at the moment?*/
    std::vector<unsigned> mIncompleteNodeIndices; /**< Vector of node indices for which the data file does not contain data. */

    hid_t mFileId;                                /**< The data file ID. */
    hid_t mUnlimitedDatasetId;                    /**< The dataset ID for the unlimited (independent) variable - usually time. */
    hid_t mVariablesDatasetId;                    /**< The dataset ID for the dependent variables. */
    hsize_t mDatasetDims[DATASET_DIMS];           /**< The sizes of each dimension of the dataset (variable, node, time). */

    /**
     * Check for the existence of a dataset in an HDF5 file.
     *
     * @param rDatasetName  the name of the dataset.
     * @return whether or not the dataset already exists in the file.
     */
    bool DoesDatasetExist(const std::string& rDatasetName);

    /**
     * This method sets #mUnlimitedDatasetId, the unlimited dataset ID. It does this
     * by looking for a given unlimited dataset name in the file.
     *
     * Either one of the new format of [#mDatasetName + "_Unlimited"],
     * or one of the old format of ["Time"].
     */
    void SetUnlimitedDatasetId();

    /**
     * Sets the raw dataset chunk cache for our main dataset (#mVariablesDatasetId).
     * The default in HDF5 is 1 MB, which is too small for many problems, so we've
     * bumped it up to 128 M.
     *
     * Note: this cache is not currently used with the parallel (MPIO and MPIPOSIX)
     * drivers in read/write mode (as we have in #Hdf5DataWriter). However it should
     * be used by the #Hdf5DataReader (which uses the default driver).
     *
     * If system memory is at a premium you may find you need to reduce the size of
     * this cache by reducing max_bytes_in_cache.
     */
    void SetMainDatasetRawChunkCache();

public:
    /**
     * Constructor for directory given as string
     *
     * @param rDirectory  the directory in which to read/write the data to file
     * @param rBaseName  The base name of the HDF5 file (with no extension)
     * @param rDatasetName  The dataset name - default is "Data" for voltage and extracellular potential.
     * @param makeAbsolute Whether the h5 file should be treated as relative to Chaste test output, otherwise treated as CWD or absolute.
     */
    AbstractHdf5Access(const std::string& rDirectory,
                       const std::string& rBaseName,
                       const std::string& rDatasetName,
                       bool makeAbsolute = true);

    /**
     * Constructor for directory given as FileFinder
     *
     * @param rDirectory  the directory in which to read/write the data to file
     * @param rBaseName  The base name of the HDF5 file (with no extension)
     * @param rDatasetName  The dataset name - default is "Data" for voltage and extracellular potential.
     */
    AbstractHdf5Access(const FileFinder& rDirectory,
                       const std::string& rBaseName,
                       const std::string& rDatasetName);

    /**
     * Destructor.
     */
    virtual ~AbstractHdf5Access();

    /**
     * Get method for mIsDataComplete.
     * If this is true the data file contains entries for all nodes.
     * If this is false the data file only contains entries for some nodes.
     *
     * @return #mIsDataComplete.
     */
    bool IsDataComplete();

    /**
     * @return #mIncompleteNodeIndices.
     * Vector of node indices for which the data file does not contain data.
     */
    std::vector<unsigned> GetIncompleteNodeMap();

    /**
     * @return the name of the Unlimited dimension (usually "Time").
     */
    std::string GetUnlimitedDimensionName();

    /**
     * @return the unit of the Unlimited dimension (usually "msec").
     */
    std::string GetUnlimitedDimensionUnit();
};

#endif // ABSTRACTHDF5ACCESS_HPP_
