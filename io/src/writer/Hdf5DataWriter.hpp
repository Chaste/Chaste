/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef HDF5DATAWRITER_HPP_
#define HDF5DATAWRITER_HPP_


// Not sure why this seems to need to be included here as well as in AbstractHdf5Access.hpp,
// but it does on some of the build machines.
#ifndef H5_USE_16_API
#define H5_USE_16_API 1
#endif

#include <vector>

#include "AbstractHdf5Access.hpp"
#include "DataWriterVariable.hpp"
#include "DistributedVectorFactory.hpp"

/**
 * A concrete HDF5 data writer class.
 */
class Hdf5DataWriter : public AbstractHdf5Access //: public AbstractDataWriter
{
private:

    /** The factory to use in creating PETSc Vec and DistributedVector objects. */
    DistributedVectorFactory& mrVectorFactory;

    bool mCleanDirectory;                           /**< Whether to wipe the output directory */
    bool mUseExistingFile;                          /**< Whether we are using an existing file (for extending existing dataset, or adding a new one)*/
    bool mIsInDefineMode;                           /**< Is the DataWriter in define mode or not */
    bool mIsFixedDimensionSet;                      /**< Is the fixed dimension set */

    unsigned mEstimatedUnlimitedLength;             /**< An estimate of the unlimited dimension length for performance reasons. */
    unsigned mFileFixedDimensionSize;               /**< The size of the fixed dimension (number of rows)*/
    unsigned mDataFixedDimensionSize;               /**< The size of the fixed dimension (size of the vector of nodes)*/
    unsigned mLo;                                   /**< Local ownership of a PETSc vector of size #mFileFixedDimensionSize*/
    unsigned mHi;                                   /**< Local ownership of a PETSc vector of size #mFileFixedDimensionSize*/
    unsigned mNumberOwned;                          /**< mNumberOwned=#mHi-#mLo; except with incomplete data*/
    unsigned mOffset;                               /**< mOffset=#mLo; except with incomplete data*/

    bool mNeedExtend;                               /**< Used so that the data set is only extended when data is written*/
    bool mUseMatrixForIncompleteData;               /**< Whether to use a matrix format for incomplete data */

    std::vector<DataWriterVariable> mVariables;     /**< The data variables */

    long unsigned mCurrentTimeStep;                          /**< The current time step. */

    Mat mSinglePermutation;                         /**< Stores a permutation as a matrix */
    Mat mDoublePermutation;                         /**< Stores a permutation of a striped structure (u_0 v_0 u_1 v_1) as a matrix */

    Mat mSingleIncompleteOutputMatrix;              /**< Stores nodes to be output as a matrix */
    Mat mDoubleIncompleteOutputMatrix;              /**< Stores striped nodes to be output as a matrix */

    bool mUseOptimalChunkSizeAlgorithm;             /**< Whether to use the built-in algorithm for optimal chunk size */
    unsigned mFixedChunkSize;                       /**< User-provided chunk size */

    /**
     * Check name of variable is allowed, i.e. contains only alphanumeric & _, and isn't blank.
     *
     * @param rName variable name
     */
    void CheckVariableName(const std::string& rName);

    /**
     * Check name of unit is allowed, i.e. contains only alphanumeric & _, and isn't blank.
     *
     * @param rName unit name
     */
    void CheckUnitsName(const std::string& rName);

    /**
     * Compute #mOffset and #mNumberOwned from #mIncompleteNodeIndices when we have incomplete data.
     */
    void ComputeIncompleteOffset();

    /**
     * Opens an existing file or creates a new file, depending on #mUseExistingFile.
     *
     * This method sets #mFileId.
     */
    void OpenFile();

public:

    /**
     * Constructor.
     *
     * @param rVectorFactory the factory to use in creating PETSc Vec and DistributedVector objects.
     * @param rDirectory  the directory in which to write the data to file, relative to chaste test output.
     * @param rBaseName  the name of the file in which to write the data
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     * @param extendData  whether to try opening an existing file and appending to it.
     * @param datasetName The name of the HDF5 dataset to write, defaults to "Data".
     *
     * The extendData parameter allows us to add to an existing dataset.  It only really makes
     * sense if the existing file has an unlimited dimension which we can extend.  It also only
     * makes sense if cleanDirectory is false, otherwise there won't be a file there to read...
     */
    Hdf5DataWriter(DistributedVectorFactory& rVectorFactory,
                   const std::string& rDirectory,
                   const std::string& rBaseName,
                   bool cleanDirectory=true,
                   bool extendData=false,
                   std::string datasetName="Data");

    /**
     * Destructor.
     */
    virtual ~Hdf5DataWriter();

    /**
     * Define the fixed dimension, assuming complete data output (all the nodes).
     *
     * @param dimensionSize The size of the dimension
     */
    void DefineFixedDimension(long dimensionSize);

    /**
     * Define the fixed dimension, assuming incomplete data output (subset of the nodes).
     *
     * @param rNodesToOuput Node indexes to be output (precondition: to be monotonic increasing)
     * @param vecSize
     */
    void DefineFixedDimension(const std::vector<unsigned>& rNodesToOuput, long vecSize);

    /**
     * Define a variable.
     *
     * @param rVariableName The name of the dimension
     * @param rVariableUnits The physical units of the dimension
     * @param estimatedLength An estimate of the unlimited dimension length for performance reasons
     *
     */
    void DefineUnlimitedDimension(const std::string& rVariableName, const std::string& rVariableUnits, unsigned estimatedLength = 1);

    /**
     * Advance along the unlimited dimension. Normally this will be called
     * when all variables in a row have been input.
     */
    void AdvanceAlongUnlimitedDimension();

    /**
     * Define a variable.
     *
     * @param rVariableName The name of the dimension
     * @param rVariableUnits The physical units of the dimension
     *
     * @return The identifier of the variable
     */
    int DefineVariable(const std::string& rVariableName, const std::string& rVariableUnits);

    /**
     * End the define mode of the DataWriter.
     */
    virtual void EndDefineMode();

    /**
     * Extend the dataset to the correct to the correct dimensions if needed.
     */
    void PossiblyExtend();

    /**
     * Write data for a given variable from a PETSc vector to the dataset.
     *
     * @param variableID the variable id of the column in the HDF5 file
     * @param petscVector the data
     */
    void PutVector(int variableID, Vec petscVector);

    /**
     * Write data for two variables from a PETSc vector to the dataset.
     *
     * @param variableIDs vector of IDs of variables of the columns in the HDF5 file
     * @param petscVector the data
     */
    void PutStripedVector(std::vector<int> variableIDs, Vec petscVector);

    /**
     * Write a single value for the unlimited variable (e.g. time) to the dataset.
     *
     * @param value the data
     */
    void PutUnlimitedVariable(double value);

    /**
     * Close any open files.
     */
    void Close();

    /**
     * Get the id of the given variable, the variable must already
     * exist or an exception will be thrown.
     *
     * @param rVariableName  variable name to look up
     * @return  HDF5 id for the given variable.
     */
    int GetVariableByName(const std::string& rVariableName);


    /**
     * Apply a permutation to all occurences of PutVector
     * Should be called when in define mode
     * @param rPermutation  a forward/?reverse permutation
     * @return success value.  A value "false" indictates that the permutation was empty or was the identity and was not applied
     */
    bool ApplyPermutation(const std::vector<unsigned>& rPermutation);

     /**
     * Define the fixed dimension, assuming incomplete data output (subset of the nodes) and using a matrix
     * to convert from full to incomplete output (rather than picking required data values out one at a time).
     *
     * @param rNodesToOuput Node indexes to be output (precondition: to be monotonic increasing)
     * @param vecSize
     */
    void DefineFixedDimensionUsingMatrix(const std::vector<unsigned>& rNodesToOuput, long vecSize);

    /**
     * Use a particular chunk size, ignoring the algorithm that figures out the optimal value. Useful for reducing
     * size of checkpoints. USE WITH CAUTION.
     *
     * @param chunkSize user provided chunk size
     */
    void SetFixedChunkSize(unsigned chunkSize);


};

#endif /*HDF5DATAWRITER_HPP_*/
