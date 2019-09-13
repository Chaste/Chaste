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

#ifndef HDF5DATAWRITER_HPP_
#define HDF5DATAWRITER_HPP_

#include <vector>

#include "AbstractHdf5Access.hpp"
#include "DataWriterVariable.hpp"
#include "DistributedVectorFactory.hpp"

/**
 * A concrete HDF5 data writer class.
 */
class Hdf5DataWriter : public AbstractHdf5Access //: public AbstractDataWriter
{
    friend class TestHdf5DataWriter;
private:

    /** The factory to use in creating PETSc Vec and DistributedVector objects. */
    DistributedVectorFactory& mrVectorFactory;

    const bool mCleanDirectory;                     /**< Whether to wipe the output directory */
    const bool mUseExistingFile;                    /**< Whether we are using an existing file (for extending existing dataset, or adding a new one)*/
    bool mIsInDefineMode;                           /**< Is the DataWriter in define mode or not */
    bool mIsFixedDimensionSet;                      /**< Is the fixed dimension set */

    unsigned mEstimatedUnlimitedLength;             /**< An estimate of the unlimited dimension length for performance reasons. */
    unsigned mFileFixedDimensionSize;               /**< The size of the fixed dimension (number of rows)*/
    unsigned mDataFixedDimensionSize;               /**< The size of the fixed dimension (size of the vector of nodes)*/
    unsigned mLo;                                   /**< Local ownership of a PETSc vector of size #mFileFixedDimensionSize*/
    unsigned mHi;                                   /**< Local ownership of a PETSc vector of size #mFileFixedDimensionSize*/
    unsigned mNumberOwned;                          /**< mNumberOwned=#mHi-#mLo; except with incomplete data*/
    unsigned mOffset;                               /**< mOffset=#mLo; except with incomplete data*/
    std::vector<unsigned> mIncompletePermIndices;   /**< The actual runtime indices of the incomplete data nodes (used when running in parallel and a permutation has been applied)*/
    bool mNeedExtend;                               /**< Used so that the data set is only extended when data is written*/

    std::vector<DataWriterVariable> mVariables;     /**< The data variables */

    long unsigned mCurrentTimeStep;                 /**< The current time step. */

    Mat mSinglePermutation;                         /**< Stores a permutation as a matrix */
    Mat mDoublePermutation;                         /**< Stores a permutation of a striped structure (u_0 v_0 u_1 v_1) as a matrix */

    bool mUseOptimalChunkSizeAlgorithm;             /**< Whether to use the built-in algorithm for optimal chunk size */
    hsize_t mChunkSize[DATASET_DIMS];               /**< Stores chunk dimensions */
    hsize_t mNumberOfChunks;                        /**< The total number of chunks in the dataset */
    hsize_t mFixedChunkSize[DATASET_DIMS];          /**< User-provided chunk size */
    hsize_t mChunkTargetSize;                       /**< User-provided target chunk size (for the algorithm) */

    hsize_t mAlignment;                             /**< User-provided alignment parameter */

    bool mUseCache;                                 /**< Whether to use a cache */
    long unsigned mCacheFirstTimeStep;              /**< Coordinate to keep track of cache writes */
    std::vector<double> mDataCache;                 /**< Cache results here before writing */

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

    /**
     * Little method to calculate the number of chunks resulting from given chunk dimensions.
     * @return The number of chunks resulting from chunk dimensions (stored in member variable).
     */
    hsize_t CalculateNumberOfChunks();

    /**
     * Calculate (and save in the member variables) chunk dimensions based on a target size.
     * Chunks are kept as "square" as possible while wasting as little space around the edges
     * of the dataset as possible.
     * @param[in] targetSize The target number of entries in each dimension in the chunk.
     * @param[out] pChunkSizeInBytes The size of the resulting chunk in bytes.
     * @param[out] pAllOneChunk Whether the dataset is spanned by a single chunk.
     */
    void CalculateChunkDims( unsigned targetSize, unsigned* pChunkSizeInBytes, bool* pAllOneChunk );

    /**
     * This method sets the chunk size by building up in each dimension until a threshold is reached,
     * unless user-specified values have been set using #SetFixedChunkSize.
     * By default, chunks of 128 K are used, which seems to be a good compromise. For large problems
     * performance will usually improve by increasing this value (to e.g. 1 M).
     *
     * Note: The public method for altering the algorithm's target chunk size is #SetTargetChunkSize.
     */
    void SetChunkSize();

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
     * @param useCache  Whether to cache writes so only whole chunks are written to disk.
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
                   std::string datasetName="Data",
                   bool useCache=false);

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
     * @param rNodesToOuputOriginalIndices Node indices to be output as in the original mesh (precondition: to be monotonic increasing)
     * @param vecSize
     * @param rNodesToOuputPermutedIndices Node indices to be output as they are in the runtime parallel partitioned mesh
     */
    void DefineFixedDimension(const std::vector<unsigned>& rNodesToOuputOriginalIndices, const std::vector<unsigned>& rNodesToOuputPermutedIndices, long vecSize);

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
     * Check whether writer is in define mode.
     *
     * When in define mode variables can be defined but no data may be written, and vice versa.
     * When extending, the writer will be in define mode if the dataset does not exist, and vice
     * versa.
     *
     * @return whether in define mode.
     */
    bool IsInDefineMode();

    /**
     * End the define mode of the DataWriter.
     */
    virtual void EndDefineMode();

    /**
     * Extend the dataset to the correct to the correct dimensions if needed.
     */
    void PossiblyExtend();

    /**
     * Reset mCurrentTimeStep to 0 and resize the dataset and unlimited dataset to size 1
     * in the first dimension. Future writes will therefore overwrite the original contents.
     *
     * * THIS METHOD (EFFECTIVELY) DELETES THE DATASET! * *
     */
    void EmptyDataset();

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
     * Whether we're caching writes
     * @return whether we're caching writes
     */
    bool GetUsingCache();

    /**
     * Write the cache to disk.
     */
    void WriteCache();

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
     * @param unsafeExtendingMode is true when we are extending a file which requires a permutation to be applied to it.
     *        In particular we are extending a cardiac simulation with "original node ordering"
     * @return success value.  A value "false" indictates that the permutation was empty or was the identity and was not applied
     */
    bool ApplyPermutation(const std::vector<unsigned>& rPermutation, bool unsafeExtendingMode=false);

    /**
     * Use a particular chunk size, ignoring the algorithm that figures out a sensible value.
     *
     * This method may be useful in very specific circumstances, e.g. if you want chunks to align
     * perfectly with stripes on a Lustre file system, or if you have a parallel problem where
     * every process can be assigned the same size and shape hyperslab of the HDF5 file (thus
     * getting rid of contention between processes when reading/writing).
     *
     * USE WITH CAUTION - as it can degrade performance.
     *
     * @param rTimestepsPerChunk  The number of unlimited dimension steps (usually time steps) per chunk.
     * @param rNodesPerChunk  The number of objects (usually nodes) per chunk in the second dimension.
     * @param rVariablesPerChunk  The number of objects (usually output variables) per chunk in the third dimension.
     */
    void SetFixedChunkSize(const unsigned& rTimestepsPerChunk,
                           const unsigned& rNodesPerChunk,
                           const unsigned& rVariablesPerChunk);

    /*
     * * NOTES ON CHUNK SIZE AND ALIGNMENT *
     *
     * The default target chunk size is 128 K, which seems to be a good compromise
     * for small problems (e.g. on a desktop PC). For larger problems, I/O
     * performance often improves with increased chunk size. A sweet spot seems to
     * be 1 M chunks.
     *
     * On a striped filesystem, for best performance set the chunk size and
     * alignment (using `H5Pset_alignment` above) to the file stripe size. With
     * `H5Pset_alignment`, every chunk starts at a multiple of the alignment value.
     *
     * To avoid wasting space, the chunk size should be an integer multiple of the
     * alignment value. Note that the algorithm below automatically goes back one
     * step after exceeding the chunk size, which minimises wasted space. To see
     * why, consider the examples below.
     *
     * (Example 1) Say our file system uses 1 M stripes. If we set
     *     target_size_in_bytes = 1024*1024;
     * below and uncomment
     *     H5Pset_alignment(fapl, 0, 1024*1024);
     * above, i.e. aim for (slightly under) 1 M chunks and align them to 1 M
     * boundaries, then the algorithm below will get as close as possible to 1 M
     * chunks but not exceed it, so each chunk will be padded slightly to sit on
     * the 1 M boundaries. Each chunk will therefore have its own stripe on the
     * file system, which should give us the best bandwidth and least contention.
     * Conclusion: this is optimal!
     *
     * Note: In general the algorithm can get very close to the target so the
     * waste isn't bad. Typical utilization is 99.99% (check with "h5ls -v ...").
     *
     * (Example 2) We set
     *     target_size_in_bytes = 128*1024;
     * and uncomment
     *     H5Pset_alignment(fapl, 0, 1024*1024);
     * i.e. 128 K chunks aligned to 1 M boundaries. This would pad every chunk to
     * 1 M boundaries, wasting 7/8 of the space in the file! A file which might be
     * 5 G with an efficient layout would be more like 40 G! Conclusion: setting
     * the chunk size to less than the alignment value is very bad!
     *
     * (Example 3) Say our file system uses 1 M stripes. We set
     *     target_size_in_bytes = 2*1024*1024;
     * and uncomment
     *     H5Pset_alignment(fapl, 0, 1024*1024);
     * i.e. 2 M chunks aligned to 1 M boundaries. This might not be optimal, but
     * it's OK, since the chunk size is (slightly under) twice the alignment, as in
     * Example 1 the amount of padding would be very small. Each read/write would
     * require accessing 2 stripes on the file system. Conclusion: a chunk size of
     * an integer multiple of the alignment value is fine (but not optimal).
     */

    /**
     * Adjust the target (max) chunk size in the chunking algorithm. Useful
     * when one knows roughly how big a chunk should be but doesn't care about
     * the exact dimensions or shape. Default is 128 K.
     *
     * Especially useful with SetAlignment for ensuring each chunk gets its own
     * stripe on striped file systems.
     *
     * This method only has an effect when creating a NEW DATASET. Must be
     * called in define mode.
     *
     * @param targetSize Max chunk size (bytes)
     */
    void SetTargetChunkSize(hsize_t targetSize);

    /**
     * Set the alignment parameter to pass through to H5Pset_alignment. Every
     * file object will be aligned on the disk to a multiple of this parameter.
     * See the H5P docs for more information.
     *
     * Especially useful with SetTargetChunkSize to ensure each chunk (each
     * chunk is an 'file object') is aligned to a disk stripe on a striped file
     * system with minimal wastage.
     *
     * This method only has an effect when creating a NEW HDF5 FILE. Must be
     * called in define mode with extendData set to False.
     *
     * @param alignment Alignment (bytes)
     */
    void SetAlignment(hsize_t alignment);
};

#endif /*HDF5DATAWRITER_HPP_*/
