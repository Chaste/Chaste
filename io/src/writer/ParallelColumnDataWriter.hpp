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

#ifndef PARALLELCOLUMNDATAWRITER_HPP_
#define PARALLELCOLUMNDATAWRITER_HPP_

#include <string>

#include "ColumnDataWriter.hpp"
#include "DistributedVector.hpp"

#include <petscvec.h>

/**
 * A parallelised column data writer class.
 */
class ParallelColumnDataWriter  : public ColumnDataWriter
{
private:

    bool mIsParallel;     /**< Set to true in constructor if running in parallel*/
    Vec mConcentrated;    /**< Vector to hold concentrated copy of distributed vector on the master process*/
    VecScatter mToMaster; /**< variable holding information for concentrating a vector*/

public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the data to file
     * @param rBaseName  the name of the file in which to write the data
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     */
    ParallelColumnDataWriter(const std::string& rDirectory,
                             const std::string& rBaseName,
                             bool cleanDirectory=true);

    /**
     * Destructor.
     */
    virtual ~ParallelColumnDataWriter();

    /**
     * Write data for a given variable from a PETSc vector to the dataset.
     *
     * @param variableID the variable
     * @param petscVector the data
     */
    void PutVector(int variableID, Vec petscVector);

    /**
     * Write data for a given variable from a stripe to the dataset.
     *
     * @param variableId the variable
     * @param rStripe the data
     * \todo allow this to be a const-reference
     */
    void PutVectorStripe(int variableId, DistributedVector::Stripe& rStripe);

    /**
     * Input the variable value to the output file or ancillary file
     *
     * @param variableID
     * @param variableValue
     * @param dimensionPosition  The position in column (defaults to -1). This is required if
     *      there is a fixed dimension, and will be the position along that dimension
     */
    void PutVariable(int variableID, double variableValue, long dimensionPosition = -1);

    /**
     * End the define mode of the DataWriter.
     */
    void EndDefineMode();

    /**
     * Advance along the unlimited dimension. Normally this will be called
     * when all variables in a row have been input.
     */
    void AdvanceAlongUnlimitedDimension();

    /**
     * Close any open files.
     */
    void Close();
};

#endif /*PARALLELCOLUMNDATAWRITER_HPP_*/
