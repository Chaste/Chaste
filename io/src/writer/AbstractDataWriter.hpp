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

#ifndef ABSTRACTDATAWRITER_HPP
#define ABSTRACTDATAWRITER_HPP

#include <string>

/**
 * An abstract base class for data output.
 * Loosely based on NetCDF api.
 */
class AbstractDataWriter
{
public:

    /**
     * Define the fixed dimension.
     * This method must be overridden in concrete classes.
     *
     * @param rDimensionName The name of the dimension
     * @param rDimensionUnits The physical units of the dimension
     * @param dimensionSize The size of the dimension
     *
     * @return The identifier of the variable
     */
    virtual int DefineFixedDimension(const std::string& rDimensionName,
                                     const std::string& rDimensionUnits,
                                     long dimensionSize)=0;

    /**
     * Define the unlimited dimension, i.e. the dimension that increases as the simulation progresses.
     * This method must be overridden in concrete classes.
     *
     * @param rDimensionName The name of the unlimited dimension
     * @param rDimensionUnits The physical units of the unlimited dimension
     *
     * @return The identifier of the variable
     */
    virtual int  DefineUnlimitedDimension(const std::string& rDimensionName,
                                          const std::string& rDimensionUnits)=0;

    /**
     * Define a variable.
     * This method must be overridden in concrete classes.
     *
     * @param rVariableName The name of the dimension
     * @param rVariableUnits The physical units of the dimension
     *
     * @return The identifier of the variable
     */
    virtual int DefineVariable(const std::string& rVariableName,
                               const std::string& rVariableUnits)=0;

    /**
     * End the define mode of the DataWriter.
     * This method must be overridden in concrete classes.
     */
    virtual void EndDefineMode()=0;

    /**
     * Dummy function for DoAdvanceAlongUnlimitedDimension.
     * This method must be overridden in concrete classes.
     */
    virtual void AdvanceAlongUnlimitedDimension()=0;

    /**
     * Input the variable value to the output file or ancillary file
     * This method must be overridden in concrete classes.
     *
     * @param variableID
     * @param variableValue
     * @param dimensionPosition  The position in column (defaults to -1). This is required if
     *      there is a fixed dimension, and will be the position along that dimension
     */
    virtual void PutVariable(int variableID, double variableValue, long dimensionPosition=-1)=0;

    /**
     * Close any open files.
     * This method must be overridden in concrete classes.
     */
    virtual void Close()=0;

    /**
     * Base class with virtual methods needs a virtual destructor.
     */
    virtual ~AbstractDataWriter()
    {}
};

#endif //ABSTRACTDATAWRITER_HPP
