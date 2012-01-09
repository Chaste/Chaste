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
