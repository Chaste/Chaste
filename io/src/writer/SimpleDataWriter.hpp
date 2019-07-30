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

#ifndef SIMPLEDATAWRITER_HPP_
#define SIMPLEDATAWRITER_HPP_

#include <string>
#include <vector>

/**
 * A basic data writer that is easier to use than ColumnDataWriter but has less
 * functionality. NOTE: this is not an efficient writer.
 *
 * This class does not writer header lines so is ideal for immediately reading
 * with MATLAB or Gnuplot.
 */
class SimpleDataWriter
{
public:

    /**
     * Write the provided data out to the given file in columns
     *
     * @param rDirectory The directory, relative to TEST_OUTPUT
     * @param rFileName  The full file name (no format will be apended)
     * @param rData     The data. data[0] will written as the first column, data[1] the
     *                  second, and so on. An exception is thrown if they are not the same size
     * @param cleanDirectory Whether to clean the directory (defaults to true)
     */
    SimpleDataWriter(const std::string& rDirectory,
                     const std::string& rFileName,
                     const std::vector<std::vector<double> >& rData,
                     bool cleanDirectory=true);

    /**
     * Write the provided data out to the given file in 2 columns
     *
     * @param rDirectory The directory, relative to TEST_OUTPUT
     * @param rFileName  The full file name (no format will be apended)
     * @param rT        The first column of data
     * @param rX        The second column of data. An exception is thrown if the size
     *             of x is not the same as the size of t.
     * @param cleanDirectory Whether to clean the directory (defaults to true)
     */
    SimpleDataWriter(const std::string& rDirectory,
                     const std::string& rFileName,
                     const std::vector<double>& rT,
                     const std::vector<double>& rX,
                     bool cleanDirectory=true);

    /**
     * Write the provided data out to the given file in one column
     *
     * @param rDirectory The directory, relative to TEST_OUTPUT
     * @param rFileName  The full file name (no format will be apended)
     * @param rData      A std::vec of data
     * @param cleanDirectory Whether to clean the directory (defaults to true)
     */
    SimpleDataWriter(const std::string& rDirectory,
                     const std::string& rFileName,
                     const std::vector<double>& rData,
                     bool cleanDirectory=true);
};

#endif /*SIMPLEDATAWRITER_HPP_*/
