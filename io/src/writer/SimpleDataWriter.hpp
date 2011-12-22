/*

Copyright (C) University of Oxford, 2005-2011

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
