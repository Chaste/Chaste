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
#ifndef COMPAREHDF5RESULTSFILES_HPP_
#define COMPAREHDF5RESULTSFILES_HPP_

#include <string>

/**
 * Method to compare datasets (and associated unlimited variables) from two HDF5 files.
 *
 * @param pathname1  Path to first file
 * @param filename1  Filename of first file
 * @param makeAbsolute1  Whether the h5 file should be treated as relative to Chaste test output, otherwise treated as CWD or absolute.
 * @param pathname2  Path to second file
 * @param filename2  Filename of second file
 * @param makeAbsolute2  Whether the h5 file should be treated as relative to Chaste test output, otherwise treated as CWD or absolute.
 * @param tol  Numerical tolerance to use when comparing data values
 * @param datasetName  the dataset to compare (defaults to "Data")
 * @return True if the file contents are the same, False if they differ.
 */
bool CompareFilesViaHdf5DataReader(std::string pathname1, std::string filename1, bool makeAbsolute1,
                                   std::string pathname2, std::string filename2, bool makeAbsolute2,
                                   double tol=1e-10, std::string datasetName = "Data");

/**
 * Alternative, weaker, method for comparing two files. It computes the global norms of the two data
 * vectors and checks whether the difference is less than 1e-10
 *
 * @param pathname1  Path to first file
 * @param filename1  Filename of first file
 * @param makeAbsolute1  Whether the h5 file should be treated as relative to Chaste test output, otherwise treated as CWD or absolute.
 * @param pathname2  Path to second file
 * @param filename2  Filename of second file
 * @param makeAbsolute2  Whether the h5 file should be treated as relative to Chaste test output, otherwise treated as CWD or absolute.
 * @param tol  Numerical tolerance to use when comparing data values
 * @param datasetName  the dataset to compare (defaults to "Data")
 * @return true if the global norms differ by less than 1e-10 for each time step
 * */
bool CompareFilesViaHdf5DataReaderGlobalNorm(std::string pathname1, std::string filename1, bool makeAbsolute1,
                                             std::string pathname2, std::string filename2, bool makeAbsolute2,
                                             double tol=1e-10, std::string datasetName = "Data");

#endif /*COMPAREHDF5RESULTSFILES_HPP_*/
