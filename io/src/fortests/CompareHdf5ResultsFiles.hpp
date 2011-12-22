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
#ifndef COMPAREHDF5RESULTSFILES_HPP_
#define COMPAREHDF5RESULTSFILES_HPP_

#include <string>

bool CompareFilesViaHdf5DataReader(std::string pathname1, std::string filename1, bool makeAbsolute1,
                                   std::string pathname2, std::string filename2, bool makeAbsolute2,
                                   double tol=1e-10);

/**
 * Alternative, weaker, method for comparing two files. It computes the global norms of the two data
 * vectors and checks whether the difference is less than 1e-10
 *
 * @return true if the global norms differ by less than 1e-10 for each time step
 * */
bool CompareFilesViaHdf5DataReaderGlobalNorm(std::string pathname1, std::string filename1, bool makeAbsolute1,
                                             std::string pathname2, std::string filename2, bool makeAbsolute2,
                                             double tol=1e-10);

#endif /*COMPAREHDF5RESULTSFILES_HPP_*/
