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

#ifndef DISTRIBUTEDTETRAHEDRALMESHPARTITIONTYPE_HPP_
#define DISTRIBUTEDTETRAHEDRALMESHPARTITIONTYPE_HPP_

/** Definition of partition types.
 * "DUMB" is using natural mesh ordering with PETSC_DECIDE.
 * "PARMETIS_LIBRARY" is a call to the parallel parMETIS library
 * "METIS_LIBRARY" is a call to the sequential METIS library
 * "PETSC_MAT_PARTITION" is a call to parMETIS (or whatever) via PETSc functionality
 */
struct DistributedTetrahedralMeshPartitionType
{
    /** The actual type enumeration */
    typedef enum
    {
        DUMB=0,
        PARMETIS_LIBRARY=1,
        METIS_LIBRARY=2,
        PETSC_MAT_PARTITION=3
    } type;
};

#endif /*DISTRIBUTEDTETRAHEDRALMESHPARTITIONTYPE_HPP_*/
