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

#ifndef DISTRIBUTEDTETRAHEDRALMESHPARTITIONTYPE_HPP_
#define DISTRIBUTEDTETRAHEDRALMESHPARTITIONTYPE_HPP_

/** Definition of partition types.
 * "DUMB" is using natural mesh ordering with PETSC_DECIDE.
 * "PARMETIS_LIBRARY" is a call to the parallel parMETIS library
 * "METIS_LIBRARY" used to be a call to the sequential METIS library.  (Now deprecated in favour of a drop through call to parMETIS.)
 * "PETSC_MAT_PARTITION" is a call to parMETIS (or whatever) via PETSc functionality.  This is not always available on a given installation.
 * "GEOMETRIC" requires user to define which region of space is owned by each process.
 */
struct DistributedTetrahedralMeshPartitionType
{
    /** The actual type enumeration */
    typedef enum
    {
        DUMB=0,
        PARMETIS_LIBRARY=1,  // Deprecated
        METIS_LIBRARY=2,
        PETSC_MAT_PARTITION=3,
        GEOMETRIC=4
    } type;
};

#endif /*DISTRIBUTEDTETRAHEDRALMESHPARTITIONTYPE_HPP_*/
