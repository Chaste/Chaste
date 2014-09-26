/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef CITATIONS_HPP_
#define CITATIONS_HPP_

#include <vector>
#include <iostream>
#include <fstream>
#include "PetscTools.hpp"
#include "CommandLineArguments.hpp"

/**
 * A class to register citations placed throughout the codebase, and pass them
 * through to PETSc (3.5+) or save them and print to screen/disk at program exit
 * (pre-PETSc 3.5).
 *
 * This feature may be switched on by passing the "-citations" flag at runtime
 * to print to screen, or to file.txt using "-citations file.txt". Paths for the
 * latter can be relative or absolute, but have NO error checking, so care is
 * advised.
 * 
 * See r22760 for an example of adding a new citation to the trunk and registering
 * it.
 */ 
class Citations
{
public:
    /** 
     * Method to register a citation. Just passes through the PETSc for v 3.5+.
     * 
     * @param cit  The citation in bibtex format.
     * @param set  A PetscBool that starts as false, used to track whether the citation
     * has been added before.
     */
    static void Register(const char * cit, PetscBool * set);

    /**
     * Print the list of citations, called automatically when finalising the PETSc
     * environment.
     */
    static void Print();

private:
#if ( PETSC_VERSION_MAJOR<3 || PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<5 )
    /** The list of citations if using Chaste's built-in manager (pre-PETSc 3.5). */
    static std::vector<const char *> mCitations;
#endif

};

#endif // CITATIONS_HPP_
