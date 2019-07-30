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

#ifndef CITATIONS_HPP_
#define CITATIONS_HPP_

#include <fstream>
#include <iostream>
#include <vector>
#include "CommandLineArguments.hpp"
#include "PetscTools.hpp"

/**
 * A class to register citations placed throughout the codebase, and pass them
 * through to PETSc (3.5+) or save them and print to screen/disk at program exit
 * (pre-PETSc 3.5). The behaviour from the user's point of view should be the
 * same regardless of PETSc version.
 *
 * This feature may be switched on by passing the "-citations" flag at runtime
 * to print to screen, or to file.txt using "-citations file.txt". Paths for the
 * latter can be relative or absolute, if the file can't be written (because the
 * directory doesn't exist, or otherwise) you'll get an exception.
 *
 * See r22760 for an example of adding a new citation to the trunk and registering
 * it.
 */
class Citations
{
public:
    /**
     * Method to register a citation. Just passes through to PETSc for v3.5+.
     *
     * @param pCitation  the citation in BibTeX format.
     * @param pSet  pointer to a flag used to track whether the citation
     *     has been added before; will be set to true by this method.
     */
    static void Register(const char* pCitation, PetscBool* pSet);

    /**
     * Print the list of citations, called automatically when finalising the PETSc
     * environment (see PetscSetupUtils).
     */
    static void Print();

private:
    // We need to manually switch on mCitations for testing.
    friend class TestCitations;

    /** The list of citations if using Chaste's built-in manager (pre-PETSc 3.5 or PETSc not initialised). */
    static std::vector<const char*> mCitations;

    /** Whether to use Chaste's own implementation rather than the PETSc one (e.g. if PETSc hasn't been initialised, or is too old). */
    static bool mUseChasteImplementation;
};

#endif // CITATIONS_HPP_
