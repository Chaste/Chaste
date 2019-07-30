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

#include "Citations.hpp"
#include "PetscTools.hpp"

// Initialise 'member' variables
std::vector<const char*> Citations::mCitations;
bool Citations::mUseChasteImplementation = false;

void Citations::Register(const char pCitation[], PetscBool* pSet)
{
    // Figure out if we can use PETSc's implementation.
    // Note that the effect of this if() is to ensure we don't switch implementation if PETSc is initialised after the first citation is registered.
    // (If PETSc is finalised after the first citation is registered we will probably switch implementation, but there's no way to avoid that.)
    if (!mUseChasteImplementation)
    {
#if ((PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5) || PETSC_VERSION_MAJOR > 3)
        mUseChasteImplementation = !PetscTools::IsInitialised();
#else
        mUseChasteImplementation = true;
#endif
    }

    if (!mUseChasteImplementation)
    {
#if ((PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5) || PETSC_VERSION_MAJOR > 3)
        PetscCitationsRegister(pCitation, pSet);
#endif
    }
    else if (!(*pSet))
    {
        // Not yet added this one
        mCitations.push_back(pCitation);
        *pSet = PETSC_TRUE;
    }
}

void Citations::Print()
{
    if (mUseChasteImplementation)
    {
        if (PetscTools::AmMaster() && CommandLineArguments::Instance()->OptionExists("-citations"))
        {
            std::ostream* p_output = &(std::cout);
            bool writing_to_file = false;

            if (CommandLineArguments::Instance()->GetNumberOfArgumentsForOption("-citations") > 0)
            {
                // We've got a file to write to - assume the user has given us something sensible!
                writing_to_file = true;
                std::string out_file_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-citations");
                p_output = new std::ofstream(out_file_path.c_str(), std::ios::out);
                EXCEPT_IF_NOT(p_output->good());
            }

            /* Write header */
            (*p_output) << "If you publish results based on this computation please cite the following:" << std::endl;
            (*p_output) << "===========================================================================" << std::endl;
            /* Write citations */
            for (unsigned i = 0; i < mCitations.size(); ++i)
            {
                (*p_output) << mCitations[i];
            }
            /* Write footer */
            (*p_output) << "===========================================================================" << std::endl;

            if (writing_to_file)
            {
                delete p_output;
            }
        }
    }
}
