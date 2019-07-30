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

#ifndef COMMANDLINEARGUMENTSMOCKER_HPP_
#define COMMANDLINEARGUMENTSMOCKER_HPP_


/**
 * This is a little helper wrapper for CommandLineArguments that
 * allows a test to replace ALL of the command line arguments
 * by calling the constructor.
 *
 * It is useful for allowing more tests on code that reads
 * CommandLineArguments without resorting to acceptance tests of executables.
 */
class CommandLineArgumentsMocker
{
public:
    /**
     * Only available constructor.
     * This replaces all of the command line arguments that were being used (but remembers them).
     *
     * @param newArguments  a string of command line arguments to put into the singleton CommandLineArguments class.
     */
    CommandLineArgumentsMocker(std::string newArguments);

    /**
     * Destructor.
     *
     * Cleans up memory and restores the original CommandLineArguments.
     */
    ~CommandLineArgumentsMocker();

private:

    /** The new number of arguments */
    int mNumArgs;
    /** The new arguments */
    char** mpArgs;

    /** The location of the old number of arguments */
    int* mpNumOldArgs;
    /** The location of the old arguments */
    char*** mpOldArgs;
};

#endif // COMMANDLINEARGUMENTSMOCKER_HPP_
