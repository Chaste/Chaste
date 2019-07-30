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

#ifndef EXECUTABLESUPPORT_HPP_
#define EXECUTABLESUPPORT_HPP_

#include <string>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp" // For out_stream

/**
 * Various helpful static methods for people writing their own executables
 * within the Chaste framework.
 *
 * Most executables will just need to call StandardStartup as the first thing in their
 * main() function, and FinalizePetsc before quitting.  The other methods allow you to
 * fine-tune what output is presented to users.
 */
class ExecutableSupport
{
public:
    /**
     * Initialise PETSc from the command line arguments.
     *
     * @param pArgc  pointer to the number of arguments
     * @param pArgv  pointer to the argument list
     */
    static void InitializePetsc(int* pArgc, char*** pArgv);

    /**
     * Display Chaste's copyright information.
     */
    static void ShowCopyright();

    /**
     * Output extra diagnostics when Chaste is launched in parallel.
     */
    static void ShowParallelLaunching();

    /**
     * Set the directory to which the files created by WriteMachineInfoFile,
     * and WriteProvenanceInfoFile will be written.
     * By default they will write to the CHASTE_TEST_OUTPUT folder.
     *
     * @param rOutputDirectory  the directory to write to
     */
    static void SetOutputDirectory(const std::string& rOutputDirectory);

    /**
     * Write to log files information about the machine that ran the code.
     * Each process will output its own file.
     *
     * @param fileBaseName base name of the file to write to
     */
    static void WriteMachineInfoFile(std::string fileBaseName);

    /**
     * Write to log files provenance information about this executable.
     * Each process will output its own file, named provenance_info.n.txt.
     */
    static void WriteProvenanceInfoFile();

    /**
     * Get information about library and compiler versions in XML-like format.
     *
     * @param rInfo a string to populate with library info.
     */
    static void GetBuildInfo(std::string& rInfo);

    /**
     * Call InitializePetsc, ShowCopyright, then ShowParallelLaunching.
     *
     * @param pArgc  pointer to the number of arguments
     * @param pArgv  pointer to the argument list
     */
    static void StandardStartup(int* pArgc, char*** pArgv);

    /**
     * Call InitializePetsc, then ShowParallelLaunching. Omit ShowCopyright.
     *
     * @param pArgc  pointer to the number of arguments
     * @param pArgv  pointer to the argument list
     */
    static void StartupWithoutShowingCopyright(int* pArgc, char*** pArgv);

    /**
     * Display an error message to the user, on stderr.
     *
     * @param rMessage  the message to display
     * @param masterOnly  whether only the master process should display the error
     */
    static void PrintError(const std::string& rMessage, bool masterOnly=false);

    /**
     * Display an informative message to the user, on stdout.
     * Message is only displayed by master process
     * @param rMessage  the message to display
     */
    static void Print(const std::string& rMessage);

    /**
     * Shut down PETSc so we exit cleanly.
     */
    static void FinalizePetsc();

    /**
     * Standard exit codes for executables to return from main():
     * successful termination.
     */
    static const int EXIT_OK = 0;

    /**
     * Standard exit codes for executables to return from main():
     * exception thrown during execution.
     */
    static const int EXIT_ERROR = 1;

    /**
     * Standard exit codes for executables to return from main():
     * bad arguments passed on command line.
     */
    static const int EXIT_BAD_ARGUMENTS = 2;

private:
    /** The output directory to put machine provenance information into. */
    static FileFinder mOutputDirectory;
};

#endif /* EXECUTABLESUPPORT_HPP_ */
