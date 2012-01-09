/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef EXECUTABLESUPPORT_HPP_
#define EXECUTABLESUPPORT_HPP_

#include <string>
#include "OutputFileHandler.hpp"

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
     * Write information about library and compiler versions to the provenance_info.txt
     * output file.
     *
     * @param outFile the provenance_info.txt file.
     */
    static void WriteLibraryInfo( out_stream &outFile );

    /**
     * Call InitializePetsc, ShowCopyright, then ShowParallelLaunching.
     *
     * @param pArgc  pointer to the number of arguments
     * @param pArgv  pointer to the argument list
     */
    static void StandardStartup(int* pArgc, char*** pArgv);

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
    static std::string mOutputDirectory;
};

#endif /* EXECUTABLESUPPORT_HPP_ */
