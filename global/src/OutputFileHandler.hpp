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

#ifndef OUTPUTFILEHANDLER_HPP_
#define OUTPUTFILEHANDLER_HPP_

#include <fstream>
#include <ios>
#include <memory>
#include <string>

class FileFinder; // Avoid circular includes

/** Type of our output streams; a managed pointer to an std::ofstream. */
typedef std::shared_ptr<std::ofstream> out_stream;

/**
 * This file abstracts stuff that needs to be done when creating output files for tests.
 * It defines some helpful functions, so that things that are otherwise repeated in lots
 * of places are just done here.
 */
class OutputFileHandler
{
public:
    /**
     * Create an OutputFileHandler that will create output files in the given directory.
     * The directory name should be relative to the place where Chaste test output is
     * stored.  If the user needs to know where this is, use the GetChasteTestOutputDirectory
     * method.
     *
     * Will check that the directory exists and create it if needed.
     *
     * @note This MUST be called collectively, since the directory creation routine contains a barrier call.
     *
     * @param rDirectory  the directory to put output files in.
     * @param cleanOutputDirectory  whether to remove any existing files in the output directory
     */
    OutputFileHandler(const std::string& rDirectory,
                      bool cleanOutputDirectory = true);

    /**
     * Create an OutputFileHandler that will create output files in the given directory.
     * This must be located inside the folder where Chaste test output is stored, and will typically
     * be obtained with the FindFile method of this class from a parent handler.
     *
     * Will check that the directory exists and create it if needed.
     *
     * @note This MUST be called collectively, since it contains a barrier call.
     *
     * @param rDirectory  the directory to put output files in
     * @param cleanOutputDirectory  whether to remove any existing files in the output directory
     */
    OutputFileHandler(const FileFinder& rDirectory,
                      bool cleanOutputDirectory = true);

    /**
     * Static method for getting the test output directory (the directory where
     * chaste stores test output).  This is set from the environment variable
     * CHASTE_TEST_OUTPUT, and defaults to "/tmp/$USER/testoutput" if it is not set.
     *
     * Attempts to return an absolute path, but may get confused by odd setups.
     *
     * Static so an output file handler does not have to be created if the test output
     * directory is wanted for, say, reading a file.
     * @return the test output directory name (CHASTE_TEST_OUTPUT)
     */
    static std::string GetChasteTestOutputDirectory();

    /**
     * @return the full pathname to the directory this object will create files
     * in.
     */
    std::string GetOutputDirectoryFullPath() const;

    /**
     * @return the path to this output directory, relative to the Chaste root output folder.
     */
    std::string GetRelativePath() const;

    /**
     * Helper method to set up ArchiveLocationInfo.
     */
    void SetArchiveDirectory() const;

    /**
     * Open an output file in our directory, and check it was opened successfully.
     * Throws an Exception if not.
     *
     * @param rFileName  the name of the file to open, relative to the output directory.
     * @param mode  optionally, flags to use when opening the file (defaults are as for
     *         std::ofstream).
     * @return  a managed pointer to the opened file stream.
     */
    out_stream OpenOutputFile(const std::string& rFileName,
                              std::ios_base::openmode mode = std::ios::out | std::ios::trunc) const;

    /**
     * This just calls the other OpenOutputFile after concatenating the first three arguments
     * together to make the full filename. For example OpenOutputFile("results_", 3, ".dat")
     * creates results_3.dat. See documentation for
     * OpenOutputFile(std::string, std::ios_base::openmode).
     *
     * @param rFileName  the root name of the file to open
     * @param number  the number to append to the root name of the file
     * @param rFileFormat  the file format (extension)
     * @param mode  optionally, flags to use when opening the file (defaults are as for
     *         std::ofstream).
     * @return  a managed pointer to the opened file stream.
     */
    out_stream OpenOutputFile(const std::string& rFileName,
                              unsigned number,
                              const std::string& rFileFormat,
                              std::ios_base::openmode mode = std::ios::out | std::ios::trunc) const;

    /**
     * Copy the given file to this output directory.
     *
     * @note This MUST be called collectively, since it contains a barrier call.
     *
     * @param rSourceFile  the file to copy
     * @return the copied file
     */
    FileFinder CopyFileTo(const FileFinder& rSourceFile) const;

    /**
     * @return a FileFinder for a file in this output directory.
     *
     * @param leafName  the name of the file to find
     */
    FileFinder FindFile(std::string leafName) const;

    /** The name of the Chaste signature file added to folders we create. */
    static const std::string SIG_FILE_NAME;

private:
    std::string mDirectory; ///< The directory to store output files in (always ends in "/")

    /**
     * Functionality common to both constructors.
     *
     * @param rDirectory  relative path to the directory to put output files in
     * @param cleanOutputDirectory  whether to remove any existing files in the output directory
     */
    void CommonConstructor(const std::string& rDirectory,
                           bool cleanOutputDirectory);

    /**
     * Takes a reference to a string and adds a trailing slash if one is not already present
     *
     * @param rDirectory  The directory name to add a trailing slash to.
     */
    static void AddTrailingSlash(std::string& rDirectory);

    /**
     * Check that the desired output directory exists and is writable by us.
     * Create it if needed.
     * Return the full pathname of the output directory.
     *
     * The environment variable CHASTE_TEST_OUTPUT will be examined.  If it is set
     * and non-empty it is taken to be a directory where test output should be stored.
     * Otherwise the current directory is used.
     *
     * @note Contains a barrier, so must be called collectively.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste
     *         output will be stored (user shouldn't care about this).
     * @return full pathname to the output directory
     */
    std::string MakeFoldersAndReturnFullPath(const std::string& rDirectory) const;
};

#endif /*OUTPUTFILEHANDLER_HPP_*/
