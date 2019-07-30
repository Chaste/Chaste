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
#ifndef ARCHIVELOCATIONINFO_HPP_
#define ARCHIVELOCATIONINFO_HPP_

#include <string>

#include "FileFinder.hpp"
#include "PetscTools.hpp"

/**
 * Mini-class to help with 'archiving' various classes that don't write their
 * data directly to the archive file.  They thus need to know information
 * about where the archive is being written to, in order to write their own
 * files into the same folder.  The main methods are GetArchiveDirectory and
 * SetArchiveDirectory.
 *
 * This functionality is used by the meshes, LinearSystem and HeartConfig.
 *
 * For the benefit of the meshes (and the cell_based code), there are also
 * shortcut methods SetMeshPathname and GetMeshFilename, allowing you to
 * specify the base file name for the mesh.  This is needed because the cell_based
 * code adds timestamp information to the file name.
 */
class ArchiveLocationInfo
{
private:

    /** Absolute path for directory that archives are being written to. */
    static std::string mDirAbsPath;

    /** Mesh filename (relative to #mDirAbsPath). */
    static std::string mMeshFilename;

public:

    /**
     * Set the location to write mesh files.
     *
     * @param rDirectory  the directory to write to.
     * @param rFilename  the base name (minus extension) for the mesh files.
     */
    static void SetMeshPathname(const FileFinder& rDirectory, const std::string& rFilename);

    /**
     * Set the location to write mesh files.
     *
     * @param rDirectory  the directory to write to (if relative, assumes relative to CHASTE_TEST_OUTPUT).
     * @param rFilename  the base name (minus extension) for the mesh files.
     */
    static void SetMeshPathname(const std::string& rDirectory, const std::string& rFilename);

    /**
     * Set the filename for mesh files.
     *
     * @param rFilename  the base name (minus extension) for the mesh files, used to put on a timestamp.
     */
    static void SetMeshFilename(const std::string& rFilename);

    /**
     * Get the filename that the mesh should be written to
     *
     * @return mesh filename
     */
    static std::string GetMeshFilename();

    /**
     * Get the directory that archives are being written to.
     * Will always end in a '/'.
     *
     * @return full path to directory
     */
    static std::string GetArchiveDirectory();

    /**
     * Get the full path to an output file which has a name unique to the current
     * process.  Useful for ensuring that each process writes to / reads from a
     * separate file when running in parallel.
     *
     * The path will have the form "path_to_output_dir/rFileName.process_rank"
     *
     * @param rFileName  the base file name
     * @param procId  the process id number (defaults to current process)
     * @return a full path to the file for this process
     */
    static std::string GetProcessUniqueFilePath(const std::string& rFileName,
                                                unsigned procId=PetscTools::GetMyRank());

    /**
     * Set the directory that archives are being written to.
     *
     * @param rDirectory  the directory in question.
     */
    static void SetArchiveDirectory(const FileFinder& rDirectory);

    /**
     * Get the directory to which the archives are being written.
     * Remove CHASTE_TEST_OUTPUT prefix (assuming it exists); if it doesn't, returns the absolute path.
     * Will always end in a '/'.
     *
     * @return relative path to directory
     */
    static std::string GetArchiveRelativePath();

    /**
     * @return true if the directory provided is relative to CHASTE_TEST_OUTPUT.
     */
    static bool GetIsDirRelativeToChasteTestOutput();
};

#endif /*ARCHIVELOCATIONINFO_HPP_*/
