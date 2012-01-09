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
     * Get whether the directory provided is relative to CHASTE_TEST_OUTPUT.
     */
    static bool GetIsDirRelativeToChasteTestOutput();
};

#endif /*ARCHIVELOCATIONINFO_HPP_*/
