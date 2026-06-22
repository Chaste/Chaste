/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef FILESYSTEMPERMISSIONS_HPP_
#define FILESYSTEMPERMISSIONS_HPP_

#include "FileFinder.hpp"

/**
 * Helper methods for creating and copying filesystem entries with Chaste's
 * standard output permissions.
 */
class FilesystemPermissions
{
public:
    /**
     * @return standard Chaste permissions for created files.
     */
    static fs::perms GetFilePermissions()
    {
        return fs::perms::owner_read | fs::perms::owner_write | fs::perms::group_read;
    }

    /**
     * @return standard Chaste permissions for created directories.
     */
    static fs::perms GetDirectoryPermissions()
    {
        return fs::perms::owner_all | fs::perms::group_read | fs::perms::group_exec
               | fs::perms::others_read | fs::perms::others_exec;
    }

    /**
     * Set standard Chaste file permissions.
     *
     * @param rPath  path to a file
     */
    static void SetFilePermissions(const fs::path& rPath)
    {
        fs::permissions(rPath, GetFilePermissions(), fs::perm_options::replace);
    }

    /**
     * Set standard Chaste directory permissions.
     *
     * @param rPath  path to a directory
     */
    static void SetDirectoryPermissions(const fs::path& rPath)
    {
        fs::permissions(rPath, GetDirectoryPermissions(), fs::perm_options::replace);
    }

    /**
     * Create a directory and set standard Chaste directory permissions if it
     * was created.
     *
     * @param rPath  path to create
     * @return true if the directory was created
     */
    static bool CreateDirectoryWithPermissions(const fs::path& rPath)
    {
        bool created_dir = fs::create_directory(rPath);
        if (created_dir)
        {
            SetDirectoryPermissions(rPath);
        }
        return created_dir;
    }

    /**
     * Create directories and set standard Chaste directory permissions on each
     * directory that is created.
     *
     * @param rPath  path to create
     * @return true if any directory was created
     */
    static bool CreateDirectoriesWithPermissions(const fs::path& rPath)
    {
        if (fs::exists(rPath))
        {
            return false;
        }

        fs::path parent_path = rPath.parent_path();
        bool created_parent = false;
        if (!parent_path.empty() && parent_path != rPath && !fs::exists(parent_path))
        {
            created_parent = CreateDirectoriesWithPermissions(parent_path);
        }

        return CreateDirectoryWithPermissions(rPath) || created_parent;
    }

    /**
     * Copy a file and set standard Chaste file permissions on the destination.
     *
     * @param rFromPath  source file
     * @param rToPath  destination file
     */
    static void CopyFileWithPermissions(const fs::path& rFromPath, const fs::path& rToPath)
    {
        fs::copy_file(rFromPath, rToPath);
        SetFilePermissions(rToPath);
    }
};

#endif /*FILESYSTEMPERMISSIONS_HPP_*/
