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

#include "ArchiveLocationInfo.hpp"

#include <sstream>

#include "Exception.hpp"
#include "OutputFileHandler.hpp"

std::string ArchiveLocationInfo::mDirAbsPath = "";
std::string ArchiveLocationInfo::mMeshFilename = "mesh";

void ArchiveLocationInfo::SetMeshPathname(const FileFinder& rDirectory, const std::string& rFilename)
{
    SetArchiveDirectory(rDirectory);
    mMeshFilename = rFilename;
}

void ArchiveLocationInfo::SetMeshPathname(const std::string& rDirectory, const std::string& rFilename)
{
    bool is_absolute = FileFinder::IsAbsolutePath(rDirectory);
    RelativeTo::Value relative_to;
    if (!is_absolute)
    {
        relative_to = RelativeTo::ChasteTestOutput;
    }
    else
    {
        relative_to = RelativeTo::Absolute;
    }
    FileFinder dir(rDirectory, relative_to);
    SetArchiveDirectory(dir);
    mMeshFilename = rFilename;
}

void ArchiveLocationInfo::SetMeshFilename(const std::string& rFilename)
{
    mMeshFilename = rFilename;
}

std::string ArchiveLocationInfo::GetMeshFilename()
{
    if (mMeshFilename == "")
    {
        EXCEPTION("ArchiveLocationInfo::mMeshFilename has not been set");
    }
    return mMeshFilename;
}

std::string ArchiveLocationInfo::GetArchiveDirectory()
{
    if (mDirAbsPath == "")
    {
        EXCEPTION("ArchiveLocationInfo::mDirAbsPath has not been set");
    }
    return mDirAbsPath;
}

std::string ArchiveLocationInfo::GetProcessUniqueFilePath(
        const std::string& rFileName,
        unsigned procId)
{
    std::stringstream filepath;
    filepath << GetArchiveDirectory() << rFileName << "." << procId;
    return filepath.str();
}

void ArchiveLocationInfo::SetArchiveDirectory(const FileFinder& rDirectory)
{
    mDirAbsPath = rDirectory.GetAbsolutePath();
    if (!(*(mDirAbsPath.end()-1) == '/'))
    {
        mDirAbsPath = mDirAbsPath + "/";
    }
}

std::string ArchiveLocationInfo::GetArchiveRelativePath()
{
    std::string chaste_output = OutputFileHandler::GetChasteTestOutputDirectory();

    // Find occurrence of CHASTE_TEST_OUTPUT in string
    std::string::size_type pos = mDirAbsPath.find(chaste_output, 0);
    if (pos == 0)
    {
        return mDirAbsPath.substr(chaste_output.length());
    }
    else
    {
        return mDirAbsPath;
    }
}

bool ArchiveLocationInfo::GetIsDirRelativeToChasteTestOutput()
{
    std::string chaste_output = OutputFileHandler::GetChasteTestOutputDirectory();
    std::string::size_type pos = mDirAbsPath.find(chaste_output, 0);
    return (pos == 0);
}
