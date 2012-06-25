/*

Copyright (c) 2005-2012, University of Oxford.
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

#include "CellMLToSharedLibraryConverter.hpp"

#include <sstream>
#include <unistd.h> // For getpid()
#include <sys/stat.h> // For mkdir()
#include <ctime>
#include <cstring> // For strerror()
#include <cerrno> // For errno

#include "Exception.hpp"
#include "ChasteBuildRoot.hpp"
#include "PetscTools.hpp"
#include "DynamicModelLoaderRegistry.hpp"
#include "GetCurrentWorkingDirectory.hpp"

CellMLToSharedLibraryConverter::CellMLToSharedLibraryConverter(bool preserveGeneratedSources,
                                                               std::string component)
    : mPreserveGeneratedSources(preserveGeneratedSources),
      mComponentName(component)
{
}

DynamicCellModelLoaderPtr CellMLToSharedLibraryConverter::Convert(const FileFinder& rFilePath,
                                                                  bool isCollective)
{
    DynamicCellModelLoaderPtr p_loader;
    std::string absolute_path = rFilePath.GetAbsolutePath();
    // Check the file exists
    if (!rFilePath.Exists())
    {
        EXCEPTION("Dynamically loadable cell model '" + absolute_path + "' does not exist.");
    }
    // Find out whether rFilePath is a .cellml or .so
    size_t dot_position = absolute_path.find_last_of(".");
    if (dot_position == std::string::npos)
    {
        EXCEPTION("File does not have an extension: " + absolute_path);
    }
    std::string extension = absolute_path.substr(dot_position+1);
    if (extension == "cellml")
    {
        // Split the path into folder and leaf
        size_t slash_position = absolute_path.find_last_of("/\\");
        assert(slash_position != std::string::npos);
        std::string folder = absolute_path.substr(0, slash_position+1); // Include trailing slash
        std::string leaf = absolute_path.substr(slash_position+1, dot_position-slash_position); // Include dot
        std::string so_path = folder + "lib" + leaf + "so";
        // Does the .so file already exist (and was it modified after the .cellml?)
        FileFinder so_file(so_path, RelativeTo::Absolute);
        if (!so_file.Exists() || rFilePath.IsNewerThan(so_file))
        {
            if (!isCollective)
            {
                EXCEPTION("Unable to convert .cellml to .so unless called collectively, due to possible race conditions.");
            }
            ConvertCellmlToSo(absolute_path, folder, leaf);
        }
        // Load the .so
        p_loader = DynamicModelLoaderRegistry::Instance()->GetLoader(so_file);
    }
    else if (extension == "so")
    {
        // Just load the .so
        p_loader = DynamicModelLoaderRegistry::Instance()->GetLoader(rFilePath);
    }
    else
    {
        EXCEPTION("Unsupported extension '." + extension + "' of file '" + absolute_path + "'; must be .so or .cellml");
    }

    return p_loader;
}

void CellMLToSharedLibraryConverter::ConvertCellmlToSo(const std::string& rCellmlFullPath,
                                                       const std::string& rCellmlFolder,
                                                       const std::string& rModelLeafName)
{
    std::string tmp_folder, build_folder;
    std::string old_cwd = GetCurrentWorkingDirectory();
    // Check that the Chaste source tree exists
    FileFinder chaste_root("", RelativeTo::ChasteSourceRoot);
    if (!chaste_root.IsDir())
    {
        EXCEPTION("No Chaste source tree found at '" << chaste_root.GetAbsolutePath()
                  << "' - you need the source to use CellML models directly in Chaste.");
    }
    FileFinder component_dir(mComponentName, RelativeTo::ChasteSourceRoot);
    if (!component_dir.IsDir())
    {
        EXCEPTION("Unable to convert CellML model: required Chaste component '" << mComponentName
                  << "' does not exist in '" << ChasteBuildRootDir() << "'.");
    }
    // Try the conversion
    try
    {
        // Need to create a .so file from the CellML...
        if (PetscTools::AmMaster())
        {
            // Create a temporary folder within heart/dynamic
            std::stringstream folder_name;
            folder_name << "dynamic/tmp_" << getpid() << "_" << time(NULL);
            tmp_folder = component_dir.GetAbsolutePath() + "/" + folder_name.str();
            build_folder = component_dir.GetAbsolutePath() + "/build/" + ChasteBuildDirName() + "/" + folder_name.str();
            int ret = mkdir(tmp_folder.c_str(), 0700);
            if (ret != 0)
            {
                EXCEPTION("Failed to create temporary folder '" << tmp_folder << "' for CellML conversion: "
                          << strerror(errno));
            }
            // Copy the .cellml file (and any relevant others) into the temporary folder
            size_t dot_pos = rCellmlFullPath.rfind('.');
            std::string cellml_base = rCellmlFullPath.substr(0, dot_pos);
            EXPECT0(system, "cp " + cellml_base + "* " + tmp_folder);
            // If there's a config file, copy that too
            std::string config_path = rCellmlFullPath.substr(0, rCellmlFullPath.length() - 7) + "-conf.xml";
            if (FileFinder(config_path, RelativeTo::Absolute).Exists())
            {
                EXPECT0(system, "cp " + config_path + " " + tmp_folder);
            }
            // Change to Chaste source folder
            EXPECT0(chdir, ChasteBuildRootDir());
            // Run scons to generate C++ code and compile it to a .so
            EXPECT0(system, "scons --warn=no-all dyn_libs_only=1 build=" + ChasteBuildType() + " " + tmp_folder);
            EXCEPT_IF_NOT(FileFinder(tmp_folder + "/lib" + rModelLeafName + "so", RelativeTo::Absolute).Exists());
            // CD back
            EXPECT0(chdir, old_cwd);
            // Copy the .so to the same folder as the original .cellml file
            EXPECT0(system, "cp " + tmp_folder + "/lib" + rModelLeafName + "so " + rCellmlFolder);
            if (mPreserveGeneratedSources)
            {
                // Copy generated source code as well
                EXPECT0(system, "cp " + build_folder + "/*.?pp " + rCellmlFolder);
            }
            // Delete the temporary folders
            EXPECT0(system, "rm -r " + build_folder);
            EXPECT0(system, "rm -r " + tmp_folder);
        }
    }
    catch (Exception& e)
    {
        PetscTools::ReplicateException(true);
        if (FileFinder(tmp_folder, RelativeTo::Absolute).Exists())
        {
            if (mPreserveGeneratedSources)
            {
                // Copy any temporary files
                IGNORE_RET(system, "cp -r " + build_folder + " " + rCellmlFolder + "/build/");
                IGNORE_RET(system, "cp -r " + tmp_folder + " " + rCellmlFolder + "/tmp/");
            }
            // Delete the temporary folders
            IGNORE_RET(system, "rm -rf " + build_folder); // -f because folder might not exist
            IGNORE_RET(system, "rm -r " + tmp_folder);
        }
        IGNORE_RET(chdir, old_cwd);
        EXCEPTION("Conversion of CellML to Chaste shared object failed. Error was: " + e.GetMessage());
    }
    // This also has the effect of a barrier, ensuring all processes wait for the
    // shared library to be created.
    PetscTools::ReplicateException(false);
}

void CellMLToSharedLibraryConverter::CreateOptionsFile(const OutputFileHandler& rHandler,
                                                       const std::string& rModelName,
                                                       const std::vector<std::string>& rArgs,
                                                       const std::string& rExtraXml)
{
    if (PetscTools::AmMaster())
    {
        out_stream p_optfile = rHandler.OpenOutputFile(rModelName + "-conf.xml");
        (*p_optfile) << "<?xml version='1.0'?>" << std::endl
                     << "<pycml_config>" << std::endl;
        if (!rArgs.empty())
        {
            (*p_optfile) << "<command_line_args>" << std::endl;
            for (unsigned i=0; i<rArgs.size(); i++)
            {
                (*p_optfile) << "<arg>" << rArgs[i] << "</arg>" << std::endl;
            }
            (*p_optfile) << "</command_line_args>" << std::endl;
        }
        (*p_optfile) << rExtraXml << "</pycml_config>" << std::endl;
        p_optfile->close();
    }
    PetscTools::Barrier("CellMLToSharedLibraryConverter::CreateOptionsFile");
}
