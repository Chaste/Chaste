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

#include "CellMLToSharedLibraryConverter.hpp"

#include <sstream>
#include <fstream>      // for std::ofstream
#include <sys/stat.h> // For mkdir()
#include <ctime>
#include <cstring> // For strerror()
#include <cerrno> // For errno

#include <boost/foreach.hpp>

#include "ChasteSyscalls.hpp"
#include "Exception.hpp"
#include "Warnings.hpp"
#include "ChasteBuildRoot.hpp"
#include "PetscTools.hpp"
#include "DynamicModelLoaderRegistry.hpp"
#include "GetCurrentWorkingDirectory.hpp"

#define IGNORE_EXCEPTIONS(code) \
    try {                       \
        code;                   \
    } catch (...) {}

/** Set the .so suffix */
#ifdef __APPLE__
//Mac OSX
const std::string CellMLToSharedLibraryConverter::msSoSuffix = "dylib";
#else
//Normal behaviour
const std::string CellMLToSharedLibraryConverter::msSoSuffix = "so";
#endif

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
    // Find out whether rFilePath is a .cellml or .so
    size_t dot_position = absolute_path.find_last_of(".");
    if (dot_position == std::string::npos || absolute_path[dot_position+1] =='/')
    {
        // Either there is no dot in the absolute path or there is one along the path (/home/user/chaste/../temp/stuff)
        EXCEPTION("File does not have an extension: " + absolute_path);
    }
    std::string extension = absolute_path.substr(dot_position+1);
    // We make a modifiable version of the const FileFinder just incase we feel like
    // amending the suffix
    FileFinder file_path_copy(rFilePath);
#ifdef __APPLE__
    if (extension == "so")
    {
        WARN_ONCE_ONLY("CellMLToSharedLibraryConverter asked to load a \".so\" file.  On this architecture it should be \".dylib\"");
        extension = "dylib";
        absolute_path.replace(dot_position+1, 5, extension);
        file_path_copy.SetPath(absolute_path,  RelativeTo::Absolute);
    }
#endif
    // Check the file exists
    if (!file_path_copy.Exists())
    {
        EXCEPTION("Dynamically loadable cell model '" + absolute_path + "' does not exist.");
    }
    if (extension == "cellml")
    {
        // Split the path into folder and leaf
        size_t slash_position = absolute_path.find_last_of("/\\");
        assert(slash_position != std::string::npos);
        std::string folder = absolute_path.substr(0, slash_position+1); // Include trailing slash
        std::string leaf = absolute_path.substr(slash_position+1, dot_position-slash_position); // Include dot
        std::string so_path = folder + "lib" + leaf + msSoSuffix;
        // Does the .so file already exist (and was it modified after the .cellml?)
        FileFinder so_file(so_path, RelativeTo::Absolute);
        if (!so_file.Exists() || rFilePath.IsNewerThan(so_file))
        {
            if (!isCollective)
            {
                EXCEPTION("Unable to convert .cellml to .so unless called collectively, due to possible race conditions.");
            }
            ConvertCellmlToSo(absolute_path, folder);
        }
        // Load the .so
        p_loader = DynamicModelLoaderRegistry::Instance()->GetLoader(so_file);
    }
    else if (extension == msSoSuffix)
    {
        // Just load the .so
        // Note that this may have been modified to .dylib
        p_loader = DynamicModelLoaderRegistry::Instance()->GetLoader(file_path_copy);
    }
    else
    {
        EXCEPTION("Unsupported extension '." + extension + "' of file '" + absolute_path + "'; must be .so, .dylib or .cellml");
    }

    return p_loader;
}

void CellMLToSharedLibraryConverter::ConvertCellmlToSo(const std::string& rCellmlFullPath,
                                                       const std::string& rCellmlFolder)
{
    FileFinder tmp_folder;
    FileFinder build_folder;

    std::string old_cwd = GetCurrentWorkingDirectory();
    // Check that the Chaste build tree exists
    FileFinder chaste_root("", RelativeTo::ChasteBuildRoot);

    if (!chaste_root.IsDir())
    {
        EXCEPTION("No Chaste build tree found at '" << chaste_root.GetAbsolutePath()
                  << "' - you need the source to use CellML models directly in Chaste.");
    }
    FileFinder component_dir(mComponentName, RelativeTo::ChasteBuildRoot);
    if (!component_dir.IsDir())
    {
        EXCEPTION("Unable to convert CellML model: required Chaste component '" << mComponentName
                  << "' does not exist in '" << chaste_root.GetAbsolutePath() << "'.");
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

#ifdef CHASTE_CMAKE
            tmp_folder.SetPath(component_dir.GetAbsolutePath() + "/" + folder_name.str(), RelativeTo::Absolute);
            build_folder.SetPath(component_dir.GetAbsolutePath() + "/" + folder_name.str(), RelativeTo::Absolute);
#else
            tmp_folder.SetPath(component_dir.GetAbsolutePath() + "/" + folder_name.str(), RelativeTo::Absolute);
            build_folder.SetPath(component_dir.GetAbsolutePath() + "/build/" + ChasteBuildDirName() + "/" + folder_name.str(), RelativeTo::Absolute);
#endif

            int ret = mkdir((tmp_folder.GetAbsolutePath()).c_str(), 0700);
            if (ret != 0)
            {
                EXCEPTION("Failed to create temporary folder '" << tmp_folder.GetAbsolutePath() << "' for CellML conversion: "
                          << strerror(errno));
            }

            // Copy the .cellml file (and any relevant others) into the temporary folder
            FileFinder cellml_file(rCellmlFullPath, RelativeTo::Absolute);
            FileFinder cellml_folder = cellml_file.GetParent();
            std::string cellml_leaf_name = cellml_file.GetLeafNameNoExtension();
            std::vector<FileFinder> cellml_files = cellml_folder.FindMatches(cellml_leaf_name + "*");

            BOOST_FOREACH(const FileFinder& r_cellml_file, cellml_files)
            {
                r_cellml_file.CopyTo(tmp_folder);
            }

#ifdef CHASTE_CMAKE
            std::string cmake_lists_filename = tmp_folder.GetAbsolutePath() + "/CMakeLists.txt";
            std::ofstream cmake_lists_filestream(cmake_lists_filename.c_str());
            cmake_lists_filestream << "cmake_minimum_required(VERSION 2.8.12)\n" <<
                                      "add_compile_options(-std=c++14)\n" <<
                                      "find_package(Chaste COMPONENTS " << mComponentName << ")\n" <<
                                      "chaste_do_cellml(sources " << cellml_file.GetAbsolutePath() << " " << "ON)\n" <<
                                      "set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})\n" <<
                                      "include_directories(${Chaste_THIRD_PARTY_INCLUDE_DIRS} ${Chaste_INCLUDE_DIRS})\n" <<
                                      "add_library(" << cellml_leaf_name << " SHARED " << "${sources})\n" <<
                                      "if (${CMAKE_SYSTEM_NAME} MATCHES \"Darwin\")\n" <<
                                      "   target_link_libraries(" << cellml_leaf_name << " \"-Wl,-undefined,dynamic_lookup\")\n" <<
                                      "endif()\n"
                                      //"target_link_libraries(" << cellml_leaf_name << " ${Chaste_LIBRARIES})\n"
                                      ;
            cmake_lists_filestream.close();
            std::string cmake_args = " -DCMAKE_PREFIX_PATH=" + chaste_root.GetAbsolutePath() +
                                     " -DCMAKE_BUILD_TYPE=" + ChasteBuildType() +
                                     " -DChaste_ENABLE_TESTING=OFF" +
                                     " -DBUILD_SHARED_LIBS=ON";
            EXPECT0(chdir, tmp_folder.GetAbsolutePath());
            EXPECT0(system, "cmake" + cmake_args + " .");
            EXPECT0(system, "cmake --build . --config " + ChasteBuildType());
#else
            // Change to Chaste source folder
            EXPECT0(chdir, chaste_root.GetAbsolutePath());
            // Run scons to generate C++ code and compile it to a .so
            EXPECT0(system, "scons --warn=no-all dyn_libs_only=1 build=" + ChasteBuildType() + " " + tmp_folder.GetAbsolutePath());
            if (mPreserveGeneratedSources)
            {
                // Copy the generated source (.hpp and .cpp) to the same place as the .so file is going.
                // NB. CMake does this by default
                FileFinder destination_folder_for_sources(rCellmlFolder, RelativeTo::Absolute);
                // Copy generated source code as well
                std::vector<FileFinder> generated_files = build_folder.FindMatches("*.?pp");
                BOOST_FOREACH(const FileFinder& r_generated_file, generated_files)
                {
                    r_generated_file.CopyTo(destination_folder_for_sources);
                }
            }
#endif

            FileFinder so_file(tmp_folder.GetAbsolutePath() + "/lib" + cellml_leaf_name + "." + msSoSuffix, RelativeTo::Absolute);
            EXCEPT_IF_NOT(so_file.Exists());
            // CD back
            EXPECT0(chdir, old_cwd);

            // Copy the .so to the same folder as the original .cellml file
            FileFinder destination_folder(rCellmlFolder, RelativeTo::Absolute);
            so_file.CopyTo(destination_folder);

            // Delete the temporary folders
            build_folder.DangerousRemove();
            tmp_folder.DangerousRemove();
        }
    }
    catch (Exception& e)
    {
        PetscTools::ReplicateException(true);
        if (tmp_folder.IsPathSet() && tmp_folder.Exists())
        {
            if (mPreserveGeneratedSources)
            {
                // Copy any temporary files
                IGNORE_EXCEPTIONS(build_folder.CopyTo(FileFinder(rCellmlFolder + "/build/", RelativeTo::Absolute)));
                IGNORE_EXCEPTIONS(tmp_folder.CopyTo(FileFinder(rCellmlFolder + "/tmp/", RelativeTo::Absolute)));
            }
            // Delete the temporary folders
            IGNORE_EXCEPTIONS(build_folder.DangerousRemove()); // rm -r under source
            IGNORE_EXCEPTIONS(tmp_folder.DangerousRemove()); // rm -r under source
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
