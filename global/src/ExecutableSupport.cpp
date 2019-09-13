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

#include "ExecutableSupport.hpp"

#include "Version.hpp"

#include <iostream>
#include <sstream>
#ifdef _MSC_VER
#include <windows.h> // For system info on Windows
//GetCurrentTime clashes with a similarly-named API in <windows.h>
#define GetCurrentTime() ChasteBuildInfo::GetCurrentTime()
#define ChasteGetCurrentTime() GetCurrentTime()
#else
#include <sys/utsname.h> // For uname
#define ChasteGetCurrentTime() ChasteBuildInfo::GetCurrentTime()
#endif
#include <hdf5.h>
#include <parmetis.h>


#include "ChasteSerialization.hpp"
#include "CommandLineArguments.hpp"
#include "PetscException.hpp"
#include "PetscSetupUtils.hpp"
#include "PetscTools.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

#ifdef CHASTE_CVODE
#include <sundials/sundials_config.h>
#if CHASTE_SUNDIALS_VERSION >= 20600
#if CHASTE_SUNDIALS_VERSION >= 30000
// SUNDIALS 3.0 upwards uses SUNDIALS_VERSION instead of SUNDIALS_PACKAGE_VERSION.
#define CHASTE_SUNDIALS_PACKAGE_VERSION SUNDIALS_VERSION
#else
// SUNDIALS 2.6 upwards defines SUNDIALS_PACKAGE_VERSION with quotes...
#include <boost/preprocessor/stringize.hpp>
#define CHASTE_SUNDIALS_PACKAGE_VERSION BOOST_PP_STRINGIZE(SUNDIALS_PACKAGE_VERSION)
#endif // SUNDIALS >= 3.0.0
#else
#define CHASTE_SUNDIALS_PACKAGE_VERSION SUNDIALS_PACKAGE_VERSION
#endif // SUNDIALS >= 2.6.0
#endif

// Note: the following are not a requirement for cell-based Chaste so may not be present!
//#include <xsd/cxx/version.hxx>
#ifdef CHASTE_XERCES
#include <xercesc/util/XercesVersion.hpp>
#endif

// Check whether the version of ParMETIS being used is the one we wanted
#ifdef CHASTE_PARMETIS_REQUIRED
#if PARMETIS_MAJOR_VERSION != CHASTE_PARMETIS_REQUIRED
#error "Wrong ParMETIS version found: " #CHASTE_PARMETIS_REQUIRED " requested but " #PARMETIS_MAJOR_VERSION " present"
#endif
#endif

FileFinder ExecutableSupport::mOutputDirectory;

void ExecutableSupport::SetOutputDirectory(const std::string& rOutputDirectory)
{
    if (FileFinder::IsAbsolutePath(rOutputDirectory))
    {
        mOutputDirectory.SetPath(rOutputDirectory, RelativeTo::Absolute);
    }
    else
    {
        mOutputDirectory.SetPath(rOutputDirectory, RelativeTo::ChasteTestOutput);
    }
}

void ExecutableSupport::InitializePetsc(int* pArgc, char*** pArgv)
{
    // Store the arguments in case other code needs them
    CommandLineArguments::Instance()->p_argc = pArgc;
    CommandLineArguments::Instance()->p_argv = pArgv;
    // Initialise PETSc
    PETSCEXCEPT(PetscInitialize(pArgc, pArgv, PETSC_NULL, PETSC_NULL));
    // Set default output folder
    if (!mOutputDirectory.IsPathSet())
    {
        // LCOV_EXCL_START
        //depends on order of calls.  Extract to method?
        mOutputDirectory.SetPath("", RelativeTo::ChasteTestOutput);
        // LCOV_EXCL_STOP
    }
}

void ExecutableSupport::ShowCopyright()
{
    // Compilation information
    std::stringstream provenance_msg;
    provenance_msg << "This version of Chaste was compiled on:\n";
    provenance_msg << ChasteBuildInfo::GetBuildTime() << " by " << ChasteBuildInfo::GetBuilderUnameInfo() << " (uname)\n";
    provenance_msg << "from revision number " << std::hex << ChasteBuildInfo::GetRevisionNumber() << std::dec << " with build type " << ChasteBuildInfo::GetBuildInformation() << ".\n\n";

    // Only show one copy of copyright/header
    if (PetscTools::AmMaster())
    {
        std::cout << ChasteBuildInfo::GetLicenceText() << std::endl;

        // Write provenance information to stdout
        std::cout << provenance_msg.str() << std::flush;
    }
}

void ExecutableSupport::ShowParallelLaunching()
{
    if (PetscTools::IsParallel())
    {
        // Information to show that Chaste is being run in parallel
        PetscTools::BeginRoundRobin();
        std::cout << "Chaste launched on process " << PetscTools::GetMyRank()
                  << " of " << PetscTools::GetNumProcs() << "." << std::endl
                  << std::flush;
        PetscTools::EndRoundRobin();
    }
}

void ExecutableSupport::WriteMachineInfoFile(std::string fileBaseName)
{
    if (!mOutputDirectory.IsPathSet())
    {
        // LCOV_EXCL_START
        //depends on order of calls.  Extract to method?
        mOutputDirectory.SetPath("", RelativeTo::ChasteTestOutput);
        // LCOV_EXCL_STOP
    }
    OutputFileHandler out_file_handler(mOutputDirectory, false);
    std::stringstream file_name;
    file_name << fileBaseName << "_" << PetscTools::GetMyRank() << ".txt";
    out_stream out_file = out_file_handler.OpenOutputFile(file_name.str());
    *out_file << "Process " << PetscTools::GetMyRank() << " of "
              << PetscTools::GetNumProcs() << "." << std::endl
              << std::flush;

#ifdef _MSC_VER
    //use native windows equivalent of system info. from uname and /proc/cpuinfo

    //operating system and machine name
    *out_file << "uname sysname  = "
              << "Microsoft Windows" << std::endl;
#define INFO_BUFFER_SIZE 32767
    TCHAR info_buffer[INFO_BUFFER_SIZE];
    DWORD buffer_char_count = INFO_BUFFER_SIZE;
    if (!GetComputerName(info_buffer, &buffer_char_count))
        *out_file << "uname nodename = "
                  << "Windows machine name is unknown" << std::endl;
    else
        *out_file << "uname nodename = " << info_buffer << std::endl;
    //more detailed operating system version information
    OSVERSIONINFOEX os_info;
    ZeroMemory(&os_info, sizeof(OSVERSIONINFO));
    os_info.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    //copy os version info into the structure
    GetVersionEx((OSVERSIONINFO*)&os_info);
    //Pivot around Windows Vista (dwMajorVersion >= 6)
    //See http://msdn.microsoft.com/en-us/library/ms724834%28v=vs.85%29.aspx for details
    if (os_info.dwMajorVersion < 6)
    { //earlier than Windows Vista
        *out_file << "uname release  = "
                  << "Microsoft Windows Server 2003 R2 (or earlier)" << std::endl;
    }
    else
    {
        //reverse chronological order (simply add newer OS version to the top)
        if (os_info.dwMajorVersion > 6)
        {
            *out_file << "uname release  = "
                      << "Microsoft Windows (Later than Microsoft Windows 8)" << std::endl;
        }
        else //os_info.dwMajorVersion == 6
        {
            if (os_info.dwMinorVersion == 2)
            {
                if (os_info.wProductType == VER_NT_WORKSTATION)
                    *out_file << "uname release  = "
                              << "Microsoft Windows 8" << std::endl;
                else
                    *out_file << "uname release  = "
                              << "Microsoft Windows Server 2012" << std::endl;
            }
            else if (os_info.dwMinorVersion == 1)
            {
                if (os_info.wProductType == VER_NT_WORKSTATION)
                    *out_file << "uname release  = "
                              << "Microsoft Windows 7" << std::endl;
                else
                    *out_file << "uname release  = "
                              << "Microsoft Windows Server 2008 R2" << std::endl;
            }
            else if (os_info.dwMinorVersion == 0)
            {
                if (os_info.wProductType == VER_NT_WORKSTATION)
                    *out_file << "uname release  = "
                              << "Microsoft Windows Server 2008" << std::endl;
                else
                    *out_file << "uname release  = "
                              << "Microsoft Windows Vista" << std::endl;
            }
        }
    }
    *out_file << "uname version  = " << os_info.dwMajorVersion << "." << os_info.dwMinorVersion << std::endl;

    //hardware information
    // See http://msdn.microsoft.com/en-us/library/ms724958%28v=vs.85%29.aspx for details of the SYSTEM_INFO structure
    SYSTEM_INFO sys_info;
    GetSystemInfo(&sys_info);
    switch (sys_info.wProcessorArchitecture)
    {
        case PROCESSOR_ARCHITECTURE_AMD64:
            *out_file << "uname machine  = "
                      << "x64 (AMD or Intel)" << std::endl;
            break;
        case PROCESSOR_ARCHITECTURE_ARM:
            *out_file << "uname machine  = "
                      << "ARM" << std::endl;
            break;
        case PROCESSOR_ARCHITECTURE_IA64:
            *out_file << "uname machine  = "
                      << "Intel Itanium-based" << std::endl;
            break;
        case PROCESSOR_ARCHITECTURE_INTEL:
            *out_file << "uname machine  = "
                      << "x86" << std::endl;
            break;
        case PROCESSOR_ARCHITECTURE_UNKNOWN:
            *out_file << "uname machine  = "
                      << "Unknown Architecture" << std::endl;
            break;
        default:
            *out_file << "uname machine  = "
                      << "Other Architecture. Code = " << sys_info.wProcessorArchitecture << std::endl;
            break;
    }

    *out_file << "\nInformation on number and type of processors:" << std::endl;
    *out_file << sys_info.dwNumberOfProcessors;
    *out_file << "\nInformation on processor caches, in the same order as above:" << std::endl;
    *out_file << "Unknown" << std::endl; ///\todo #2016 Get CPU info on Windows?

    //get physical memory
    MEMORYSTATUSEX mem_status;
    mem_status.dwLength = sizeof(mem_status);
    GlobalMemoryStatusEx(&mem_status);
    *out_file << "\nInformation on system memory:" << std::endl;
    *out_file << mem_status.ullTotalPhys / 1024 << " kB" << std::endl;
#else
    struct utsname uts_info;
    uname(&uts_info);

    *out_file << "uname sysname  = " << uts_info.sysname << std::endl
              << std::flush;
    *out_file << "uname nodename = " << uts_info.nodename << std::endl
              << std::flush;
    *out_file << "uname release  = " << uts_info.release << std::endl
              << std::flush;
    *out_file << "uname version  = " << uts_info.version << std::endl
              << std::flush;
    *out_file << "uname machine  = " << uts_info.machine << std::endl
              << std::flush;
    char buffer[100];
    FILE* system_info;

#ifdef __APPLE__
    *out_file << "\nInformation on number and type processors, and cache and memory sizes (in bytes)\n";
    system_info = popen("sysctl hw.ncpu hw.physicalcpu machdep.cpu.brand_string hw.l1icachesize hw.l1dcachesize hw.l2cachesize hw.l3cachesize hw.memsize", "r");
    while (fgets(buffer, 100, system_info) != NULL)
    {
        *out_file << buffer;
    }
    fclose(system_info);
#else
    //GNU
    *out_file << "\nInformation on number and type of processors:\n";
    system_info = popen("grep ^model.name /proc/cpuinfo", "r");
    while (fgets(buffer, 100, system_info) != nullptr)
    {
        *out_file << buffer;
    }
    fclose(system_info);

    *out_file << "\nInformation on processor caches, in the same order as above:\n";
    system_info = popen("grep ^cache.size /proc/cpuinfo", "r");
    while (fgets(buffer, 100, system_info) != nullptr)
    {
        *out_file << buffer;
    }
    fclose(system_info);

    *out_file << "\nInformation on system memory:\n";
    system_info = popen("grep ^MemTotal /proc/meminfo", "r");
    while (fgets(buffer, 100, system_info) != nullptr)
    {
        *out_file << buffer;
    }
    fclose(system_info);
#endif //end of __APPLE__ not defined
#endif //end of _MSC_VER not defined

    out_file->close();
}

void ExecutableSupport::WriteProvenanceInfoFile()
{
    if (!mOutputDirectory.IsPathSet())
    {
        // LCOV_EXCL_START
        //depends on order of calls.  Extract to method?
        mOutputDirectory.SetPath("", RelativeTo::ChasteTestOutput);
        // LCOV_EXCL_STOP
    }
    OutputFileHandler out_file_handler(mOutputDirectory, false);
    out_stream out_file = out_file_handler.OpenOutputFile("provenance_info_", PetscTools::GetMyRank(), ".txt");

    // Compilation information
    std::stringstream provenance_msg;
    provenance_msg << "This version of Chaste was compiled on:\n";
    provenance_msg << ChasteBuildInfo::GetBuildTime() << " by " << ChasteBuildInfo::GetBuilderUnameInfo() << " (uname)\n";
    provenance_msg << "from revision number " << std::hex << ChasteBuildInfo::GetRevisionNumber() << std::dec << " with build type " << ChasteBuildInfo::GetBuildInformation() << ".\n\n";
    *out_file << provenance_msg.str();

    std::string output;
    GetBuildInfo(output);
    *out_file << output;

    out_file->close();
}

void ExecutableSupport::GetBuildInfo(std::string& rInfo)
{
    std::stringstream output;
    output << "<ChasteBuildInfo>\n";

    output << "\t<ProvenanceInfo>\n";
    output << "\t\t<VersionString>" << ChasteBuildInfo::GetVersionString() << "</VersionString> <!-- build specific -->\n";
    output << "\t\t<IsWorkingCopyModified>" << ChasteBuildInfo::IsWorkingCopyModified() << "</IsWorkingCopyModified>\n";
    output << "\t\t<BuildInformation>" << ChasteBuildInfo::GetBuildInformation() << "</BuildInformation>\n";
    output << "\t\t<BuildTime>" << ChasteBuildInfo::GetBuildTime() << "</BuildTime>\n";
    output << "\t\t<CurrentTime>" << ChasteGetCurrentTime() << "</CurrentTime>\n";
    output << "\t\t<BuilderUnameInfo>" << ChasteBuildInfo::GetBuilderUnameInfo() << "</BuilderUnameInfo>\n";

    output << "\t\t<Projects>\n";
    {
        const std::map<std::string, std::string>& r_projects_modified = ChasteBuildInfo::rGetIfProjectsModified();
        const std::map<std::string, std::string>& r_projects_versions = ChasteBuildInfo::rGetProjectVersions();
        for (const auto& r_project_version : r_projects_versions)
        {
            // LCOV_EXCL_START
            // No projects are checked out for continuous builds normally!
            output << "\t\t\t<Project>" << std::endl;
            output << "\t\t\t\t<Name>" << r_project_version.first << "</Name>" << std::endl;
            output << "\t\t\t\t<Version>" << r_project_version.second << "</Version>" << std::endl;
            output << "\t\t\t\t<Modified>" << r_projects_modified.at(r_project_version.first) << "</Modified>" << std::endl;
            output << "\t\t\t</Project>" << std::endl;
            // LCOV_EXCL_STOP
        }
    }
    output << "\t\t</Projects>\n";

    output << "\t</ProvenanceInfo>\n";

    output << "\t<Compiler>\n";
    output << "\t\t<NameAndVersion>" << ChasteBuildInfo::GetCompilerType() << ", version " << ChasteBuildInfo::GetCompilerVersion() << "</NameAndVersion>\n";
    output << "\t\t<Flags>" << ChasteBuildInfo::GetCompilerFlags() << "</Flags>\n";
    output << "\t</Compiler>\n";

    output << "\t<Libraries>\n";

    output << "\t\t<CompiledIn>\n";
    output << "\t\t\t<PETSc>" << PETSC_VERSION_MAJOR << "." << PETSC_VERSION_MINOR << "." << PETSC_VERSION_SUBMINOR << "</PETSc>\n";
    output << "\t\t\t<Boost>" << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << "</Boost>\n";
    output << "\t\t\t<HDF5>" << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE << "</HDF5>\n";
    output << "\t\t\t<Parmetis>" << PARMETIS_MAJOR_VERSION << "." << PARMETIS_MINOR_VERSION;
#ifdef PARMETIS_SUBMINOR_VERSION // they only added this in v4.? !!
    output << "." << PARMETIS_SUBMINOR_VERSION;
#endif
    output << "</Parmetis>" << std::endl;
    output << "\t\t</CompiledIn>\n";

    output << "\t\t<Binaries>\n";
    output << "\t\t\t<XSD>" << ChasteBuildInfo::GetXsdVersion() << "</XSD>\n";
    output << "\t\t</Binaries>\n";

    output << "\t\t<Optional>\n";
#ifdef CHASTE_CVODE
    output << "\t\t\t<SUNDIALS>" << CHASTE_SUNDIALS_PACKAGE_VERSION << "</SUNDIALS>";
#if CHASTE_SUNDIALS_VERSION < 30000
    output << "<!-- includes Cvode of a different version number -->";
#endif
    output << std::endl;
#else
    output << "\t\t\t<SUNDIALS>no</SUNDIALS>\n";
#endif
#ifdef CHASTE_VTK
    output << "\t\t\t<VTK>" << VTK_MAJOR_VERSION << "." << VTK_MINOR_VERSION << "</VTK>\n";
#else
    output << "\t\t\t<VTK>no</VTK>\n";
#endif
#ifdef CHASTE_XERCES
    output << "\t\t\t<Xerces>" << XERCES_FULLVERSIONDOT << "</Xerces>\n"; // Note: not a requirement for cell-based so may not be present!
#else
    output << "\t\t\t<Xerces>no</Xerces>\n";
#endif
    output << "\t\t</Optional>\n";

    output << "\t</Libraries>\n";

    output << "</ChasteBuildInfo>" << std::endl;
    rInfo = output.str();
}

void ExecutableSupport::StandardStartup(int* pArgc, char*** pArgv)
{
    InitializePetsc(pArgc, pArgv);
    ShowCopyright();
    ShowParallelLaunching();
}

void ExecutableSupport::StartupWithoutShowingCopyright(int* pArgc, char*** pArgv)
{
    InitializePetsc(pArgc, pArgv);
    ShowParallelLaunching();
}

void ExecutableSupport::PrintError(const std::string& rMessage, bool masterOnly)
{
    if (!masterOnly || PetscTools::AmMaster())
    {
        // Write the error message to stderr
        std::cerr << rMessage << std::endl;
    }

    // Write the error message to file
    if (!mOutputDirectory.IsPathSet())
    {
        // LCOV_EXCL_START
        //depends on order of calls.  Extract to method?
        mOutputDirectory.SetPath("", RelativeTo::ChasteTestOutput);
        // LCOV_EXCL_STOP
    }
    OutputFileHandler out_file_handler(mOutputDirectory, false);
    out_stream out_file = out_file_handler.OpenOutputFile("chaste_errors_", PetscTools::GetMyRank(), ".txt", std::ios::out | std::ios::app);
    *out_file << rMessage << std::endl;
    out_file->close();
}

void ExecutableSupport::Print(const std::string& rMessage)
{
    if (PetscTools::AmMaster())
    {
        // Write the error message to stdout
        std::cout << rMessage << std::endl
                  << std::flush;
    }
}

void ExecutableSupport::FinalizePetsc()
{
    // Causes memory failure (and seg fault) in PETSc 3.2 with MPICH-1
    PetscSetupUtils::CommonFinalize();
}
