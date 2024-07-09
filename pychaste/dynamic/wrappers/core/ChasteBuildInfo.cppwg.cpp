#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Version.hpp"

#include "ChasteBuildInfo.cppwg.hpp"

namespace py = pybind11;
typedef ChasteBuildInfo ChasteBuildInfo;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_ChasteBuildInfo_class(py::module &m){
py::class_<ChasteBuildInfo  , boost::shared_ptr<ChasteBuildInfo >   >(m, "ChasteBuildInfo")
        .def(py::init< >())
        .def_static(
            "GetLicenceText",
            (::std::string(*)()) &ChasteBuildInfo::GetLicenceText,
            " "  )
        .def_static(
            "GetRootDir",
            (char const *(*)()) &ChasteBuildInfo::GetRootDir,
            " "  , py::return_value_policy::reference)
        .def_static(
            "GetVersionString",
            (::std::string(*)()) &ChasteBuildInfo::GetVersionString,
            " "  )
        .def_static(
            "GetMajorReleaseNumber",
            (unsigned int(*)()) &ChasteBuildInfo::GetMajorReleaseNumber,
            " "  )
        .def_static(
            "GetMinorReleaseNumber",
            (unsigned int(*)()) &ChasteBuildInfo::GetMinorReleaseNumber,
            " "  )
        .def_static(
            "GetRevisionNumber",
            (long long unsigned int(*)()) &ChasteBuildInfo::GetRevisionNumber,
            " "  )
        .def_static(
            "IsWorkingCopyModified",
            (bool(*)()) &ChasteBuildInfo::IsWorkingCopyModified,
            " "  )
        .def_static(
            "GetBuildTime",
            (char const *(*)()) &ChasteBuildInfo::GetBuildTime,
            " "  , py::return_value_policy::reference)
        .def_static(
            "GetCurrentTime",
            (char const *(*)()) &ChasteBuildInfo::GetCurrentTime,
            " "  , py::return_value_policy::reference)
        .def_static(
            "GetBuilderUnameInfo",
            (char const *(*)()) &ChasteBuildInfo::GetBuilderUnameInfo,
            " "  , py::return_value_policy::reference)
        .def_static(
            "GetBuildInformation",
            (char const *(*)()) &ChasteBuildInfo::GetBuildInformation,
            " "  , py::return_value_policy::reference)
        .def_static(
            "GetCompilerType",
            (char const *(*)()) &ChasteBuildInfo::GetCompilerType,
            " "  , py::return_value_policy::reference)
        .def_static(
            "GetCompilerVersion",
            (char const *(*)()) &ChasteBuildInfo::GetCompilerVersion,
            " "  , py::return_value_policy::reference)
        .def_static(
            "GetCompilerFlags",
            (char const *(*)()) &ChasteBuildInfo::GetCompilerFlags,
            " "  , py::return_value_policy::reference)
        .def_static(
            "GetXsdVersion",
            (char const *(*)()) &ChasteBuildInfo::GetXsdVersion,
            " "  , py::return_value_policy::reference)
        .def_static(
            "rGetProjectVersions",
            (::std::map<std::basic_string<char>, std::basic_string<char>> const &(*)()) &ChasteBuildInfo::rGetProjectVersions,
            " "  , py::return_value_policy::reference_internal)
        .def_static(
            "rGetIfProjectsModified",
            (::std::map<std::basic_string<char>, std::basic_string<char>> const &(*)()) &ChasteBuildInfo::rGetIfProjectsModified,
            " "  , py::return_value_policy::reference_internal)
        .def_static(
            "GetProvenanceString",
            (::std::string(*)()) &ChasteBuildInfo::GetProvenanceString,
            " "  )
        .def_static(
            "GetChasteCodegenVersion",
            (::std::string(*)()) &ChasteBuildInfo::GetChasteCodegenVersion,
            " "  )
    ;
}
