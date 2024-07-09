#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "FileFinder.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "OutputFileHandler.hpp"

#include "OutputFileHandler.cppwg.hpp"

namespace py = pybind11;
typedef OutputFileHandler OutputFileHandler;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_OutputFileHandler_class(py::module &m){
py::class_<OutputFileHandler  , boost::shared_ptr<OutputFileHandler >   >(m, "OutputFileHandler")
        .def(py::init<::std::string const &, bool >(), py::arg("rDirectory"), py::arg("cleanOutputDirectory") = true)
        .def(py::init<::FileFinder const &, bool >(), py::arg("rDirectory"), py::arg("cleanOutputDirectory") = true)
        .def_static(
            "GetChasteTestOutputDirectory",
            (::std::string(*)()) &OutputFileHandler::GetChasteTestOutputDirectory,
            " "  )
        .def(
            "GetOutputDirectoryFullPath",
            (::std::string(OutputFileHandler::*)() const ) &OutputFileHandler::GetOutputDirectoryFullPath,
            " "  )
        .def(
            "GetRelativePath",
            (::std::string(OutputFileHandler::*)() const ) &OutputFileHandler::GetRelativePath,
            " "  )
        .def(
            "SetArchiveDirectory",
            (void(OutputFileHandler::*)() const ) &OutputFileHandler::SetArchiveDirectory,
            " "  )
        .def(
            "CopyFileTo",
            (::FileFinder(OutputFileHandler::*)(::FileFinder const &) const ) &OutputFileHandler::CopyFileTo,
            " " , py::arg("rSourceFile") )
        .def(
            "FindFile",
            (::FileFinder(OutputFileHandler::*)(::std::string) const ) &OutputFileHandler::FindFile,
            " " , py::arg("leafName") )
    ;
}
