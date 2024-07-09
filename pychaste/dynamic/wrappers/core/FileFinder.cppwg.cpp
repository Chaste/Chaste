#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FileFinder.hpp"

#include "FileFinder.cppwg.hpp"

namespace py = pybind11;
typedef FileFinder FileFinder;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class FileFinder_Overrides : public FileFinder{
    public:
    using FileFinder::FileFinder;
    void SetPath(::std::string const & rPath, ::RelativeTo::Value relativeTo) override {
        PYBIND11_OVERRIDE(
            void,
            FileFinder,
            SetPath,
                    rPath,
        relativeTo);
    }
    void SetPath(::std::string const & rLeafName, ::FileFinder const & rParentOrSibling) override {
        PYBIND11_OVERRIDE(
            void,
            FileFinder,
            SetPath,
                    rLeafName,
        rParentOrSibling);
    }

};
void register_FileFinder_class(py::module &m){
py::class_<FileFinder , FileFinder_Overrides , boost::shared_ptr<FileFinder >   >(m, "FileFinder")
        .def(py::init< >())
        .def(py::init<::std::string const &, ::RelativeTo::Value >(), py::arg("rPath"), py::arg("relativeTo"))
        .def(py::init<::std::string const &, ::FileFinder const & >(), py::arg("rLeafName"), py::arg("rParentOrSibling"))
        .def(py::init<::std::filesystem::path const & >(), py::arg("rPath"))
        .def(
            "SetPath",
            (void(FileFinder::*)(::std::string const &, ::RelativeTo::Value)) &FileFinder::SetPath,
            " " , py::arg("rPath"), py::arg("relativeTo") )
        .def(
            "SetPath",
            (void(FileFinder::*)(::std::string const &, ::FileFinder const &)) &FileFinder::SetPath,
            " " , py::arg("rLeafName"), py::arg("rParentOrSibling") )
        .def(
            "IsPathSet",
            (bool(FileFinder::*)() const ) &FileFinder::IsPathSet,
            " "  )
        .def(
            "Exists",
            (bool(FileFinder::*)() const ) &FileFinder::Exists,
            " "  )
        .def(
            "IsFile",
            (bool(FileFinder::*)() const ) &FileFinder::IsFile,
            " "  )
        .def(
            "IsDir",
            (bool(FileFinder::*)() const ) &FileFinder::IsDir,
            " "  )
        .def(
            "IsEmpty",
            (bool(FileFinder::*)() const ) &FileFinder::IsEmpty,
            " "  )
        .def(
            "GetAbsolutePath",
            (::std::string(FileFinder::*)() const ) &FileFinder::GetAbsolutePath,
            " "  )
        .def(
            "IsNewerThan",
            (bool(FileFinder::*)(::FileFinder const &) const ) &FileFinder::IsNewerThan,
            " " , py::arg("rOtherEntity") )
        .def(
            "GetLeafName",
            (::std::string(FileFinder::*)() const ) &FileFinder::GetLeafName,
            " "  )
        .def(
            "GetLeafNameNoExtension",
            (::std::string(FileFinder::*)() const ) &FileFinder::GetLeafNameNoExtension,
            " "  )
        .def(
            "GetExtension",
            (::std::string(FileFinder::*)() const ) &FileFinder::GetExtension,
            " "  )
        .def(
            "GetParent",
            (::FileFinder(FileFinder::*)() const ) &FileFinder::GetParent,
            " "  )
        .def(
            "GetRelativePath",
            (::std::string(FileFinder::*)(::FileFinder const &) const ) &FileFinder::GetRelativePath,
            " " , py::arg("rBasePath") )
        .def(
            "CopyTo",
            (::FileFinder(FileFinder::*)(::FileFinder const &) const ) &FileFinder::CopyTo,
            " " , py::arg("rDest") )
        .def(
            "Remove",
            (void(FileFinder::*)() const ) &FileFinder::Remove,
            " "  )
        .def(
            "DangerousRemove",
            (void(FileFinder::*)() const ) &FileFinder::DangerousRemove,
            " "  )
        .def(
            "FindMatches",
            (::std::vector<FileFinder>(FileFinder::*)(::std::string const &) const ) &FileFinder::FindMatches,
            " " , py::arg("rPattern") )
        .def_static(
            "IsAbsolutePath",
            (bool(*)(::std::string const &)) &FileFinder::IsAbsolutePath,
            " " , py::arg("rPath") )
        .def_static(
            "ReplaceSpacesWithUnderscores",
            (void(*)(::std::string &)) &FileFinder::ReplaceSpacesWithUnderscores,
            " " , py::arg("rPath") )
        .def_static(
            "ReplaceUnderscoresWithSpaces",
            (void(*)(::std::string &)) &FileFinder::ReplaceUnderscoresWithSpaces,
            " " , py::arg("rPath") )
        .def_static(
            "FakePath",
            (void(*)(::RelativeTo::Value, ::std::string const &)) &FileFinder::FakePath,
            " " , py::arg("fakeWhat"), py::arg("rFakePath") )
        .def_static(
            "StopFaking",
            (void(*)()) &FileFinder::StopFaking,
            " "  )
    ;
}
