#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include "PythonPetscObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PetscTools.hpp"

#include "PetscTools.cppwg.hpp"

namespace py = pybind11;
typedef PetscTools PetscTools;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_PetscTools_class(py::module &m){
py::class_<PetscTools  , boost::shared_ptr<PetscTools >   >(m, "PetscTools")
        .def(py::init< >())
        .def_static(
            "ResetCache",
            (void(*)()) &PetscTools::ResetCache,
            " "  )
        .def_static(
            "IsInitialised",
            (bool(*)()) &PetscTools::IsInitialised,
            " "  )
        .def_static(
            "IsSequential",
            (bool(*)()) &PetscTools::IsSequential,
            " "  )
        .def_static(
            "IsParallel",
            (bool(*)()) &PetscTools::IsParallel,
            " "  )
        .def_static(
            "IsIsolated",
            (bool(*)()) &PetscTools::IsIsolated,
            " "  )
        .def_static(
            "GetNumProcs",
            (unsigned int(*)()) &PetscTools::GetNumProcs,
            " "  )
        .def_static(
            "GetMyRank",
            (unsigned int(*)()) &PetscTools::GetMyRank,
            " "  )
        .def_static(
            "AmMaster",
            (bool(*)()) &PetscTools::AmMaster,
            " "  )
        .def_static(
            "AmTopMost",
            (bool(*)()) &PetscTools::AmTopMost,
            " "  )
        .def_static(
            "Barrier",
            (void(*)(::std::string const)) &PetscTools::Barrier,
            " " , py::arg("callerId") = "" )
        .def_static(
            "BeginRoundRobin",
            (void(*)()) &PetscTools::BeginRoundRobin,
            " "  )
        .def_static(
            "EndRoundRobin",
            (void(*)()) &PetscTools::EndRoundRobin,
            " "  )
        .def_static(
            "IsolateProcesses",
            (void(*)(bool)) &PetscTools::IsolateProcesses,
            " " , py::arg("isolate") = true )
        .def_static(
            "CreateVec",
            (::Vec(*)(int, int, bool)) &PetscTools::CreateVec,
            " " , py::arg("size"), py::arg("localSize") = -1, py::arg("ignoreOffProcEntries") = true , py::return_value_policy::reference)
        .def_static(
            "CreateVec",
            (::Vec(*)(::std::vector<double>)) &PetscTools::CreateVec,
            " " , py::arg("data") , py::return_value_policy::reference)
        .def_static(
            "CreateAndSetVec",
            (::Vec(*)(int, double)) &PetscTools::CreateAndSetVec,
            " " , py::arg("size"), py::arg("value") , py::return_value_policy::reference)
        .def_static(
            "ReplicateBool",
            (bool(*)(bool)) &PetscTools::ReplicateBool,
            " " , py::arg("flag") )
        .def_static(
            "ReplicateException",
            (void(*)(bool)) &PetscTools::ReplicateException,
            " " , py::arg("flag") )
        .def_static(
            "DumpPetscObject",
            (void(*)(::Mat const &, ::std::string const &)) &PetscTools::DumpPetscObject,
            " " , py::arg("rMat"), py::arg("rOutputFileFullPath") )
        .def_static(
            "DumpPetscObject",
            (void(*)(::Vec const &, ::std::string const &)) &PetscTools::DumpPetscObject,
            " " , py::arg("rVec"), py::arg("rOutputFileFullPath") )
        .def_static(
            "HasParMetis",
            (bool(*)()) &PetscTools::HasParMetis,
            " "  )
        .def_static(
            "SetOption",
            (void(*)(char const *, char const *)) &PetscTools::SetOption,
            " " , py::arg("pOptionName"), py::arg("pOptionValue") )
        .def_static(
            "ChasteMatCopy",
            (::PetscErrorCode(*)(::Mat, ::Mat, ::MatStructure)) &PetscTools::ChasteMatCopy,
            " " , py::arg("A"), py::arg("B"), py::arg("str") )
    ;
}
