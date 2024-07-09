#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "EdgeOperation.hpp"

#include "EdgeOperation.cppwg.hpp"

namespace py = pybind11;
typedef EdgeOperation EdgeOperation;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_EdgeOperation_class(py::module &m){
py::class_<EdgeOperation  , boost::shared_ptr<EdgeOperation >   >(m, "EdgeOperation")
        .def(py::init< >())
        .def(py::init<::EDGE_OPERATION, unsigned int, ::EdgeRemapInfo, bool const >(), py::arg("operation"), py::arg("elementIndex"), py::arg("remapInfo"), py::arg("isIndexRemapped") = false)
        .def(py::init<unsigned int, unsigned int, ::EdgeRemapInfo, ::EdgeRemapInfo >(), py::arg("elementIndex"), py::arg("elementIndex2"), py::arg("remapInfo"), py::arg("remapInfo2"))
        .def(
            "GetOperation",
            (::EDGE_OPERATION(EdgeOperation::*)() const ) &EdgeOperation::GetOperation,
            " "  )
        .def(
            "GetElementIndex",
            (unsigned int(EdgeOperation::*)() const ) &EdgeOperation::GetElementIndex,
            " "  )
        .def(
            "SetElementIndex",
            (void(EdgeOperation::*)(unsigned int const)) &EdgeOperation::SetElementIndex,
            " " , py::arg("index") )
        .def(
            "GetElementIndex2",
            (unsigned int(EdgeOperation::*)() const ) &EdgeOperation::GetElementIndex2,
            " "  )
        .def(
            "SetElementIndex2",
            (void(EdgeOperation::*)(unsigned int const)) &EdgeOperation::SetElementIndex2,
            " " , py::arg("index") )
        .def(
            "rGetRemapInfo",
            (::EdgeRemapInfo const &(EdgeOperation::*)() const ) &EdgeOperation::rGetRemapInfo,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetRemapInfo2",
            (::EdgeRemapInfo const &(EdgeOperation::*)() const ) &EdgeOperation::rGetRemapInfo2,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "IsElementIndexRemapped",
            (bool(EdgeOperation::*)() const ) &EdgeOperation::IsElementIndexRemapped,
            " "  )
    ;
}
