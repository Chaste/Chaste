#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PottsElement.hpp"

#include "PottsElement3.cppwg.hpp"

namespace py = pybind11;
typedef PottsElement<3 > PottsElement3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_PottsElement3_class(py::module &m){
py::class_<PottsElement3  , boost::shared_ptr<PottsElement3 >  , MutableElement<3, 3>  >(m, "PottsElement3")
        .def(py::init<unsigned int, ::std::vector<Node<3> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(
            "AddNode",
            (void(PottsElement3::*)(::Node<3> *, unsigned int const &)) &PottsElement3::AddNode,
            " " , py::arg("pNode"), py::arg("rIndex") = (2147483647 * 2U + 1U) )
        .def(
            "GetAspectRatio",
            (double(PottsElement3::*)()) &PottsElement3::GetAspectRatio,
            " "  )
    ;
}
