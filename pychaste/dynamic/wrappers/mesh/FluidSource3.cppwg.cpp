#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FluidSource.hpp"

#include "FluidSource3.cppwg.hpp"

namespace py = pybind11;
typedef FluidSource<3 > FluidSource3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_FluidSource3_class(py::module &m){
py::class_<FluidSource3  , boost::shared_ptr<FluidSource3 >   >(m, "FluidSource3")
        .def(py::init<unsigned int, ::ChastePoint<3> >(), py::arg("index"), py::arg("point"))
        .def(py::init<unsigned int, ::boost::numeric::ublas::c_vector<double, 3> >(), py::arg("index"), py::arg("location"))
        .def(py::init<unsigned int, double, double, double >(), py::arg("index"), py::arg("v1") = 0., py::arg("v2") = 0., py::arg("v3") = 0.)
        .def(
            "GetIndex",
            (unsigned int(FluidSource3::*)() const ) &FluidSource3::GetIndex,
            " "  )
        .def(
            "SetIndex",
            (void(FluidSource3::*)(unsigned int)) &FluidSource3::SetIndex,
            " " , py::arg("index") )
        .def(
            "GetPoint",
            (::ChastePoint<3>(FluidSource3::*)() const ) &FluidSource3::GetPoint,
            " "  )
        .def(
            "rGetLocation",
            (::boost::numeric::ublas::c_vector<double, 3> const &(FluidSource3::*)() const ) &FluidSource3::rGetLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetModifiableLocation",
            (::boost::numeric::ublas::c_vector<double, 3> &(FluidSource3::*)()) &FluidSource3::rGetModifiableLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetStrength",
            (double(FluidSource3::*)() const ) &FluidSource3::GetStrength,
            " "  )
        .def(
            "SetStrength",
            (void(FluidSource3::*)(double)) &FluidSource3::SetStrength,
            " " , py::arg("strength") )
        .def(
            "SetIfSourceIsAssociatedWithElement",
            (void(FluidSource3::*)(bool)) &FluidSource3::SetIfSourceIsAssociatedWithElement,
            " " , py::arg("associated") )
        .def(
            "IsSourceAssociatedWithElement",
            (bool(FluidSource3::*)()) &FluidSource3::IsSourceAssociatedWithElement,
            " "  )
        .def(
            "GetAssociatedElementIndex",
            (unsigned int(FluidSource3::*)() const ) &FluidSource3::GetAssociatedElementIndex,
            " "  )
        .def(
            "SetAssociatedElementIndex",
            (void(FluidSource3::*)(unsigned int)) &FluidSource3::SetAssociatedElementIndex,
            " " , py::arg("associatedElementIndex") )
    ;
}
