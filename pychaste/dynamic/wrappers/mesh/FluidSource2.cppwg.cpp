#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FluidSource.hpp"

#include "FluidSource2.cppwg.hpp"

namespace py = pybind11;
typedef FluidSource<2 > FluidSource2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_FluidSource2_class(py::module &m){
py::class_<FluidSource2  , boost::shared_ptr<FluidSource2 >   >(m, "FluidSource2")
        .def(py::init<unsigned int, ::ChastePoint<2> >(), py::arg("index"), py::arg("point"))
        .def(py::init<unsigned int, ::boost::numeric::ublas::c_vector<double, 2> >(), py::arg("index"), py::arg("location"))
        .def(py::init<unsigned int, double, double, double >(), py::arg("index"), py::arg("v1") = 0., py::arg("v2") = 0., py::arg("v3") = 0.)
        .def(
            "GetIndex",
            (unsigned int(FluidSource2::*)() const ) &FluidSource2::GetIndex,
            " "  )
        .def(
            "SetIndex",
            (void(FluidSource2::*)(unsigned int)) &FluidSource2::SetIndex,
            " " , py::arg("index") )
        .def(
            "GetPoint",
            (::ChastePoint<2>(FluidSource2::*)() const ) &FluidSource2::GetPoint,
            " "  )
        .def(
            "rGetLocation",
            (::boost::numeric::ublas::c_vector<double, 2> const &(FluidSource2::*)() const ) &FluidSource2::rGetLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetModifiableLocation",
            (::boost::numeric::ublas::c_vector<double, 2> &(FluidSource2::*)()) &FluidSource2::rGetModifiableLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetStrength",
            (double(FluidSource2::*)() const ) &FluidSource2::GetStrength,
            " "  )
        .def(
            "SetStrength",
            (void(FluidSource2::*)(double)) &FluidSource2::SetStrength,
            " " , py::arg("strength") )
        .def(
            "SetIfSourceIsAssociatedWithElement",
            (void(FluidSource2::*)(bool)) &FluidSource2::SetIfSourceIsAssociatedWithElement,
            " " , py::arg("associated") )
        .def(
            "IsSourceAssociatedWithElement",
            (bool(FluidSource2::*)()) &FluidSource2::IsSourceAssociatedWithElement,
            " "  )
        .def(
            "GetAssociatedElementIndex",
            (unsigned int(FluidSource2::*)() const ) &FluidSource2::GetAssociatedElementIndex,
            " "  )
        .def(
            "SetAssociatedElementIndex",
            (void(FluidSource2::*)(unsigned int)) &FluidSource2::SetAssociatedElementIndex,
            " " , py::arg("associatedElementIndex") )
    ;
}
