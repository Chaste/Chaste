#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ChastePoint.hpp"

#include "ChastePoint3.cppwg.hpp"

namespace py = pybind11;
typedef ChastePoint<3 > ChastePoint3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_ChastePoint3_class(py::module &m){
py::class_<ChastePoint3  , boost::shared_ptr<ChastePoint3 >   >(m, "ChastePoint3")
        .def(py::init<double, double, double >(), py::arg("v1") = 0, py::arg("v2") = 0, py::arg("v3") = 0)
        .def(py::init<::std::vector<double> >(), py::arg("coords"))
        .def(py::init<::boost::numeric::ublas::c_vector<double, 3> >(), py::arg("location"))
        .def(
            "rGetLocation",
            (::boost::numeric::ublas::c_vector<double, 3> &(ChastePoint3::*)()) &ChastePoint3::rGetLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetLocation",
            (::boost::numeric::ublas::c_vector<double, 3> const &(ChastePoint3::*)() const ) &ChastePoint3::rGetLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetWithDefault",
            (double(ChastePoint3::*)(unsigned int, double) const ) &ChastePoint3::GetWithDefault,
            " " , py::arg("i"), py::arg("def") = 0. )
        .def(
            "SetCoordinate",
            (void(ChastePoint3::*)(unsigned int, double)) &ChastePoint3::SetCoordinate,
            " " , py::arg("i"), py::arg("value") )
        .def(
            "IsSamePoint",
            (bool(ChastePoint3::*)(::ChastePoint<3> const &) const ) &ChastePoint3::IsSamePoint,
            " " , py::arg("rPoint") )
    ;
}
