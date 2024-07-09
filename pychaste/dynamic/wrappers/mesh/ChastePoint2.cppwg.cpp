#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ChastePoint.hpp"

#include "ChastePoint2.cppwg.hpp"

namespace py = pybind11;
typedef ChastePoint<2 > ChastePoint2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_ChastePoint2_class(py::module &m){
py::class_<ChastePoint2  , boost::shared_ptr<ChastePoint2 >   >(m, "ChastePoint2")
        .def(py::init<double, double, double >(), py::arg("v1") = 0, py::arg("v2") = 0, py::arg("v3") = 0)
        .def(py::init<::std::vector<double> >(), py::arg("coords"))
        .def(py::init<::boost::numeric::ublas::c_vector<double, 2> >(), py::arg("location"))
        .def(
            "rGetLocation",
            (::boost::numeric::ublas::c_vector<double, 2> &(ChastePoint2::*)()) &ChastePoint2::rGetLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetLocation",
            (::boost::numeric::ublas::c_vector<double, 2> const &(ChastePoint2::*)() const ) &ChastePoint2::rGetLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetWithDefault",
            (double(ChastePoint2::*)(unsigned int, double) const ) &ChastePoint2::GetWithDefault,
            " " , py::arg("i"), py::arg("def") = 0. )
        .def(
            "SetCoordinate",
            (void(ChastePoint2::*)(unsigned int, double)) &ChastePoint2::SetCoordinate,
            " " , py::arg("i"), py::arg("value") )
        .def(
            "IsSamePoint",
            (bool(ChastePoint2::*)(::ChastePoint<2> const &) const ) &ChastePoint2::IsSamePoint,
            " " , py::arg("rPoint") )
    ;
}
