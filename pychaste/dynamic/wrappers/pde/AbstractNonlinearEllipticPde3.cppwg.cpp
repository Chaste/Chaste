#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractNonlinearEllipticPde.hpp"

#include "AbstractNonlinearEllipticPde3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractNonlinearEllipticPde<3 > AbstractNonlinearEllipticPde3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class AbstractNonlinearEllipticPde3_Overrides : public AbstractNonlinearEllipticPde3{
    public:
    using AbstractNonlinearEllipticPde3::AbstractNonlinearEllipticPde;
    double ComputeLinearSourceTerm(::ChastePoint<3> const & rX) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractNonlinearEllipticPde3,
            ComputeLinearSourceTerm,
                    rX);
    }
    double ComputeNonlinearSourceTerm(::ChastePoint<3> const & rX, double u) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractNonlinearEllipticPde3,
            ComputeNonlinearSourceTerm,
                    rX,
        u);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTerm(::ChastePoint<3> const & rX, double u) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            AbstractNonlinearEllipticPde3,
            ComputeDiffusionTerm,
                    rX,
        u);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTermPrime(::ChastePoint<3> const & rX, double u) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            AbstractNonlinearEllipticPde3,
            ComputeDiffusionTermPrime,
                    rX,
        u);
    }
    double ComputeNonlinearSourceTermPrime(::ChastePoint<3> const & rX, double u) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractNonlinearEllipticPde3,
            ComputeNonlinearSourceTermPrime,
                    rX,
        u);
    }

};
void register_AbstractNonlinearEllipticPde3_class(py::module &m){
py::class_<AbstractNonlinearEllipticPde3 , AbstractNonlinearEllipticPde3_Overrides , boost::shared_ptr<AbstractNonlinearEllipticPde3 >   >(m, "AbstractNonlinearEllipticPde3")
        .def(py::init< >())
        .def(
            "ComputeLinearSourceTerm",
            (double(AbstractNonlinearEllipticPde3::*)(::ChastePoint<3> const &)) &AbstractNonlinearEllipticPde3::ComputeLinearSourceTerm,
            " " , py::arg("rX") )
        .def(
            "ComputeNonlinearSourceTerm",
            (double(AbstractNonlinearEllipticPde3::*)(::ChastePoint<3> const &, double)) &AbstractNonlinearEllipticPde3::ComputeNonlinearSourceTerm,
            " " , py::arg("rX"), py::arg("u") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(AbstractNonlinearEllipticPde3::*)(::ChastePoint<3> const &, double)) &AbstractNonlinearEllipticPde3::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("u") )
        .def(
            "ComputeDiffusionTermPrime",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(AbstractNonlinearEllipticPde3::*)(::ChastePoint<3> const &, double)) &AbstractNonlinearEllipticPde3::ComputeDiffusionTermPrime,
            " " , py::arg("rX"), py::arg("u") )
        .def(
            "ComputeNonlinearSourceTermPrime",
            (double(AbstractNonlinearEllipticPde3::*)(::ChastePoint<3> const &, double)) &AbstractNonlinearEllipticPde3::ComputeNonlinearSourceTermPrime,
            " " , py::arg("rX"), py::arg("u") )
    ;
}
