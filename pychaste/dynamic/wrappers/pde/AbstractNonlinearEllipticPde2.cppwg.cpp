#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractNonlinearEllipticPde.hpp"

#include "AbstractNonlinearEllipticPde2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractNonlinearEllipticPde<2 > AbstractNonlinearEllipticPde2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class AbstractNonlinearEllipticPde2_Overrides : public AbstractNonlinearEllipticPde2{
    public:
    using AbstractNonlinearEllipticPde2::AbstractNonlinearEllipticPde;
    double ComputeLinearSourceTerm(::ChastePoint<2> const & rX) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractNonlinearEllipticPde2,
            ComputeLinearSourceTerm,
                    rX);
    }
    double ComputeNonlinearSourceTerm(::ChastePoint<2> const & rX, double u) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractNonlinearEllipticPde2,
            ComputeNonlinearSourceTerm,
                    rX,
        u);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX, double u) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            AbstractNonlinearEllipticPde2,
            ComputeDiffusionTerm,
                    rX,
        u);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTermPrime(::ChastePoint<2> const & rX, double u) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            AbstractNonlinearEllipticPde2,
            ComputeDiffusionTermPrime,
                    rX,
        u);
    }
    double ComputeNonlinearSourceTermPrime(::ChastePoint<2> const & rX, double u) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractNonlinearEllipticPde2,
            ComputeNonlinearSourceTermPrime,
                    rX,
        u);
    }

};
void register_AbstractNonlinearEllipticPde2_class(py::module &m){
py::class_<AbstractNonlinearEllipticPde2 , AbstractNonlinearEllipticPde2_Overrides , boost::shared_ptr<AbstractNonlinearEllipticPde2 >   >(m, "AbstractNonlinearEllipticPde2")
        .def(py::init< >())
        .def(
            "ComputeLinearSourceTerm",
            (double(AbstractNonlinearEllipticPde2::*)(::ChastePoint<2> const &)) &AbstractNonlinearEllipticPde2::ComputeLinearSourceTerm,
            " " , py::arg("rX") )
        .def(
            "ComputeNonlinearSourceTerm",
            (double(AbstractNonlinearEllipticPde2::*)(::ChastePoint<2> const &, double)) &AbstractNonlinearEllipticPde2::ComputeNonlinearSourceTerm,
            " " , py::arg("rX"), py::arg("u") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(AbstractNonlinearEllipticPde2::*)(::ChastePoint<2> const &, double)) &AbstractNonlinearEllipticPde2::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("u") )
        .def(
            "ComputeDiffusionTermPrime",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(AbstractNonlinearEllipticPde2::*)(::ChastePoint<2> const &, double)) &AbstractNonlinearEllipticPde2::ComputeDiffusionTermPrime,
            " " , py::arg("rX"), py::arg("u") )
        .def(
            "ComputeNonlinearSourceTermPrime",
            (double(AbstractNonlinearEllipticPde2::*)(::ChastePoint<2> const &, double)) &AbstractNonlinearEllipticPde2::ComputeNonlinearSourceTermPrime,
            " " , py::arg("rX"), py::arg("u") )
    ;
}
