#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "UniformSourceEllipticPde.hpp"

#include "UniformSourceEllipticPde2.cppwg.hpp"

namespace py = pybind11;
typedef UniformSourceEllipticPde<2 > UniformSourceEllipticPde2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class UniformSourceEllipticPde2_Overrides : public UniformSourceEllipticPde2{
    public:
    using UniformSourceEllipticPde2::UniformSourceEllipticPde;
    double ComputeConstantInUSourceTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            UniformSourceEllipticPde2,
            ComputeConstantInUSourceTerm,
                    rX,
        pElement);
    }
    double ComputeLinearInUCoeffInSourceTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            UniformSourceEllipticPde2,
            ComputeLinearInUCoeffInSourceTerm,
                    rX,
        pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            UniformSourceEllipticPde2,
            ComputeDiffusionTerm,
                    rX);
    }

};
void register_UniformSourceEllipticPde2_class(py::module &m){
py::class_<UniformSourceEllipticPde2 , UniformSourceEllipticPde2_Overrides , boost::shared_ptr<UniformSourceEllipticPde2 >  , AbstractLinearEllipticPde<2, 2>  >(m, "UniformSourceEllipticPde2")
        .def(py::init<double >(), py::arg("sourceCoefficient") = 0.)
        .def(
            "GetCoefficient",
            (double(UniformSourceEllipticPde2::*)() const ) &UniformSourceEllipticPde2::GetCoefficient,
            " "  )
        .def(
            "ComputeConstantInUSourceTerm",
            (double(UniformSourceEllipticPde2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &UniformSourceEllipticPde2::ComputeConstantInUSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeLinearInUCoeffInSourceTerm",
            (double(UniformSourceEllipticPde2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &UniformSourceEllipticPde2::ComputeLinearInUCoeffInSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(UniformSourceEllipticPde2::*)(::ChastePoint<2> const &)) &UniformSourceEllipticPde2::ComputeDiffusionTerm,
            " " , py::arg("rX") )
    ;
}
