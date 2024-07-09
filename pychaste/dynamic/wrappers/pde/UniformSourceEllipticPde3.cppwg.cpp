#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "UniformSourceEllipticPde.hpp"

#include "UniformSourceEllipticPde3.cppwg.hpp"

namespace py = pybind11;
typedef UniformSourceEllipticPde<3 > UniformSourceEllipticPde3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class UniformSourceEllipticPde3_Overrides : public UniformSourceEllipticPde3{
    public:
    using UniformSourceEllipticPde3::UniformSourceEllipticPde;
    double ComputeConstantInUSourceTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            UniformSourceEllipticPde3,
            ComputeConstantInUSourceTerm,
                    rX,
        pElement);
    }
    double ComputeLinearInUCoeffInSourceTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            UniformSourceEllipticPde3,
            ComputeLinearInUCoeffInSourceTerm,
                    rX,
        pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTerm(::ChastePoint<3> const & rX) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            UniformSourceEllipticPde3,
            ComputeDiffusionTerm,
                    rX);
    }

};
void register_UniformSourceEllipticPde3_class(py::module &m){
py::class_<UniformSourceEllipticPde3 , UniformSourceEllipticPde3_Overrides , boost::shared_ptr<UniformSourceEllipticPde3 >  , AbstractLinearEllipticPde<3, 3>  >(m, "UniformSourceEllipticPde3")
        .def(py::init<double >(), py::arg("sourceCoefficient") = 0.)
        .def(
            "GetCoefficient",
            (double(UniformSourceEllipticPde3::*)() const ) &UniformSourceEllipticPde3::GetCoefficient,
            " "  )
        .def(
            "ComputeConstantInUSourceTerm",
            (double(UniformSourceEllipticPde3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &UniformSourceEllipticPde3::ComputeConstantInUSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeLinearInUCoeffInSourceTerm",
            (double(UniformSourceEllipticPde3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &UniformSourceEllipticPde3::ComputeLinearInUCoeffInSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(UniformSourceEllipticPde3::*)(::ChastePoint<3> const &)) &UniformSourceEllipticPde3::ComputeDiffusionTerm,
            " " , py::arg("rX") )
    ;
}
