#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "UniformSourceParabolicPde.hpp"

#include "UniformSourceParabolicPde3.cppwg.hpp"

namespace py = pybind11;
typedef UniformSourceParabolicPde<3 > UniformSourceParabolicPde3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class UniformSourceParabolicPde3_Overrides : public UniformSourceParabolicPde3{
    public:
    using UniformSourceParabolicPde3::UniformSourceParabolicPde;
    double ComputeSourceTerm(::ChastePoint<3> const & rX, double u, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            UniformSourceParabolicPde3,
            ComputeSourceTerm,
                    rX,
        u,
        pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            UniformSourceParabolicPde3,
            ComputeDiffusionTerm,
                    rX,
        pElement);
    }
    double ComputeDuDtCoefficientFunction(::ChastePoint<3> const & rX) override {
        PYBIND11_OVERRIDE(
            double,
            UniformSourceParabolicPde3,
            ComputeDuDtCoefficientFunction,
                    rX);
    }

};
void register_UniformSourceParabolicPde3_class(py::module &m){
py::class_<UniformSourceParabolicPde3 , UniformSourceParabolicPde3_Overrides , boost::shared_ptr<UniformSourceParabolicPde3 >  , AbstractLinearParabolicPde<3>  >(m, "UniformSourceParabolicPde3")
        .def(py::init<double >(), py::arg("sourceCoefficient") = 0.)
        .def(
            "GetCoefficient",
            (double(UniformSourceParabolicPde3::*)() const ) &UniformSourceParabolicPde3::GetCoefficient,
            " "  )
        .def(
            "ComputeSourceTerm",
            (double(UniformSourceParabolicPde3::*)(::ChastePoint<3> const &, double, ::Element<3, 3> *)) &UniformSourceParabolicPde3::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("u"), py::arg("pElement") = __null )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(UniformSourceParabolicPde3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &UniformSourceParabolicPde3::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pElement") = __null )
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(UniformSourceParabolicPde3::*)(::ChastePoint<3> const &)) &UniformSourceParabolicPde3::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX") )
    ;
}
