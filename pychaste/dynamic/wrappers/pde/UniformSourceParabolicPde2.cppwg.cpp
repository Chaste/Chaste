#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "UniformSourceParabolicPde.hpp"

#include "UniformSourceParabolicPde2.cppwg.hpp"

namespace py = pybind11;
typedef UniformSourceParabolicPde<2 > UniformSourceParabolicPde2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class UniformSourceParabolicPde2_Overrides : public UniformSourceParabolicPde2{
    public:
    using UniformSourceParabolicPde2::UniformSourceParabolicPde;
    double ComputeSourceTerm(::ChastePoint<2> const & rX, double u, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            UniformSourceParabolicPde2,
            ComputeSourceTerm,
                    rX,
        u,
        pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            UniformSourceParabolicPde2,
            ComputeDiffusionTerm,
                    rX,
        pElement);
    }
    double ComputeDuDtCoefficientFunction(::ChastePoint<2> const & rX) override {
        PYBIND11_OVERRIDE(
            double,
            UniformSourceParabolicPde2,
            ComputeDuDtCoefficientFunction,
                    rX);
    }

};
void register_UniformSourceParabolicPde2_class(py::module &m){
py::class_<UniformSourceParabolicPde2 , UniformSourceParabolicPde2_Overrides , boost::shared_ptr<UniformSourceParabolicPde2 >  , AbstractLinearParabolicPde<2>  >(m, "UniformSourceParabolicPde2")
        .def(py::init<double >(), py::arg("sourceCoefficient") = 0.)
        .def(
            "GetCoefficient",
            (double(UniformSourceParabolicPde2::*)() const ) &UniformSourceParabolicPde2::GetCoefficient,
            " "  )
        .def(
            "ComputeSourceTerm",
            (double(UniformSourceParabolicPde2::*)(::ChastePoint<2> const &, double, ::Element<2, 2> *)) &UniformSourceParabolicPde2::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("u"), py::arg("pElement") = __null )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(UniformSourceParabolicPde2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &UniformSourceParabolicPde2::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pElement") = __null )
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(UniformSourceParabolicPde2::*)(::ChastePoint<2> const &)) &UniformSourceParabolicPde2::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX") )
    ;
}
