#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractLinearParabolicPde.hpp"

#include "AbstractLinearParabolicPde2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearParabolicPde<2,2 > AbstractLinearParabolicPde2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class AbstractLinearParabolicPde2_2_Overrides : public AbstractLinearParabolicPde2_2{
    public:
    using AbstractLinearParabolicPde2_2::AbstractLinearParabolicPde;
    double ComputeDuDtCoefficientFunction(::ChastePoint<2> const & rX) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPde2_2,
            ComputeDuDtCoefficientFunction,
                    rX);
    }
    double ComputeSourceTerm(::ChastePoint<2> const & rX, double u, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPde2_2,
            ComputeSourceTerm,
                    rX,
        u,
        pElement);
    }
    double ComputeSourceTermAtNode(::Node<2> const & rNode, double u) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearParabolicPde2_2,
            ComputeSourceTermAtNode,
                    rNode,
        u);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            AbstractLinearParabolicPde2_2,
            ComputeDiffusionTerm,
                    rX,
        pElement);
    }

};
void register_AbstractLinearParabolicPde2_2_class(py::module &m){
py::class_<AbstractLinearParabolicPde2_2 , AbstractLinearParabolicPde2_2_Overrides , boost::shared_ptr<AbstractLinearParabolicPde2_2 >  , AbstractLinearPde<2>  >(m, "AbstractLinearParabolicPde2_2")
        .def(py::init< >())
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(AbstractLinearParabolicPde2_2::*)(::ChastePoint<2> const &)) &AbstractLinearParabolicPde2_2::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX") )
        .def(
            "ComputeSourceTerm",
            (double(AbstractLinearParabolicPde2_2::*)(::ChastePoint<2> const &, double, ::Element<2, 2> *)) &AbstractLinearParabolicPde2_2::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("u"), py::arg("pElement") = nullptr )
        .def(
            "ComputeSourceTermAtNode",
            (double(AbstractLinearParabolicPde2_2::*)(::Node<2> const &, double)) &AbstractLinearParabolicPde2_2::ComputeSourceTermAtNode,
            " " , py::arg("rNode"), py::arg("u") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(AbstractLinearParabolicPde2_2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &AbstractLinearParabolicPde2_2::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pElement") = __null )
    ;
}
