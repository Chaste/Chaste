#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractLinearParabolicPde.hpp"

#include "AbstractLinearParabolicPde3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearParabolicPde<3,3 > AbstractLinearParabolicPde3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class AbstractLinearParabolicPde3_3_Overrides : public AbstractLinearParabolicPde3_3{
    public:
    using AbstractLinearParabolicPde3_3::AbstractLinearParabolicPde;
    double ComputeDuDtCoefficientFunction(::ChastePoint<3> const & rX) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPde3_3,
            ComputeDuDtCoefficientFunction,
                    rX);
    }
    double ComputeSourceTerm(::ChastePoint<3> const & rX, double u, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPde3_3,
            ComputeSourceTerm,
                    rX,
        u,
        pElement);
    }
    double ComputeSourceTermAtNode(::Node<3> const & rNode, double u) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearParabolicPde3_3,
            ComputeSourceTermAtNode,
                    rNode,
        u);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            AbstractLinearParabolicPde3_3,
            ComputeDiffusionTerm,
                    rX,
        pElement);
    }

};
void register_AbstractLinearParabolicPde3_3_class(py::module &m){
py::class_<AbstractLinearParabolicPde3_3 , AbstractLinearParabolicPde3_3_Overrides , boost::shared_ptr<AbstractLinearParabolicPde3_3 >  , AbstractLinearPde<3>  >(m, "AbstractLinearParabolicPde3_3")
        .def(py::init< >())
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(AbstractLinearParabolicPde3_3::*)(::ChastePoint<3> const &)) &AbstractLinearParabolicPde3_3::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX") )
        .def(
            "ComputeSourceTerm",
            (double(AbstractLinearParabolicPde3_3::*)(::ChastePoint<3> const &, double, ::Element<3, 3> *)) &AbstractLinearParabolicPde3_3::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("u"), py::arg("pElement") = nullptr )
        .def(
            "ComputeSourceTermAtNode",
            (double(AbstractLinearParabolicPde3_3::*)(::Node<3> const &, double)) &AbstractLinearParabolicPde3_3::ComputeSourceTermAtNode,
            " " , py::arg("rNode"), py::arg("u") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(AbstractLinearParabolicPde3_3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &AbstractLinearParabolicPde3_3::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pElement") = __null )
    ;
}
