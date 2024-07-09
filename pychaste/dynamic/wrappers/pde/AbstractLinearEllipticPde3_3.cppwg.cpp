#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractLinearEllipticPde.hpp"

#include "AbstractLinearEllipticPde3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearEllipticPde<3,3 > AbstractLinearEllipticPde3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class AbstractLinearEllipticPde3_3_Overrides : public AbstractLinearEllipticPde3_3{
    public:
    using AbstractLinearEllipticPde3_3::AbstractLinearEllipticPde;
    double ComputeConstantInUSourceTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearEllipticPde3_3,
            ComputeConstantInUSourceTerm,
                    rX,
        pElement);
    }
    double ComputeLinearInUCoeffInSourceTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearEllipticPde3_3,
            ComputeLinearInUCoeffInSourceTerm,
                    rX,
        pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTerm(::ChastePoint<3> const & rX) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            AbstractLinearEllipticPde3_3,
            ComputeDiffusionTerm,
                    rX);
    }
    double ComputeConstantInUSourceTermAtNode(::Node<3> const & rNode) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearEllipticPde3_3,
            ComputeConstantInUSourceTermAtNode,
                    rNode);
    }
    double ComputeLinearInUCoeffInSourceTermAtNode(::Node<3> const & rNode) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearEllipticPde3_3,
            ComputeLinearInUCoeffInSourceTermAtNode,
                    rNode);
    }

};
void register_AbstractLinearEllipticPde3_3_class(py::module &m){
py::class_<AbstractLinearEllipticPde3_3 , AbstractLinearEllipticPde3_3_Overrides , boost::shared_ptr<AbstractLinearEllipticPde3_3 >  , AbstractLinearPde<3>  >(m, "AbstractLinearEllipticPde3_3")
        .def(py::init< >())
        .def(
            "ComputeConstantInUSourceTerm",
            (double(AbstractLinearEllipticPde3_3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &AbstractLinearEllipticPde3_3::ComputeConstantInUSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeLinearInUCoeffInSourceTerm",
            (double(AbstractLinearEllipticPde3_3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &AbstractLinearEllipticPde3_3::ComputeLinearInUCoeffInSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(AbstractLinearEllipticPde3_3::*)(::ChastePoint<3> const &)) &AbstractLinearEllipticPde3_3::ComputeDiffusionTerm,
            " " , py::arg("rX") )
        .def(
            "ComputeConstantInUSourceTermAtNode",
            (double(AbstractLinearEllipticPde3_3::*)(::Node<3> const &)) &AbstractLinearEllipticPde3_3::ComputeConstantInUSourceTermAtNode,
            " " , py::arg("rNode") )
        .def(
            "ComputeLinearInUCoeffInSourceTermAtNode",
            (double(AbstractLinearEllipticPde3_3::*)(::Node<3> const &)) &AbstractLinearEllipticPde3_3::ComputeLinearInUCoeffInSourceTermAtNode,
            " " , py::arg("rNode") )
    ;
}
