#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractLinearEllipticPde.hpp"

#include "AbstractLinearEllipticPde2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearEllipticPde<2,2 > AbstractLinearEllipticPde2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class AbstractLinearEllipticPde2_2_Overrides : public AbstractLinearEllipticPde2_2{
    public:
    using AbstractLinearEllipticPde2_2::AbstractLinearEllipticPde;
    double ComputeConstantInUSourceTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearEllipticPde2_2,
            ComputeConstantInUSourceTerm,
                    rX,
        pElement);
    }
    double ComputeLinearInUCoeffInSourceTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearEllipticPde2_2,
            ComputeLinearInUCoeffInSourceTerm,
                    rX,
        pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            AbstractLinearEllipticPde2_2,
            ComputeDiffusionTerm,
                    rX);
    }
    double ComputeConstantInUSourceTermAtNode(::Node<2> const & rNode) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearEllipticPde2_2,
            ComputeConstantInUSourceTermAtNode,
                    rNode);
    }
    double ComputeLinearInUCoeffInSourceTermAtNode(::Node<2> const & rNode) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearEllipticPde2_2,
            ComputeLinearInUCoeffInSourceTermAtNode,
                    rNode);
    }

};
void register_AbstractLinearEllipticPde2_2_class(py::module &m){
py::class_<AbstractLinearEllipticPde2_2 , AbstractLinearEllipticPde2_2_Overrides , boost::shared_ptr<AbstractLinearEllipticPde2_2 >  , AbstractLinearPde<2>  >(m, "AbstractLinearEllipticPde2_2")
        .def(py::init< >())
        .def(
            "ComputeConstantInUSourceTerm",
            (double(AbstractLinearEllipticPde2_2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &AbstractLinearEllipticPde2_2::ComputeConstantInUSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeLinearInUCoeffInSourceTerm",
            (double(AbstractLinearEllipticPde2_2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &AbstractLinearEllipticPde2_2::ComputeLinearInUCoeffInSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(AbstractLinearEllipticPde2_2::*)(::ChastePoint<2> const &)) &AbstractLinearEllipticPde2_2::ComputeDiffusionTerm,
            " " , py::arg("rX") )
        .def(
            "ComputeConstantInUSourceTermAtNode",
            (double(AbstractLinearEllipticPde2_2::*)(::Node<2> const &)) &AbstractLinearEllipticPde2_2::ComputeConstantInUSourceTermAtNode,
            " " , py::arg("rNode") )
        .def(
            "ComputeLinearInUCoeffInSourceTermAtNode",
            (double(AbstractLinearEllipticPde2_2::*)(::Node<2> const &)) &AbstractLinearEllipticPde2_2::ComputeLinearInUCoeffInSourceTermAtNode,
            " " , py::arg("rNode") )
    ;
}
