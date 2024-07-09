#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearParabolicPdeSystemForCoupledOdeSystem<2,2,1 > AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1_Overrides : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1{
    public:
    using AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::AbstractLinearParabolicPdeSystemForCoupledOdeSystem;
    double ComputeDuDtCoefficientFunction(::ChastePoint<2> const & rX, unsigned int pdeIndex) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1,
            ComputeDuDtCoefficientFunction,
                    rX,
        pdeIndex);
    }
    double ComputeSourceTerm(::ChastePoint<2> const & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::std::vector<double> & rOdeSolution, unsigned int pdeIndex) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1,
            ComputeSourceTerm,
                    rX,
        rU,
        rOdeSolution,
        pdeIndex);
    }
    double ComputeSourceTermAtNode(::Node<2> const & rNode, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::std::vector<double> & rOdeSolution, unsigned int pdeIndex) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1,
            ComputeSourceTermAtNode,
                    rNode,
        rU,
        rOdeSolution,
        pdeIndex);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX, unsigned int pdeIndex, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1,
            ComputeDiffusionTerm,
                    rX,
        pdeIndex,
        pElement);
    }

};
void register_AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1_class(py::module &m){
py::class_<AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1 , AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1_Overrides , boost::shared_ptr<AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1 >   >(m, "AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1")
        .def(py::init< >())
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::*)(::ChastePoint<2> const &, unsigned int)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX"), py::arg("pdeIndex") )
        .def(
            "ComputeSourceTerm",
            (double(AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::*)(::ChastePoint<2> const &, ::boost::numeric::ublas::c_vector<double, 1> &, ::std::vector<double> &, unsigned int)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("rU"), py::arg("rOdeSolution"), py::arg("pdeIndex") )
        .def(
            "ComputeSourceTermAtNode",
            (double(AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::*)(::Node<2> const &, ::boost::numeric::ublas::c_vector<double, 1> &, ::std::vector<double> &, unsigned int)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::ComputeSourceTermAtNode,
            " " , py::arg("rNode"), py::arg("rU"), py::arg("rOdeSolution"), py::arg("pdeIndex") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::*)(::ChastePoint<2> const &, unsigned int, ::Element<2, 2> *)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pdeIndex"), py::arg("pElement") = __null )
    ;
}
