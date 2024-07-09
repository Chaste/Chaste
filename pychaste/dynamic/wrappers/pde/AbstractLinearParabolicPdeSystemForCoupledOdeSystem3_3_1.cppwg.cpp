#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearParabolicPdeSystemForCoupledOdeSystem<3,3,1 > AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1_Overrides : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1{
    public:
    using AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1::AbstractLinearParabolicPdeSystemForCoupledOdeSystem;
    double ComputeDuDtCoefficientFunction(::ChastePoint<3> const & rX, unsigned int pdeIndex) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1,
            ComputeDuDtCoefficientFunction,
                    rX,
        pdeIndex);
    }
    double ComputeSourceTerm(::ChastePoint<3> const & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::std::vector<double> & rOdeSolution, unsigned int pdeIndex) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1,
            ComputeSourceTerm,
                    rX,
        rU,
        rOdeSolution,
        pdeIndex);
    }
    double ComputeSourceTermAtNode(::Node<3> const & rNode, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::std::vector<double> & rOdeSolution, unsigned int pdeIndex) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1,
            ComputeSourceTermAtNode,
                    rNode,
        rU,
        rOdeSolution,
        pdeIndex);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTerm(::ChastePoint<3> const & rX, unsigned int pdeIndex, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1,
            ComputeDiffusionTerm,
                    rX,
        pdeIndex,
        pElement);
    }

};
void register_AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1_class(py::module &m){
py::class_<AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1 , AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1_Overrides , boost::shared_ptr<AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1 >   >(m, "AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1")
        .def(py::init< >())
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1::*)(::ChastePoint<3> const &, unsigned int)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX"), py::arg("pdeIndex") )
        .def(
            "ComputeSourceTerm",
            (double(AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1::*)(::ChastePoint<3> const &, ::boost::numeric::ublas::c_vector<double, 1> &, ::std::vector<double> &, unsigned int)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("rU"), py::arg("rOdeSolution"), py::arg("pdeIndex") )
        .def(
            "ComputeSourceTermAtNode",
            (double(AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1::*)(::Node<3> const &, ::boost::numeric::ublas::c_vector<double, 1> &, ::std::vector<double> &, unsigned int)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1::ComputeSourceTermAtNode,
            " " , py::arg("rNode"), py::arg("rU"), py::arg("rOdeSolution"), py::arg("pdeIndex") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1::*)(::ChastePoint<3> const &, unsigned int, ::Element<3, 3> *)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem3_3_1::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pdeIndex"), py::arg("pElement") = __null )
    ;
}
