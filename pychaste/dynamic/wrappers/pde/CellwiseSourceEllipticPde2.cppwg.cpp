#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellwiseSourceEllipticPde.hpp"

#include "CellwiseSourceEllipticPde2.cppwg.hpp"

namespace py = pybind11;
typedef CellwiseSourceEllipticPde<2 > CellwiseSourceEllipticPde2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class CellwiseSourceEllipticPde2_Overrides : public CellwiseSourceEllipticPde2{
    public:
    using CellwiseSourceEllipticPde2::CellwiseSourceEllipticPde;
    double ComputeConstantInUSourceTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            CellwiseSourceEllipticPde2,
            ComputeConstantInUSourceTerm,
                    rX,
        pElement);
    }
    double ComputeLinearInUCoeffInSourceTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            CellwiseSourceEllipticPde2,
            ComputeLinearInUCoeffInSourceTerm,
                    rX,
        pElement);
    }
    double ComputeLinearInUCoeffInSourceTermAtNode(::Node<2> const & rNode) override {
        PYBIND11_OVERRIDE(
            double,
            CellwiseSourceEllipticPde2,
            ComputeLinearInUCoeffInSourceTermAtNode,
                    rNode);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            CellwiseSourceEllipticPde2,
            ComputeDiffusionTerm,
                    rX);
    }

};
void register_CellwiseSourceEllipticPde2_class(py::module &m){
py::class_<CellwiseSourceEllipticPde2 , CellwiseSourceEllipticPde2_Overrides , boost::shared_ptr<CellwiseSourceEllipticPde2 >  , AbstractLinearEllipticPde<2, 2>  >(m, "CellwiseSourceEllipticPde2")
        .def(py::init<::AbstractCellPopulation<2> &, double >(), py::arg("rCellPopulation"), py::arg("sourceCoefficient") = 0.)
        .def(
            "rGetCellPopulation",
            (::AbstractCellPopulation<2> const &(CellwiseSourceEllipticPde2::*)() const ) &CellwiseSourceEllipticPde2::rGetCellPopulation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetCoefficient",
            (double(CellwiseSourceEllipticPde2::*)() const ) &CellwiseSourceEllipticPde2::GetCoefficient,
            " "  )
        .def(
            "ComputeConstantInUSourceTerm",
            (double(CellwiseSourceEllipticPde2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &CellwiseSourceEllipticPde2::ComputeConstantInUSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeLinearInUCoeffInSourceTerm",
            (double(CellwiseSourceEllipticPde2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &CellwiseSourceEllipticPde2::ComputeLinearInUCoeffInSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeLinearInUCoeffInSourceTermAtNode",
            (double(CellwiseSourceEllipticPde2::*)(::Node<2> const &)) &CellwiseSourceEllipticPde2::ComputeLinearInUCoeffInSourceTermAtNode,
            " " , py::arg("rNode") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(CellwiseSourceEllipticPde2::*)(::ChastePoint<2> const &)) &CellwiseSourceEllipticPde2::ComputeDiffusionTerm,
            " " , py::arg("rX") )
    ;
}
