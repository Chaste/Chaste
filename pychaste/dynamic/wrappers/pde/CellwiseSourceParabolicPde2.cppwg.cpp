#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellwiseSourceParabolicPde.hpp"

#include "CellwiseSourceParabolicPde2.cppwg.hpp"

namespace py = pybind11;
typedef CellwiseSourceParabolicPde<2 > CellwiseSourceParabolicPde2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class CellwiseSourceParabolicPde2_Overrides : public CellwiseSourceParabolicPde2{
    public:
    using CellwiseSourceParabolicPde2::CellwiseSourceParabolicPde;
    double ComputeDuDtCoefficientFunction(::ChastePoint<2> const & rX) override {
        PYBIND11_OVERRIDE(
            double,
            CellwiseSourceParabolicPde2,
            ComputeDuDtCoefficientFunction,
                    rX);
    }
    double ComputeSourceTerm(::ChastePoint<2> const & rX, double u, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            CellwiseSourceParabolicPde2,
            ComputeSourceTerm,
                    rX,
        u,
        pElement);
    }
    double ComputeSourceTermAtNode(::Node<2> const & rNode, double u) override {
        PYBIND11_OVERRIDE(
            double,
            CellwiseSourceParabolicPde2,
            ComputeSourceTermAtNode,
                    rNode,
        u);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            CellwiseSourceParabolicPde2,
            ComputeDiffusionTerm,
                    rX,
        pElement);
    }

};
void register_CellwiseSourceParabolicPde2_class(py::module &m){
py::class_<CellwiseSourceParabolicPde2 , CellwiseSourceParabolicPde2_Overrides , boost::shared_ptr<CellwiseSourceParabolicPde2 >  , AbstractLinearParabolicPde<2>  >(m, "CellwiseSourceParabolicPde2")
        .def(py::init<::AbstractCellPopulation<2> &, double, double, double >(), py::arg("rCellPopulation"), py::arg("duDtCoefficient") = 1., py::arg("diffusionCoefficient") = 1., py::arg("sourceCoefficient") = 0.)
        .def(
            "rGetCellPopulation",
            (::AbstractCellPopulation<2> const &(CellwiseSourceParabolicPde2::*)() const ) &CellwiseSourceParabolicPde2::rGetCellPopulation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(CellwiseSourceParabolicPde2::*)(::ChastePoint<2> const &)) &CellwiseSourceParabolicPde2::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX") )
        .def(
            "ComputeSourceTerm",
            (double(CellwiseSourceParabolicPde2::*)(::ChastePoint<2> const &, double, ::Element<2, 2> *)) &CellwiseSourceParabolicPde2::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("u"), py::arg("pElement") = __null )
        .def(
            "ComputeSourceTermAtNode",
            (double(CellwiseSourceParabolicPde2::*)(::Node<2> const &, double)) &CellwiseSourceParabolicPde2::ComputeSourceTermAtNode,
            " " , py::arg("rNode"), py::arg("u") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(CellwiseSourceParabolicPde2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &CellwiseSourceParabolicPde2::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pElement") = __null )
    ;
}
