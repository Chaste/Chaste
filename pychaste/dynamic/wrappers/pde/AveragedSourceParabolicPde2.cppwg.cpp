#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AveragedSourceParabolicPde.hpp"

#include "AveragedSourceParabolicPde2.cppwg.hpp"

namespace py = pybind11;
typedef AveragedSourceParabolicPde<2 > AveragedSourceParabolicPde2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class AveragedSourceParabolicPde2_Overrides : public AveragedSourceParabolicPde2{
    public:
    using AveragedSourceParabolicPde2::AveragedSourceParabolicPde;
    void SetupSourceTerms(::TetrahedralMesh<2, 2> & rCoarseMesh, ::std::map<boost::shared_ptr<Cell>, unsigned int> * pCellPdeElementMap) override {
        PYBIND11_OVERRIDE(
            void,
            AveragedSourceParabolicPde2,
            SetupSourceTerms,
                    rCoarseMesh,
        pCellPdeElementMap);
    }
    double ComputeDuDtCoefficientFunction(::ChastePoint<2> const & rX) override {
        PYBIND11_OVERRIDE(
            double,
            AveragedSourceParabolicPde2,
            ComputeDuDtCoefficientFunction,
                    rX);
    }
    double ComputeSourceTerm(::ChastePoint<2> const & rX, double u, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            AveragedSourceParabolicPde2,
            ComputeSourceTerm,
                    rX,
        u,
        pElement);
    }
    double ComputeSourceTermAtNode(::Node<2> const & rNode, double u) override {
        PYBIND11_OVERRIDE(
            double,
            AveragedSourceParabolicPde2,
            ComputeSourceTermAtNode,
                    rNode,
        u);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            AveragedSourceParabolicPde2,
            ComputeDiffusionTerm,
                    rX,
        pElement);
    }

};
void register_AveragedSourceParabolicPde2_class(py::module &m){
py::class_<AveragedSourceParabolicPde2 , AveragedSourceParabolicPde2_Overrides , boost::shared_ptr<AveragedSourceParabolicPde2 >  , AbstractLinearParabolicPde<2>  >(m, "AveragedSourceParabolicPde2")
        .def(py::init<::AbstractCellPopulation<2> &, double, double, double >(), py::arg("rCellPopulation"), py::arg("duDtCoefficient") = 1., py::arg("diffusionCoefficient") = 1., py::arg("sourceCoefficient") = 0.)
        .def(
            "rGetCellPopulation",
            (::AbstractCellPopulation<2> const &(AveragedSourceParabolicPde2::*)() const ) &AveragedSourceParabolicPde2::rGetCellPopulation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "SetupSourceTerms",
            (void(AveragedSourceParabolicPde2::*)(::TetrahedralMesh<2, 2> &, ::std::map<boost::shared_ptr<Cell>, unsigned int> *)) &AveragedSourceParabolicPde2::SetupSourceTerms,
            " " , py::arg("rCoarseMesh"), py::arg("pCellPdeElementMap") = nullptr )
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(AveragedSourceParabolicPde2::*)(::ChastePoint<2> const &)) &AveragedSourceParabolicPde2::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX") )
        .def(
            "ComputeSourceTerm",
            (double(AveragedSourceParabolicPde2::*)(::ChastePoint<2> const &, double, ::Element<2, 2> *)) &AveragedSourceParabolicPde2::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("u"), py::arg("pElement") = __null )
        .def(
            "ComputeSourceTermAtNode",
            (double(AveragedSourceParabolicPde2::*)(::Node<2> const &, double)) &AveragedSourceParabolicPde2::ComputeSourceTermAtNode,
            " " , py::arg("rNode"), py::arg("u") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(AveragedSourceParabolicPde2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &AveragedSourceParabolicPde2::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pElement") = __null )
        .def(
            "GetUptakeRateForElement",
            (double(AveragedSourceParabolicPde2::*)(unsigned int)) &AveragedSourceParabolicPde2::GetUptakeRateForElement,
            " " , py::arg("elementIndex") )
    ;
}
