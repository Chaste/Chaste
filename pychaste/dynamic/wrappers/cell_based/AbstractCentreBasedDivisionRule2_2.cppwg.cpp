#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCentreBasedDivisionRule.hpp"

#include "AbstractCentreBasedDivisionRule2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCentreBasedDivisionRule<2,2 > AbstractCentreBasedDivisionRule2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::pair<boost::numeric::ublas::c_vector<double, 2>, boost::numeric::ublas::c_vector<double, 2>> _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_2_gt__boost_numeric_ublas_c_vector_lt_double_2_gt__gt_;

class AbstractCentreBasedDivisionRule2_2_Overrides : public AbstractCentreBasedDivisionRule2_2{
    public:
    using AbstractCentreBasedDivisionRule2_2::AbstractCentreBasedDivisionRule;
    ::std::pair<boost::numeric::ublas::c_vector<double, 2>, boost::numeric::ublas::c_vector<double, 2>> CalculateCellDivisionVector(::CellPtr pParentCell, ::AbstractCentreBasedCellPopulation<2, 2> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_2_gt__boost_numeric_ublas_c_vector_lt_double_2_gt__gt_,
            AbstractCentreBasedDivisionRule2_2,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }
    void OutputCellCentreBasedDivisionRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCentreBasedDivisionRule2_2,
            OutputCellCentreBasedDivisionRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCentreBasedDivisionRule2_2_class(py::module &m){
py::class_<AbstractCentreBasedDivisionRule2_2 , AbstractCentreBasedDivisionRule2_2_Overrides , boost::shared_ptr<AbstractCentreBasedDivisionRule2_2 >   >(m, "AbstractCentreBasedDivisionRule2_2")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::std::pair<boost::numeric::ublas::c_vector<double, 2>, boost::numeric::ublas::c_vector<double, 2>>(AbstractCentreBasedDivisionRule2_2::*)(::CellPtr, ::AbstractCentreBasedCellPopulation<2, 2> &)) &AbstractCentreBasedDivisionRule2_2::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "OutputCellCentreBasedDivisionRuleInfo",
            (void(AbstractCentreBasedDivisionRule2_2::*)(::out_stream &)) &AbstractCentreBasedDivisionRule2_2::OutputCellCentreBasedDivisionRuleInfo,
            " " , py::arg("rParamsFile") )
    ;
}
