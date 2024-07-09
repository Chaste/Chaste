#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCentreBasedDivisionRule.hpp"

#include "AbstractCentreBasedDivisionRule3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCentreBasedDivisionRule<3,3 > AbstractCentreBasedDivisionRule3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>> _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_3_gt__boost_numeric_ublas_c_vector_lt_double_3_gt__gt_;

class AbstractCentreBasedDivisionRule3_3_Overrides : public AbstractCentreBasedDivisionRule3_3{
    public:
    using AbstractCentreBasedDivisionRule3_3::AbstractCentreBasedDivisionRule;
    ::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>> CalculateCellDivisionVector(::CellPtr pParentCell, ::AbstractCentreBasedCellPopulation<3, 3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_3_gt__boost_numeric_ublas_c_vector_lt_double_3_gt__gt_,
            AbstractCentreBasedDivisionRule3_3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }
    void OutputCellCentreBasedDivisionRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCentreBasedDivisionRule3_3,
            OutputCellCentreBasedDivisionRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCentreBasedDivisionRule3_3_class(py::module &m){
py::class_<AbstractCentreBasedDivisionRule3_3 , AbstractCentreBasedDivisionRule3_3_Overrides , boost::shared_ptr<AbstractCentreBasedDivisionRule3_3 >   >(m, "AbstractCentreBasedDivisionRule3_3")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>>(AbstractCentreBasedDivisionRule3_3::*)(::CellPtr, ::AbstractCentreBasedCellPopulation<3, 3> &)) &AbstractCentreBasedDivisionRule3_3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "OutputCellCentreBasedDivisionRuleInfo",
            (void(AbstractCentreBasedDivisionRule3_3::*)(::out_stream &)) &AbstractCentreBasedDivisionRule3_3::OutputCellCentreBasedDivisionRuleInfo,
            " " , py::arg("rParamsFile") )
    ;
}
