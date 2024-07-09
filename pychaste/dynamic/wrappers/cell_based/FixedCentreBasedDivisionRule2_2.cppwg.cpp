#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FixedCentreBasedDivisionRule.hpp"

#include "FixedCentreBasedDivisionRule2_2.cppwg.hpp"

namespace py = pybind11;
typedef FixedCentreBasedDivisionRule<2,2 > FixedCentreBasedDivisionRule2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::pair<boost::numeric::ublas::c_vector<double, 2>, boost::numeric::ublas::c_vector<double, 2>> _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_2_gt__boost_numeric_ublas_c_vector_lt_double_2_gt__gt_;

class FixedCentreBasedDivisionRule2_2_Overrides : public FixedCentreBasedDivisionRule2_2{
    public:
    using FixedCentreBasedDivisionRule2_2::FixedCentreBasedDivisionRule;
    ::std::pair<boost::numeric::ublas::c_vector<double, 2>, boost::numeric::ublas::c_vector<double, 2>> CalculateCellDivisionVector(::CellPtr pParentCell, ::AbstractCentreBasedCellPopulation<2, 2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_2_gt__boost_numeric_ublas_c_vector_lt_double_2_gt__gt_,
            FixedCentreBasedDivisionRule2_2,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_FixedCentreBasedDivisionRule2_2_class(py::module &m){
py::class_<FixedCentreBasedDivisionRule2_2 , FixedCentreBasedDivisionRule2_2_Overrides , boost::shared_ptr<FixedCentreBasedDivisionRule2_2 >  , AbstractCentreBasedDivisionRule<2>  >(m, "FixedCentreBasedDivisionRule2_2")
        .def(py::init<::boost::numeric::ublas::c_vector<double, 2> & >(), py::arg("rDaughterLocation"))
        .def(
            "rGetDaughterLocation",
            (::boost::numeric::ublas::c_vector<double, 2> const &(FixedCentreBasedDivisionRule2_2::*)() const ) &FixedCentreBasedDivisionRule2_2::rGetDaughterLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "CalculateCellDivisionVector",
            (::std::pair<boost::numeric::ublas::c_vector<double, 2>, boost::numeric::ublas::c_vector<double, 2>>(FixedCentreBasedDivisionRule2_2::*)(::CellPtr, ::AbstractCentreBasedCellPopulation<2, 2> &)) &FixedCentreBasedDivisionRule2_2::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
